# Imports

# Handles commandline inputs
import argparse

# Varous mathetical functions
import numpy as np

# Handles fits files
from astropy.io import fits

# Allows for importing C code
from ctypes import CDLL, c_int, c_double, POINTER, c_float, c_void_p

# The fourier transform code
import pyfftw

# File io
import os

# Saving wisdom files
import pickle

# Save the outputs as a csv file
import csv

from math import floor

# Useful for debugging
# import timeit
# import resource
# Example usage of memory usage debugging. This will print memory usage in
# kb.
# print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

# Allows me to delete C object to free up memory
from copy import deepcopy

# Run the code on multiple cores
import multiprocessing

# Load the C library
timeDifference = CDLL('/home/brent/github/timeDifference/timeDiff.so')

# Pull the relevant functions from the C library
CtimeDiff = timeDifference.timeDifference_fast
cleanup = timeDifference.cleanup

# Varous functions that I've written, but haven't used
# CinPlace  = timeDifference.timeDifference_inPlace
# Clinear = timeDifference.linearFit
# Csort = timeDifference.wrap_qsort

# Setup the simple cleanup function
cleanup.argtypes = [c_void_p]


# The function that will be run when the script is called
def main():
    """
    Run a blind search for a pulsar using Atwood's time differencing technique.

    Parameters:
        FT1_file: The photon file to search for a pulsar.

    Keyword Arguments:
        --window_size:   The photon file to search for a pulsar
        --max_freq:      Maximum frequency to search in Hz
        --min_freq:      Minimum frequency to search in Hz
        --weight_column: The name of the weights column in FT1_file
    """

    # Create the argument parser and add the docstring
    parser = argparse.ArgumentParser(description=__doc__)

    # Add arguments
    parser.add_argument('FT1_file',
                        help='The photon file to search for a pulsar')

    # Add optional arguments
    parser.add_argument('--window_size',
                        nargs='?',
                        default=524288,
                        help='Maximum time difference in seconds')

    parser.add_argument('--max_freq',
                        nargs='?',
                        default=64,
                        help='Maximum frequency to search in Hz')

    parser.add_argument('--min_freq',
                        nargs='?',
                        default=0.5,
                        help='Minimum frequency to search in Hz')

    parser.add_argument('--weight_column',
                        nargs='?',
                        default=None,
                        help='The name of the weights column in FT1_file')

    parser.add_argument('--lower_f1',
                        nargs='?',
                        default=1,
                        help='The lower f1 value to search.')

    parser.add_argument('--upper_f1',
                        nargs='?',
                        default=1,
                        help='The upper f1 value to search.')

    parser.add_argument('--lower_f2',
                        nargs='?',
                        default=None,
                        help='The lower value of f2 to search.')

    parser.add_argument('--upper_f2',
                        nargs='?',
                        default=None,
                        help='The upper value of f2 to search.')

    parser.add_argument('--wisdom_file',
                        nargs='?',
                        default=None,
                        help='Where to save and load pyfftw\'s wisdom file')

    parser.add_argument('--out_file',
                        nargs='?',
                        default=None,
                        help='Write the output to this file')

    parser.add_argument('--n_cores',
                        default=1,
                        type=int,
                        help='The number of processors to use')

    # Extract the arguments from the parser
    args = parser.parse_args()

    # If an output file was specified, create it and add column headers
    if args.out_file is not None:
        initOutFile(args.out_file)

    # I need to "lock" the output file within each process, as I can't have
    # two processes writing to it at the same time.
    lock = multiprocessing.Lock()

    # Load the photon times and weights from the input files
    times, weights = readEvents(args.FT1_file, args.weight_column)

    # Calculate the epoch as the center of the times
    epoch = (np.amax(times) + np.amin(times)) / 2.

    # Setup a basic search grid in f1/f0
    # Calculate the minimum and maximum f1/f0 ratios
    lower_f1_f0 = args.lower_f1 / args.max_freq
    upper_f1_f0 = args.upper_f1 / args.min_freq

    f1_f0_step = GetF1_F0Step(times, args.window_size, args.max_freq)
    f1_f0_list = GetF1_F0List(f1_f0_step, lower_f1_f0, upper_f1_f0)

    # Set up a search grid in F2/F0.

    # If no F2 is specified, set it to zero
    if args.lower_f2 is None and args.upper_f2 is None:
        f2_f0_list = [0]

    # If both an upper/lower F2 are specified, proceed with generating a list
    # of values to scan over.
    elif args.lower_f2 is not None and args.upper_f2 is not None:
        # Setup a basic search grid in f1/f0
        # Calculate the minimum and maximum f1/f0 ratios
        lower_f2_f0 = args.lower_f2 / args.max_freq
        upper_f2_f0 = args.upper_f2 / args.min_freq

        # Set up a search grid in F2/F0
        f2_f0_step = GetF1_F0Step(times, args.window_size, args.max_freq)
        f2_f0_list = GetF1_F0List(f2_f0_step, lower_f2_f0, upper_f2_f0)

    # If only one upper or lower f2 value was specified, then exit with an
    # explanation
    else:
        print('You must either specify BOTH an upper/lower F2, or neither.')
        return 0

    # Begin the search process
    # OverallBest = [0, 0, 1]

    # If there is a wisdom file specified, but it doesn't exist, I need to
    # create it. The easiest way to do this is to run a single core job first,
    # save the wisdom file, then jump into the multiprocessing code.

    # Assume, at first, that we didn't initalize the wisdom file using the
    # first step
    initalized = False

    if args.wisdom_file is not None and not os.path.isfile(args.wisdom_file):

        # Keep track of if we ran this over the first bin or not
        initalized = True

        # Run a single iteration of the scan, saving the wisdom file
        run_scan(times, weights, args.window_size, args.min_freq,
                 args.max_freq, epoch, f1_f0_list[0], None,
                 args.out_file, save_wisdom=args.wisdom_file)

    # Break up the f1_f0 list into smaller lists, one for each processor.

    # Initalize the list of lists
    f1_f0_master = []

    # Calculate the break points to split the process between cores.
    # This value is roughly (e.g. if you have an odd number of p1/p0, but an
    # even number of cores) the number of jobs per core.
    unit_length = floor((len(f1_f0_list) - 1) / args.n_cores)

    # A place to keep track of where we want to slice the p1_p0 list
    break_points = []

    # Check if we need to skip the first bin or not
    # We may have already calculated this and used it to generate the wisdom
    # files.
    if initalized:

        # Start by skipping the first bin
        break_points.append(1)

    # If we don't need to skip the first bin, then don't skip it...
    elif not initalized:

        # Don't skip the first bin
        break_points.append(0)

    # Fill in subsequent break points
    for ii in range(args.n_cores - 1):

        # We already added the starting index, 1, so we jump ahead
        ii = ii + 1

        # Python uses exclusive upper limits, so we need to add 1
        break_points.append(ii * unit_length + 1)

    # make sure it ends at the last bin
    break_points.append(len(f1_f0_list))

    # Create a template that will be passed to starmap.
    # This holds all the arguments that run_scan needs, other than the
    # p1/p0 values we wish for that particular core to scan over. A placeholder
    # "0" is used for those, and will be filled in in the following loop.
    template_input = [times, weights, args.window_size, args.min_freq,
                      args.max_freq, epoch, 0, 0,
                      args.wisdom_file, args.out_file]

    # Now that we have the break points, we run through them to fill in the
    # list of lists we will be iterating through.

    for jj in f2_f0_list:
        for ii in range(len(break_points) - 1):

            # Fill in the appropriate spot in the template
            template_input[6] = f1_f0_list[break_points[ii]:break_points[ii + 1]]

            # Add the f2_f0 value we want to look at
            template_input[7] = jj

            # This is good to double check that jobs are being partitioned
            # correctly.
            # print(template_input[6])

            # If we don't append a deepcopy, all cores will recieve a copy of
            # only
            # the last p1_p0_list slice and thus all to the same job.
            f1_f0_master.append(deepcopy(template_input))

        # Run through the lists in a multiprocessing way.
        # Locks can't be passed to multiprocessing as arguments (because they
        # are
        # non-pickelable), so I have to use the initalizer and initargs keyword
        # arguments to ensure that each process has access to the lock.
        pool = multiprocessing.Pool(initializer=init_lock, initargs=(lock,),
                                    processes=args.n_cores)

        # Actually runs the multiprocessing
        pool.starmap(run_scan, f1_f0_master)

        # closes the pool when we are done
        pool.close()

        # Show a summary of the search
        print("\nScan in -f1/f0 completed after %d steps" % (len(f1_f0_list)))

    return 0


# Runs through a list of p1_p0 values
# load_wisdom is not initalized to None because multiprocessing starmap cannot
# easily pass in keyword arguments. The same is true of out_file.
def run_scan(times, weights,
             window_size, min_freq, max_freq, epoch, f1_f0_list, f2_f0,
             load_wisdom, out_file, save_wisdom=None):
    """
    Runs a single search step in F1_F0

    Parameters:
        times: A list of photon times
        weights: A list of photon weights
        window_size: The maximum time difference, in seconds
        max_freq: The maximum frequency to search
        epoch: The epoch from which to correct photon times
        f1_f0_lsit: A list of F1/F0 steps over which to search
        f2_f2: A single F2/F0 value to examine

    Keyword Arguments:
        save_wisdom: A place to save the fft wisdom file
        load_wisdom: The location of a saved fft wisdom file
        out_file: A place to save outputs in .csv format

    """

    # Load the wisdom file, if requested
    if load_wisdom is not None:
        LoadWisdom(load_wisdom)

    # Initalize pyfftw
    fftw_object, fft_input, fft_output = init_FFTW(window_size, max_freq)

    step = 0

    for f1_f0 in f1_f0_list:

        step += 1
        print(multiprocessing.current_process(),
              'Step: ', step, '/', len(f1_f0_list))

        # Correct the times
        new_times = TimeWarp_F1_F2(times, f1_f0, f2_f0, epoch)

        # Find the time differences
        time_differences = call_CtimeDiff(CtimeDiff,
                                          new_times,
                                          weights,
                                          windowSize=window_size,
                                          maxFreq=max_freq)

        # run the FFTW
        power_spectrum = run_FFTW(fftw_object, fft_input, fft_output,
                                  time_differences)

        # If needed, save the wisdom file
        if save_wisdom is not None:

            # Save the wisdom file
            SaveWisdom(save_wisdom)

            # Make sure we don't continue to save the wisdom file
            save_wisdom = None

        # Extract the best candidate
        [freq, p_value] = ExtractBestCandidate(power_spectrum,
                                               min_freq, max_freq)

        # Acquire the lock
        lock.acquire()

        # Print the candidate
        DisplayCandidate([freq, f1_f0, f2_f0, p_value], out_file=out_file)

        # Release the lock
        lock.release()


# This us used in initalization of multiprocessing so that each child
# process has access to the file lock preventing them from simultaneosuly
# appending to the csv output file.
def init_lock(lock_handle):
    global lock
    lock = lock_handle


# Calculates the size of the FFT
def FFT_Size(window_size, max_freq):
    """
    Calculate the size of the FFT

    Parameters:
        window_size: Maximum time difference in seconds
        max_freq:    Maximum frequency to search in Hz
    """
    return 2 * int(np.floor(window_size * max_freq))


# Reads events from an FT1 File
# Returns sorted events
def readEvents(FT1_file, weight_column=None):
    """
    Read events from an FT1 file

    Parameters:
        FT1_file: The photon file to search for a pulsar

    Keyword Arguments:
        --weight_column: The name of the weights column in FT1_file
    """

    # Open the fits file
    hdu = fits.open(FT1_file)

    # Extract the times
    times = hdu[1].data['TIME']

    # An array that will sort times and weights
    sort_indicies = np.argsort(times)

    # Sort the times
    times = times[sort_indicies]

    # If a weight column is supplied, extract them
    if weight_column is not None:
        weights = hdu[1].data[weight_column]

        # Sort the weights
        weights = weights[sort_indicies]
    else:
        weights = 1 * len(times)

    # Close the file
    hdu.close()

    return times, weights


# A function to call the C code which performs the time differencing
def call_CtimeDiff(function, photons, weights, windowSize=524288, maxFreq=64):
    """
    Call a C function to quickly calculate time differences.

    Parameters:
        function: The handler from ctypes to work with
        photons: A list of photon times to difference
        weights: A list of photon weights corresponding to the times

    Keyword Arguments:
        WindowSize: THe maximum time difference to include
        maxFreq: The maximum frequency to search
    """

    # Specify the datatypes that will be used as inputs
    function.argtypes = [POINTER(c_double), POINTER(c_double),
                         c_int, c_int, c_int]

    # Specify the datatypes that will be used as outputs
    function.restype = POINTER(c_double * FFT_Size(windowSize, maxFreq))

    # The photons and weights need to be converted into something C can read.
    cPhotons = (c_double * len(photons))(*photons)
    cWeights = (c_double * len(weights))(*weights)

    # Calculate the time differences
    histogram = function(cPhotons, cWeights, windowSize, maxFreq,
                         len(cPhotons))

    # The C output needs to be converted back into something python can read.
    # Not creating a deepcopy will cause a memory leak
    output = deepcopy(np.frombuffer(histogram.contents, dtype=np.double))

    # C function, initalized at the start of this file, that calls "free" on
    # the input. This fixes a memory leak where if histogram.contents is
    # extracted, but never freed, it stays in memory.
    cleanup(histogram)

    return output


# A function to call the C code which performs the time differencing
def call_CinPlace(function, photons, weights, windowSize=524288, maxFreq=64):
    """
    Call a C function to quickly calculate time differences.

    Parameters:
        function: The handler from ctypes to work with
        photons: A list of photon times to difference
        weights: A list of photon weights corresponding to the times

    Keyword Arguments:
        WindowSize: The maximum time difference to include
        maxFreq: The maximum frequency to search
    """

    # Specify the datatypes that will be used as inputs
    function.argtypes = [POINTER(c_double), POINTER(c_double),
                         POINTER(c_double),
                         c_int, c_int, c_int]

    histogram = np.zeros(FFT_Size(windowSize, maxFreq))

    # The photons and weights need to be converted into something C can read.
    cPhotons = (c_double * len(photons))(*photons)
    cWeights = (c_double * len(weights))(*weights)
    cHistogram = (c_double * len(histogram))(*histogram)

    # Calculate the time differences
    function(cPhotons, cWeights, cHistogram,
             windowSize, maxFreq, len(cPhotons))

    print(sum(histogram))

    return histogram


def call_Csort(function, toSort):
    """
    Description

    Parameters:
    """

    # Define the inputs and outputs
    function.argtypes = [POINTER(c_float), c_int]
    function.restype = POINTER(c_float)

    # Convert the input into something C can read
    CtoSort = (c_float * len(toSort))(*toSort)

    # Call the funciton
    toSort = function(CtoSort, len(toSort))


# Get the step in F1/F0
def GetF1_F0Step(times, windowSize=524288, maxFreq=64):
    """
    Determine the grid step in the P1/P0 parameter space such that the maximum
    tolerated frequency drift over the full time span covered by the data
    is smaller than the FFT resolution (see eq.3 in Ziegler et all 2008) for
    the largest frequency considered int the search.

    """

    # Find the span covered by the data.
    time_span = np.amax(times) - np.amin(times)

    # The fft resolution
    FFT_resol = 1. / windowSize

    # The F1 resolution
    f1_tolerance = 1. * FFT_resol / time_span

    # At least one point in the grid is within 1/2 the grid step from the
    # correct value of the parameter (f1/f0)
    return 2. * f1_tolerance / maxFreq


# Generate a list of F1/F0 values to scan over
def GetF1_F0List(f1_f0_step, lower_f1_f0=-1.3e-11, upper_f1_f0=0.):
    """
    Generate a list of F1/F0 values to scan over

    Parameters:
        f1_f0_step: The step size of F1/F0

    Keyword Arguments:
        lower_f1_f0: The lower limit for the list values you want
        upper_f1_f0: The upper limit for the list values you want
    """

    # This is a one-liner
    # Add one last step, even if it falls outside of the range
    return np.arange(lower_f1_f0, upper_f1_f0 + f1_f0_step, f1_f0_step)


# Get the step in F2/F0
def GetF2_F0Step(times, windowSize=524288, maxFreq=64):

    # Find the span covered by the data.
    time_span = np.amax(times) - np.amin(times)

    # The fft resolution
    FFT_resol = 1. / windowSize

    # The F2 resolution
    f2_tolerance = 2. * FFT_resol / time_span**2

    # At least one point in the grid is within 1/2 the grid step from the
    # correct value of the parameter (f1/f0)
    return 2. * f2_tolerance / maxFreq


# Generate a list of F1/F0 values to scan over
def GetF2_F0List(f2_f0_step, lower_f2_f0=-1.3e-11, upper_f2_f0=0.):
    """
    Generate a list of F2/F0 values to scan over

    Parameters:
        f2_f0_step: The step size of F2/F0

    Keyword Arguments:
        lower_f2_f0: The lower limit for the list values you want
        upper_f2_f0: The upper limit for the list values you want
    """

    # This is a one-liner
    # Add one last step, even if it falls outside of the range
    return np.arange(lower_f2_f0, upper_f2_f0 + f2_f0_step, f2_f0_step)


# Find the appropriate step in P1/P0
def GetP1_P0Step(times, windowSize=524288, maxFreq=64):
    """
    Determine the grid step in the P1/P0 parameter space such that the maximum
    tolerated frequency drift over the full time span covered by the data
    is smaller than the FFT resolution (see eq.3 in Ziegler et all 2008) for
    the largest frequency considered int the search.

    Parameters:
        times: A list of photon times

    Keyword Arguments:
        windowSize: The maximum time difference
        maxFreq: The maximum frequency to search
    """

    # Find the span covered by the data.
    time_span = np.amax(times) - np.amin(times)

    # The fft resolution
    FFT_resol = 1. / windowSize

    # This value is somewhat arbitrary: a finer grid provides
    # a better sensitivity, but is more time-consuming; a coarser
    # grid is clearly faster but you risk of missing the pulsar.
    f1_tolerance = 1. * FFT_resol / time_span

    # at least one point in the grid is within 1/2 the grid
    # step from the correct value of the parameter (p1/p0)
    return 2. * f1_tolerance / maxFreq


# Generate a list of P1/P0 values to scan over
def GetP1_P0List(p1_p0_step, lower_p1_p0=0., upper_p1_p0=1.3e-11):
    """
    Generate a list of P1/P0 values to scan over

    Parameters:
        p1_p0_step: The step size of P1/P0

    Keyword Arguments:
        lower_p1_p0: The lower limit for the list values you want
        pper_p1_p1: The upper limit for the list values you want
    """

    # This is a one-liner
    # Add one last step, even if it falls outside of the range
    return np.arange(lower_p1_p0, upper_p1_p0 + p1_p0_step, p1_p0_step)


# Correct times for another step in P1/P0
def TimeWarp(times, p1_p0, epoch):
    """
    Compensate for a steady frequency drift. This way, the new time series is
    periodic and can be searched with standard FFTs. The time transform is
    described in Ransom et al. 2001. See also eq 1 in Ziegler et al. 2008.

    Parameters:
        times: A list (numpy array) of photon times
        p1_p0: The P1/P0 value to correct
        epoch: The epoch of the timing solution
    """

    # the minus sign enters because p1_p0 is p1/p0=-f1/f0
    times = times - 0.5 * p1_p0 * (times - epoch)**2

    return times


# Correct times for a step in F1/F0
def TimeWarp_F1(times, f1_f0, epoch):
    """
    Compenstate for a steady frequency drift. This way the new time series is
    periodic and can be searched with standard FFTs. The time transform is
    described in Ransom et al. 2001. See also eq 1 in Ziegler et al. 2008.

    Parameters:
        times: A list (numpy array) of photon times
        f1_f0: The F1/F0 value to correct
        epoch: The epoch of the timing solution
    """

    times = times + 0.5 * f1_f0 * (times - epoch) ** 2

    return times


# Correct times for a step in F2/F0
def TimeWarp_F2(times, f2_f0, epoch):
    """
    Compenstate for a steady frequency drift acceleration.
    This way the new time series is periodic and can be searched with standard
    FFTs.

    Paramters:
        times: A list of (numpy array) photon times
        f2_f0: The F2/F0 value to correct
        epoch: the Epoch of the timing solution
    """

    times = times + (1 / 6) * f2_f0 * (times - epoch) ** 3

    return times


# Correct times for a step in both F1 and F2
def TimeWarp_F1_F2(times, f1_f0, f2_f0, epoch):
    """
    Compensate for a F1 and F2 frequency drift.
    This way the new time series is periodic and can be searched with standard
    FFTs.

    Parameters:
        times: A list (numpy array) of photon times
        f1_f0: The F1/F0 value to correct
        f2_f0: The F2/F0 value to correct
        epoch: The epoch of the timing solution
    """

    times = times + (0.5 * f1_f0 * (times - epoch) ** 2) + \
                    ((1 / 6) * f2_f0 * (times - epoch) ** 3)

    return times


# Load a wisdom file for pyfftw
def LoadWisdom(wisdom_file):
    """
    Loads a wisdom file for use in pyfftw.

    Parameters:
        wisdom_file: The location of the wisdom file to load
    """

    with open(wisdom_file, 'rb') as f:
        wisdom = pickle.load(f)
        pyfftw.import_wisdom(wisdom)


# Saves a wisdom file for pyfftw
def SaveWisdom(wisdom_file):
    """
    Saves a wisdom file for pyfftw

    Parameters:
        wisdom_file: A wisdom file to save after generated by pyfftw
    """

    with open(wisdom_file, 'wb') as f:
        pickle.dump(pyfftw.export_wisdom(), f)
        print('Saved pickled wisdom.')

# Initalize the inputs and outputs of the fft
def init_FFTW(window_size, max_freq):

    # Calculate the size of the fft
    FFT_size = FFT_Size(window_size, max_freq)
    alignment = pyfftw.simd_alignment

    # this is tricky: it is needed to get the correct memory alignment for fftw
    input_array = pyfftw.empty_aligned(FFT_size,
                                       n=alignment,
                                       dtype='float32')

    output_array = pyfftw.empty_aligned(FFT_size // 2 + 1,
                                        n=alignment,
                                        dtype='complex64')

    # create the FFT object, BEFORE actually loading the data!!!!
    fft_object = pyfftw.FFTW(input_array, output_array, threads=1)

    # Return the object, input array, and output array
    return fft_object, input_array, output_array


def run_FFTW(fft_object, input_array, output_array, time_differences):
    # load the actual input into the allocated memory
    input_array[:] = time_differences

    # this normalization grants that, if the input array is Poisson
    # distributed,
    # the Fourier power follows a chi2 distribution with 2 degrees of freedom
    # unfortunately the time differences are NOT Poisson distributed...
    norm = np.sum(np.absolute(input_array) / 2.0, dtype=np.float32)

    # FFTW.__Call__ automatically executes the FFT and returns the output array
    # output_array = fft_object()
    fft_object()

    # return the normalized Fourier power
    return np.square(np.absolute(output_array)) / norm


# Run fftw
def FFTW_Transform(time_differences, window_size, max_freq):
    """
    Runs fftw on the data.

    Parameters:
        time_differences: The array of photon time differences
        window_size: The maximum time difference
        max_freq: The maximum frequency to scan
    """

    # Calculate the size of the fft
    FFT_size = FFT_Size(window_size, max_freq)
    alignment = pyfftw.simd_alignment

    # this is tricky: it is needed to get the correct memory alignment for fftw
    input_array = pyfftw.n_byte_align_empty(FFT_size,
                                            alignment,
                                            dtype='float32')
    output_array = pyfftw.n_byte_align_empty(FFT_size // 2 + 1,
                                             alignment,
                                             dtype='complex64')

    # create the FFT object, BEFORE actually loading the data!!!!
    fft_object = pyfftw.FFTW(input_array, output_array, threads=1)

    # load the actual input into the allocated memory
    input_array[:] = time_differences

    # this normalization grants that, if the input array is Poisson
    # distributed,
    # the Fourier power follows a chi2 distribution with 2 degrees of freedom
    # unfortunately the time differences are NOT Poisson distributed...
    norm = np.sum(np.absolute(input_array) / 2.0, dtype=np.float32)

    # FFTW.__Call__ automatically executes the FFT and returns the output array
    fft_object()

    # return the normalized Fourier power
    return np.square(np.absolute(output_array)) / norm


def ExtractBestCandidate(power_spectrum, min_freq, max_freq):
    """
    Pick the candidate frequency with the largest Fourier power, within the
    prescribed frequency range. Return the par of its frequency and the P-value
    of observing such a power as a statistical fluctuation. For a discussion,
    see Sec 4. of Ziegler et al. 2008.

    Parameters:
        power_spectrum: The output spectrum of pyfftw
        min_freq: The lowest frequency to search
        max_freq: The largest frequency to search

    """

    # This value is roughly correct, though FFT_resol := 1/window_size
    FFT_resol = float(max_freq) / (len(power_spectrum) - 1.0)

    # Ignore peaks at the lowest frequencies, in order to avoid red noise
    min_index = int(np.floor(float(min_freq) / FFT_resol))
    peak_index = min_index + np.argmax(power_spectrum[min_index:])

    power_spectrum = power_spectrum[min_index:]

    # We need this operation of CPU complexity NlogN to interpret the power
    power_spectrum = np.sort(power_spectrum[min_index:], kind='heapsort')[::-1]

    # start = timeit.default_timer()
    # call_Csort(Csort, power_spectrum)
    # power_spectrum = power_spectrum[::-1]
    # end = timeit.default_timer()
    # print('qsort done in ', end - start, ' seconds')

    # Work out the asymptotic cumulative distribution of the power values
    [slope, constant] = FitExponentialTail(power_spectrum)

    # Apply the asymptotic results to convert the power into a P-value
    P_value = PowerToPValue(power_spectrum[0], slope, constant)

    return [FFT_resol * peak_index, P_value]


def FitExponentialTail(sorted_array):
    """
    Analyze the probability distribution of values in an array and fit the tail
    with an exponential function. THe array is assumed to be already sorted
    and in decreasing order.

    Parameters:
        sorted_array: A sorted power spectrum from pyfftw
    """

    # We define the tail through an emprical approximation
    if len(sorted_array) > 2000000:
        start_index = 200
        end_index = 20000
    else:
        start_index = int(len(sorted_array) / 10000)
        end_index = int(len(sorted_array) / 100)

    # consider only the fit range and assume an exponential shape
    X_values = sorted_array[start_index:end_index]
    Y_values = np.log(np.arange(start_index + 1, end_index + 1))

    # Find the correlation between x and y
    correlation = np.corrcoef(X_values, Y_values)[0, 1]
    std_X = np.std(X_values)
    std_Y = np.std(Y_values)

    # Find the slope
    slope = correlation * (std_Y / std_X)

    intercept = np.mean(Y_values) - (slope * np.mean(X_values))

    # The below two lines were in the original code:
    # this is a bit overkilling: a line is a polynomial of order 1
    # return np.polyfit(X_values, Y_values, 1)

    # We are  back to my code/method of doing this.
    return [slope, intercept]


def PowerToPValue(power, slope, constant):
    """
    Apply the asymptotic cumulative distribution of the power under the null
    hypothesis of no pulsations to estimate the probability of getting the
    observed values due to chance (P-Value).

    Parameters:
        Power:
        Slope: The slope of the exponential tail fit
        Constant:
    """

    # this value accounts already for the trials due to FFT bins
    # it comes from an empirical fit, and is very approximative
    effective_power = power - np.sqrt(power)
    return np.min([np.exp(constant + slope * effective_power), 1])


def initOutFile(outFile):
    """
    Initalizes the output .csv file by writing column headers.

    Parameters:
        outFile: The name of the output file.
    """

    # The row that we want to write to the file
    row = ['F0', 'F1', 'P-Value']

    # Open the file for writing
    with open(outFile, 'w') as f:

        # Create the writer object
        writer = csv.writer(f)

        # Write the column headers to the file
        writer.writerow(row)


def DisplayCandidate(candidate, best=False, out_file=None):
    """
    Print the basic information about a pulsar candidate.

    Parameters:
        candidate: A candidate to print parameters about

    Keyword arguments:
        best: Prints some additional information about the candidate
        outFile: Saves the output to a .csv file
    """

    if best:
        print("\nThe best pulsar candidate is:")
    # the second entry in the candidate is the value of p1/p0=-f1/f0
    Fdot = -1. * candidate[1] * candidate[0]
    # print("F0=%.8f F1=%.3e P-Value=%.2e" % (candidate[0], Fdot,
    #       candidate[2]))

    # If outFile has been provided, save a CSV
    if out_file is not None:

        with open(out_file, 'a') as f:

            # Create the writer object.
            writer = csv.writer(f)

            # Write the row
            writer.writerow([candidate[0], Fdot, candidate[2]])

    if best:
        if candidate[1] == 0:
            print("Characteristic age and Edot not available (F1 null)")
        else:
            print("Characteristic age=%.2e years" %
                  (3.1688e-08 / candidate[1]))
            print("Edot=I45*%.2e erg/s\n" % (-3.9478e46 * Fdot * candidate[0]))


# If called from the commandline, run this script.
if __name__ == '__main__':

    # Just try deleting the below line. I dare you.
    multiprocessing.set_start_method("spawn")

    main()
