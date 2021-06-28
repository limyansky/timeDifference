# Imports

# Handles commandline inputs
import argparse

# Varous mathetical functions
import numpy as np

# Handles fits files
from astropy.io import fits

# Allows for importing C code
from ctypes import CDLL, c_int, c_double, POINTER, c_float

# The fourier transform code
import pyfftw

# File io
import os

# Savig wisdom files
import pickle

import timeit

import tracemalloc

import gc

# Prepare the C code
timeDifference = CDLL('/home/brent/github/timeDifference/timeDiff.so')
CtimeDiff = timeDifference.timeDifference_fast
CinPlace  = timeDifference.timeDifference_inPlace
#Clinear = timeDifference.linearFit
Csort = timeDifference.wrap_qsort


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

    parser.add_argument('--lower_p1_p0',
                        nargs='?',
                        default=0,
                        help='The lower bound of p1/p0 in s^-1')

    parser.add_argument('--upper_p1_p0',
                        nargs='?',
                        default=1.3e-11,
                        help='The upper bound of p1/p0 in s^-1')

    parser.add_argument('--wisdom_file',
                        nargs='?',
                        default=None,
                        help='Where to save and load pyfftw\'s wisdom file')

    # Extract the arguments from the parser
    args = parser.parse_args()

    # Read in the times and weights
    times, weights = readEvents(args.FT1_file, args.weight_column)

    # Calculate the epoch as the center of the times
    epoch = (np.amax(times) + np.amin(times)) / 2.

    # Setup a basic search grid in -f1/f0
    p1_p0_step = GetP1_P0Step(times, args.window_size, args.max_freq)
    p1_p0_list = GetP1_P0List(p1_p0_step, args.lower_p1_p0, args.upper_p1_p0)

    # Begin the search process
    # Keep track of how many steps we have done
    step = 0
    OverallBest = [0, 0, 1]

    tracemalloc.start()

    save_wisdom = False
    load_wisdom = False

    # Check for the existance of a wisdom file
    if args.wisdom_file is not None and os.path.isfile(args.wisdom_file):
        load_wisdom = True
    elif args.wisdom_file is not None and not os.path.isfile(args.wisdom_file):
        save_wisdom = True


    # Load the pyfftw wisdom file
    if load_wisdom:
        LoadWisdom(args.wisdom_file)
        load_wisdom = False

    # Initalize the pyFFTW object
    fftw_object, input_array, output_array = init_FFTW(args.window_size,
                                                       args.max_freq)


    # Step through the list of p1_p0
    for p1_p0 in p1_p0_list:

        snap_start = tracemalloc.take_snapshot()

        print('step: ', step)

        # Update the step number
        step += 1

        # Correct the times
        new_times = TimeWarp(times, p1_p0, epoch)

        print("Calculating time differences.")
        # Find the time differences
        # time_differences = call_CtimeDiff(CtimeDiff,
        #                                   new_times,
        #                                   weights,
        #                                   windowSize=args.window_size,
        #                                   maxFreq=args.max_freq)

        time_differences = call_CinPlace(CinPlace,
                                         new_times,
                                         weights,
                                         windowSize=args.window_size,
                                         maxFreq=args.max_freq)

        # # Load the pyfftw wisdom file
        # if load_wisdom:
        #     LoadWisdom(args.wisdom_file)
        #     load_wisdom = False

        print("Calculating the power spectrum")
        # run the FFTW
        # power_spectrum = FFTW_Transform(time_differences,
        #                                 args.window_size,
        #                                 args.max_freq)
        power_spectrum = run_FFTW(fftw_object, input_array, output_array,
                                  time_differences)


        if save_wisdom:
            SaveWisdom(args.wisdom_file)
            save_wisdom = False

        print("Extracting candidates")
        # Extract the best candidate
        [freq, p_value] = ExtractBestCandidate(power_spectrum,
                                               args.min_freq,
                                               args.max_freq)
        # Print the candidate
        DisplayCandidate([freq, p1_p0, p_value])

        #print(tracemalloc.take_snapshot())
        # If the candidate is the overall best, store it
        if p_value <= OverallBest[2]:
            OverallBest = [freq, p1_p0, p_value]

        snap_end = tracemalloc.take_snapshot()

        stats = snap_end.compare_to(snap_start, 'lineno')

        for ii in stats[:5]:
            print(ii)


    # Show a summary of the search
    print("\nScan in -f1/f0 completed after %d steps" % (len(p1_p0_list)))
    DisplayCandidate(OverallBest, best=True)


    return 0


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
    histogram = np.frombuffer(histogram.contents)

    return histogram


# A function to call the C code which performs the time differencing
def call_CinPlace(function, photons, weights, windowSize=524288, maxFreq=64):
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

    # The C output needs to be converted back into something python can read.
    histogram = np.frombuffer(cHistogram.contents)

    return histogram


# A function to call the C code whic performa a linear fit
def call_Clinear(function, sortedPower):
    """
    Description

    Parameters:
    """

    print('Specifying arg and res types')
    # Specify the datatypes that will be sued as inputs
    function.argtypes = [POINTER(c_double), c_int]
    function.restype = POINTER(c_double * 2)

    print('converting the array into a c double pointer')
    # Convert the input into something C can read
    CsortedPower = (c_double * len(sortedPower))(*sortedPower)

    print('calling the function.')
    # Run the function
    output = function(CsortedPower, len(sortedPower))

    # Print the output
    print(output)


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
    #output_array = fft_object()
    fft_object()

    #print("Deleting input")
    #del input_array

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
    #output_array = fft_object()
    fft_object()

    #print("Deleting input")
    #del input_array

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
    print(type(power_spectrum[0]))

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
    print('start_fit')

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
    slope = correlation * (std_Y/std_X)

    intercept = np.mean(Y_values) - (slope * np.mean(X_values))

    #print('calling call_Clinear')
    #call_Clinear(Clinear, X_values)
    #print('finished call_Clinear')

    # this is a bit overkilling: a line is a polynomial of order 1
    #return np.polyfit(X_values, Y_values, 1)
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
    print(np.exp(constant + slope * effective_power))
    return np.min([np.exp(constant + slope * effective_power), 1])


def DisplayCandidate(candidate, best=False):
    """
    Print the basic information about a pulsar candidate.

    Parameters:
        candidate: A candidate to print parameters about

    Keyword arguments:
        best: Prints some additional information about the candidate
    """
    if best:
        print("\nThe best pulsar candidate is:")
    # the second entry in the candidate is the value of p1/p0=-f1/f0
    Fdot = -1. * candidate[1] * candidate[0]
    print("F0=%.8f F1=%.3e P-Value=%.2e" % (candidate[0], Fdot, candidate[2]))
    if best:
        if candidate[1] == 0:
            print("Characteristic age and Edot not available (F1 null)")
        else:
            print("Characteristic age=%.2e years" %
                  (3.1688e-08 / candidate[1]))
            print("Edot=I45*%.2e erg/s\n" % (-3.9478e46 * Fdot * candidate[0]))


# If called from the commandline, run this script.
if __name__ == '__main__':
    main()
