# Imports

# Handles commandline inputs
import argparse

# Varous mathetical functions
import numpy as np

# Handles fits files
from astropy.io import fits

# Allows for importing C code
from ctypes import CDLL, c_int, c_double, POINTER

# Prepare the C code
timeDifference = CDLL('/home/brent/github/timeDifference/timeDiff.so')
CtimeDiff = timeDifference.timeDifference_fast


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

    # Extract the arguments from the parser
    args = parser.parse_args()

    times, weights = readEvents(args.FT1_file, args.weight_column)

    timeDiffs = call_CtimeDiff(CtimeDiff, times, weights)
    print(timeDiffs)

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
def call_CtimeDiff(function, photons, weights, windowSize=524288, maxFreq=68):
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
  """TimeWarp()"""
  """Stretch a time series in order to compensate for a"""
  """steady frequency drift. This way, the new time series"""
  """is periodic and can be searched with standard FFTs."""
  """This time transform is described in Ransom et al. 2001"""
  """See also eq 1 in Ziegler et al. 2008"""
  # the minus sign enters because p1_p0 is p1/p0=-f1/f0

    times = times - 0.5 * p1_p0 * (times - epoch)**2
    return times



def TimeDiffs(times, weights, window_size=524288, max_freq=64):
    """TimeDiffs()"""
    """Extract the binned series of time differences"""
    """The argument max_freq determines the bin size"""
    """as time_resol = 1 / (2 * max_freq) """
    """This together with window_size fixes the size"""
    """of the returned array of time differences"""
    # FFT sampling time
    time_resol = .5 / max_freq
    # directly bin the time differences
    time_diffs = [0] * FFT_Size(window_size, max_freq)
    for i1 in range(len(times) - 1):
        t1 = times[i1]
        for i2 in range(i1 + 1, len(times)):
            t2 = times[i2]
            # limit the size of the time differences
            if t2 - t1 >= window_size:
                break
            # determine the frequency bin in the array
            freq_bin = int(np.floor((t2 - t1) / time_resol))
            # combine the weights appropriately
            time_diffs[freq_bin] += weights[i1] * weights[i2]
    return time_diffs

# If called from the commandline, run this script.
if __name__ == '__main__':
    main()
