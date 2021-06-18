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
