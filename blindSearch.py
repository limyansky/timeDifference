# Imports

# Handles commandline inputs
import argparse

# Varous mathetical functions
import numpy as np

# Handles fits files
from astropy.io import fits


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

    readEvents(args.FT1_file, args.weight_column)

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

# Calculates the time binning
def timeDiffs(times, weights, window_size=524288, max_freq=64):
    """
    Calculates the time difference between events and bins the results.

    Parameters:
        times: A list of photon times
        weights: A list of photon weights

    Keyword Arguments:
        window_size:
        max_freq:
    """

    time_resol = 0.5 / max_freq

    return 0


# If called from the commandline, run this script.
if __name__ == '__main__':
    main()
