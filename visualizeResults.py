# This script can be used to visualize the results of the blind search.

# Handles commandline inputs
import argparse

# Handles csv files
import pandas as pd


# The function that will run when the script is called
def main():
    """
    Takes .csv output file and produces summary output for further analysis.

    Parameters:
        csv_file: The output of a blind search

    """

    # Create the argument parser and add the docstring
    parser = argparse.ArgumentParser(description=__doc__)

    # Add arguments
    parser.add_argument('csv_file',
                        help='The output .csv file of a blind search')

    # Extract the arguments from the parser
    args = parser.parse_args()


    return 0


# If called from the commandline, run this script.
if __name__ == '__main__':

    main()
