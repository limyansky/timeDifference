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

    parser.add_argument('--output_root',
                        default='sorted',
                        help='The name of the output files.')

    # Extract the arguments from the parser
    args = parser.parse_args()

    # Load the csv file
    data = pd.read_csv(args.csv_file, header=[0])

    # Sort the csv file
    data.sort_values(['P-Value'], inplace=True)

    # Create the name where the output sorted csv file will be saved
    csv_out = args.output_root + '.csv'

    # Save the sorted csv file
    data.to_csv(csv_out, index=False)

    return 0


# If called from the commandline, run this script.
if __name__ == '__main__':

    main()
