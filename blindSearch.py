# Imports

# Handles commandline inputs
import argparse


# The function that will be run when the script is called
def main():
    """
    Run a blind search for a pulsar using Atwood's time differencing technique.

    Parameters:
        param_1: The first parameter
        param_2: The second parameter

    Keyword Arguments:
        opt_1: The first optional parameter
        opt_2: The second optional parameter
    """

    # Create the argument parser and add the docstring
    parser = argparse.ArgumentParser(description=__doc__)

    # Add Arguments
    parser.add_argument('param_1', help='Soem description')

    # Extract the arguments from the parser
    args = parser.parse_args()

    return 0


# If called from the commandline, run this script.
if __name__ == '__main__':
    main()
