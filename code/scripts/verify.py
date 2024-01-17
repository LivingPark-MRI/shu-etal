import os, glob
import pandas as pd

def check_data(dirr):
    """
    Check the existence and structure of required data files in a specified directory.

    Parameters:
    - dirr (str): The directory path where the data files are expected to be located.

    Returns:
    - bool: True if the required data files are found and have the expected structure, False otherwise.
    
    The function checks for the presence of two files: 'cohort.csv' and 'demographics*.csv' in the specified directory.
    It further verifies that the columns in these files match the expected column names. The expected columns are:
    For 'cohort.csv': ["PATNO", "EVENT_ID", "Description"]
    For 'demographics*.csv' (first file matching the pattern): ["PATNO", "Age", "Stage", "SEX", "output"]
    """
    cohort = f"{dirr}/cohort.csv"
    demographics = glob.glob(f"{dirr}/demographics*.csv")

    # Check if demographics files exist
    if len(demographics) == 0:
        return False

    # Check if cohort and demographics files exist
    if not os.path.exists(cohort) or not os.path.exists(demographics[0]):
        return False

    # Check that columns are in the data
    if not pd.read_csv(cohort).columns.isin(["PATNO", "EVENT_ID", "Description"]).all() or \
       not pd.read_csv(demographics[0]).columns.isin(["PATNO", "Age", "Stage", "SEX", "output"]).any():
        return False

    return True
