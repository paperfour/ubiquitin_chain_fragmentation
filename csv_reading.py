# CSV file reading

import pandas as pd


SINGLE_LOOKUP_COLUMN = []
COMPOSITION_LOOKUP_COLUMN = []
MASS_LOOKUP_COLUMN = []
AMINO_ACID_DF = None

def loadRefSheet():

    global AMINO_ACID_DF

    try:
        print("Reading local csv...")
        AMINO_ACID_DF = pd.read_csv("amino_acid_lookup_info.csv")
    except:
        raise Exception("No local lookup found!")
        # URL should be combined data from https://.washington.edu/protocols06/masses.php and https://raw.githubusercontent.com/specht/proteomics-knowledge-base/master/amino-acids.csv
        url = ""
        AMINO_ACID_DF = pd.read_csv(url)
        AMINO_ACID_DF.to_csv("amino_acid_lookup.csv", index = False)
        print("File downloaded and ready for reference")

    global SINGLE_LOOKUP_COLUMN
    SINGLE_LOOKUP_COLUMN = AMINO_ACID_DF.loc[:, "single letter code"].tolist()
    
    global COMPOSITION_LOOKUP_COLUMN
    COMPOSITION_LOOKUP_COLUMN = AMINO_ACID_DF.loc[:, "composition"].tolist()

    global MASS_LOOKUP_COLUMN
    MASS_LOOKUP_COLUMN = AMINO_ACID_DF.loc[:, "average mass"].tolist()

    print("Imports finished!")

#--------------------------------------------------------------------------

import os

# Reads the last line of a csv file then removes said line from the original file
def popLastLine(file_path):
    
    # Open the file in read and write mode
    with open(file_path, 'rb+') as file:

        # Return None if the file is empty
        file.seek(0) 
        if not file.read(1):
          print("EMPTY FOUND")
          return None


        file.seek(-2, os.SEEK_END)  # Move to the second last byte
        

        while file.read(1) != b'\n':  # Until EOL is found
            try:    
                file.seek(-2, os.SEEK_CUR)  # Move backward by one byte
                # 1 byte before the EOL
            except:
                # The beginning of the file has been reached... we must truncate to 0
                file.seek(0)
                lastLine = file.readline().decode()
                file.truncate(0)  # Truncate the file, it should be empty now
                return lastLine
        
        
        end = file.tell() - 1 
        lastLine = file.readline().decode() 

        file.truncate(end)  # Truncate the file

    return lastLine

'''
import time

n = popLastLine("output/OLD_OPTIMIZED_07_08_2024_12_03_45_out.csv")

t = time.time()
for i in range(0, 2000):
    n = popLastLine("output/OLD_OPTIMIZED_07_08_2024_12_03_45_out.csv")

print(time.time() - t)
''' 


