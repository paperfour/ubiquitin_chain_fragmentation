# CSV file reading

import pandas as pd


SINGLE_LOOKUP_COLUMN = []
COMPOSITION_LOOKUP_COLUMN = []
MASS_LOOKUP_COLUMN = []
AMINO_ACID_DF = None

def loadRefSheet():

    print("Importing...")
 
    print("Loading reference URL")

    # Combined data from https://proteomicsresource.washington.edu/protocols06/masses.php and https://raw.githubusercontent.com/specht/proteomics-knowledge-base/master/amino-acids.csv
    url = "https://raw.githubusercontent.com/specht/proteomics-knowledge-base/master/amino-acids.csv"

    print("Getting the csv...")
    from os.path import exists

    global AMINO_ACID_DF
    global MASS_DF

    try:
        print("Reading amino acid names from local...")
        AMINO_ACID_DF = pd.read_csv("amino_acid_lookup.csv")
    except:
        print("No local lookup found; retreiving from online. This may take a while...")
        #AMINO_ACID_DF = pd.read_csv(url)
        #AMINO_ACID_DF.to_csv("amino_acid_lookup.csv", index = False)
        #print("File downloaded and ready for reference")

    global SINGLE_LOOKUP_COLUMN
    SINGLE_LOOKUP_COLUMN = AMINO_ACID_DF.loc[:, "single letter code"].tolist()
    
    global COMPOSITION_LOOKUP_COLUMN
    COMPOSITION_LOOKUP_COLUMN = AMINO_ACID_DF.loc[:, "composition"].tolist()

    global MASS_LOOKUP_COLUMN
    MASS_LOOKUP_COLUMN = AMINO_ACID_DF.loc[:, "average mass"].tolist()

    print("Imports finished!")

#--------------------------------------------------------------------------
