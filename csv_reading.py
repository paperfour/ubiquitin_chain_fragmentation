# CSV file reading

import pandas as pd


SINGLE_LOOKUP_COLUMN = []
COMPOSITION_LOOKUP_COLUMN = []
AMINO_ACID_DF = None

def loadRefSheet():

    print("Importing...")
 
    print("Loading reference URL")
    url = "https://raw.githubusercontent.com/specht/proteomics-knowledge-base/master/amino-acids.csv"

    print("Getting the csv...")
    from os.path import exists

    global AMINO_ACID_DF

    if (exists("amino_acid_lookup.csv")):
        print("Reading amino acids from local...")
        AMINO_ACID_DF = pd.read_csv("amino_acid_lookup.csv")
    else:
        print("No local lookup found; retreiving from online. This may take a while...")
        AMINO_ACID_DF = pd.read_csv(url)
        AMINO_ACID_DF.to_csv("amino_acid_lookup.csv", index = False)
        print("File downloaded and ready for reference")

    global SINGLE_LOOKUP_COLUMN
    SINGLE_LOOKUP_COLUMN = AMINO_ACID_DF.loc[:, "single letter code"].tolist()
    
    global COMPOSITION_LOOKUP_COLUMN
    COMPOSITION_LOOKUP_COLUMN = AMINO_ACID_DF.loc[:, "composition"].tolist()


    print("Imports finished!")

#--------------------------------------------------------------------------
