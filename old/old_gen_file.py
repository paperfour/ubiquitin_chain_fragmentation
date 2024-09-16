from datetime import datetime
from protein import *
from measurement_utils import *
from cleaving import *
from parsing import *
from fragmentation import *

import pandas as pd

import os



# Returns a file path to a list of parsable forms
# Produces duplciate forms, but should be used if its not just ubiquitin
def genUbiquitinatedFileDUPE(substrateSequence, ubiquitinCount):

  # This is the unique date and time code that will be used to identify the temp files
  dtCode = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")


  # Initiate the first level of iteration (only one form possibility for just the substrate)
  temp0df = pd.DataFrame(["form", "substrate"])
  temp0FileName = "temp/LV0_" + dtCode + "_temp.csv"
  temp0df.to_csv(temp0FileName, index=False, header=False)

  # Create the empty base temp files for all further levels of iteration
  for n in range(1, ubiquitinCount + 1):
    tempNdf = pd.DataFrame(["form"])
    tempNFileName = "temp/LV" + str(n) + "_"  + dtCode + "_temp.csv"
    tempNdf.to_csv(tempNFileName, index=False, header=False)
  
  # Add the ubiquitins

  # For every level...
  for n in range(1, ubiquitinCount + 1):
    
    print("Building level ", n)

    prevFileName = "temp/LV" + str(n - 1) + "_"  + dtCode + "_temp.csv"
    nFileName = "temp/LV" + str(n) + "_"  + dtCode + "_temp.csv"

    # Open and chunk the n-1 (previous level) file
    # TODO: Make chunk size dependent on the amount of ubiquitinatable sites (bigger chunks at the beginning)
    prevLevelChunks = pd.read_csv(prevFileName, chunksize=50000) 

    # For each chunk in the previous level...
    for chunk in prevLevelChunks:

      print("Chunking level " + str(n - 1) + " to make level " + str(n))
      
      # Read the chunk into a dataframe
      chunk_df = pd.DataFrame(chunk)
      
      # Establish a value to be saved at the end of the chunk reading
      nSavedChunk = []
      
      print("Reading lines in chunk")
      # For each line...
      for line in chunk_df["form"]:

        # Create a ubiquitin that will try every avaiable site
        newProtein = simpleUbiquitin()

        # Initialize the full protein collection complete with the new ubiquitin
        substrate = SimpleProtein(simpleUbSites(substrateSequence))
        fullProtein = simpleParse(line, substrate)

        # To add the n level ubiquitin:
          # for each protein..
          #   for each open site
          #     vary!

    
        # Add the new protein
        fullProtein.append(newProtein)
     
        # Vary the new protein for each open site

        # For each protein in the collection that isn't the new protein...
        for compositeProtein in filter(lambda x: x != newProtein, fullProtein):

          # Assign temporarily the new protein as each composite's child (this is ok even if there are no open sites because the saving is done inside the iteration over open sites)
          newProtein.parent = compositeProtein

          # For each open site...
          for site in compositeProtein.attachmentSites:
            
            # Temporarily set the site to said open site
            newProtein.parentSite = site
        
            # Unparse to text
            text = simpleUnparse(fullProtein)

            # Write to n level list saving the setup into the list as a string
            nSavedChunk.append(text)


 
      # Save the numerous variations of the chunk to the csv after every chunk
      n_df = pd.DataFrame(nSavedChunk)
      n_df.to_csv(nFileName, mode="a", index=False, header=False)

    # Delete previous temp file
    print("Deleting level " + str(n - 1))
    os.remove(prevFileName)
    


  return "temp/LV" + str(ubiquitinCount) + "_" + dtCode + "_temp.csv"
