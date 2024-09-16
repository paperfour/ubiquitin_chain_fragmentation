# !!!!!!!!!! Not updated anymore !!!!!!!!!!! use fast_gen_no_dupe.py
# 
# 
# 
# 
# Experiments to speed up generation time

# Uses the reference csv at https://github.com/specht/proteomics-knowledge-base/blob/master/amino-acids.csv for amino acid information
# Uses pandas for csv manipulation and reading

raise Exception("This is not the script you are looking for")

print("Program started")

from csv_reading import loadRefSheet
loadRefSheet()

from protein import *
from fragmentation import *
from parsing import unparse
import pandas as pd


#testCZFragmentationFast("NEW_VER")

print("================= Starting generation! =================")

while True:
  try:
    ubiquitinCount = int(input("How many ubiquitins? "))
    break
  except:
    print("That is not a number!")

# Create the substrate object
substrate = Protein(GFP_SEQ)
print("Substrate: GFP")

# Calculate the amount of different protein arrangements
structAmount = countUbiquitinatedForms(substrate, ubiquitinCount)
print(str(structAmount) + " different structures!")


print("Setting up csv...")

data = {"form" : [],
        "cleavedProtein": [],
        "ionType": [],
        "ionIndex": [],
        "mass" : []
        }

startdf = pd.DataFrame(data)

from datetime import datetime

# Set up the specific file name for generation
now = datetime.now()
dtCode = now.strftime("%d_%m_%Y_%H_%M_%S")
fileName = "output/OPTIMIZED_" + dtCode + "_out.csv"

startdf.to_csv(fileName, index=False)

print("Data will be outputted into: " + fileName)


# ================================= FRAGMENTATION BEGIN

visual = ["|", "\\", "-", "/"]

import time

startTime = time.time()

progress = 0

# Change in code to False for a 7% speed increase
updateUser = True


# Writes to the csv on the fly by directly appending the data to the file
def saveToFile():
  
  global data

  # Append to the csv
  df = pd.DataFrame(data)
  df.to_csv(fileName, mode="a", index=False, header=True)

  # Reset the variable
  data = {"form" : [],
        "cleavedProtein": [],
        "ionType": [],
        "ionIndex": [],
        "mass" : []
        }


# Save every 15k fragments generated
fragmentationChunkSize = 15000

# Adds values to the dict util it reaches a certain point
# In which case it will 'flush' the dict into the csv file and free up RAM
def calcFragment(collection):

  global progress
  progress += 1

  if (updateUser):

    percent = progress / structAmount

    timeLeft = (time.time() - startTime) * (1 - percent) / percent

    print("Generating... %" + "%.2f" % round(percent * 100, 2) + " Forms generated: " + str(progress) + " Estimated time: " + "%.2f" % round(timeLeft/60, 2) + " minutes remaining " + visual[progress % len(visual)], end="\r")

  newDataDict = genCZFragmentsFast(collection, form = unparse(collection))

  for key in data.keys():
    data[key] += newDataDict[key]

# Clear space in ram!
  if (progress % fragmentationChunkSize == 0):
    saveToFile()

print("Generating forms...")
forms = genUbiquitinatedForms(substrate, ubiquitinCount)

print("Starting fragmentation")
# Generate a fragment for every form
for form in forms:
    calcFragment(form)

# Save one final time for the last rows not encapsulated by the modulo
saveToFile()


  
print("")
print("Fragmentation finished! In " + str((time.time() - startTime)) + " seconds")

# Done!
print("See " + fileName + " for data")

# ===================================== FINISH =================================================



# ============== NON-RAM FORM MAKING... WIP

if 3 == 3:
  raise Exception("Not Done Yet!")

# Initiate the first level of recursion (only one possibility for just the substrate)
lv0List = [Protein(substrate.sequence)]
temp0df = pd.DataFrame([unparse(lv0List)])
temp0FileName = "temp/LV0_" + dtCode + "_temp.csv"
temp0df.to_csv(temp0FileName, index=False, header=False)

# Create the empty base temp files for all further levels of recursion
for n in range(0, ubiquitinCount):
  tempNdf = pd.DataFrame([])
  tempNFileName = "temp/LV" + str(n + 1) + "_"  + dtCode + "_temp.csv"
  tempNdf.to_csv(tempNFileName, mode="a", index=False, header=False)


# Instead of making a list with genUbiquitinatedForms(), the ram storage is no longer necessary by storing each level of recursion in a temp file
'''
more hardcoding = more speed

[K1, K2, K3, Nterminus] < Ubiquitin bonding sites

(sub, K1(ub))
(sub, K1(ub))
(sub, K1(ub))
(sub, K1(ub))

USE PANDAS CHUNKING INSTEAD OF READ AND POP


'''
# Level 0 contains only the substrate

# For level N... starting with 0
def rec(substList, n = 0):

  lvN1df = pd.DataFrame(list())
  lvN1Name = "temp/LV" + str(n + 1) + "_" +  dtCode + "_temp.csv"
  lvN1df.to_csv(lvN1Name)

  # For each protein in level N...
  for lvNProtein in substList:
    # And for each available site in said protein...
    for lvNSite in ubSites(lvNProtein):

      proteinCollec

      tempNdf = pd.DataFrame([unparse(proteinCollection)])
      tempNdf.to_csv(fileName, mode="a", index=False, header=False)

      # Create a copy of the list - each one will be given a unique branching structure
      lvN1List = copyProteinCollection(substList)

      # Add a different protein to each one
      lvN1List.append(Protein(UBIQUITIN_SEQ, parent = lvN1List[substList.index(lvNProtein)], site = lvNSite))

      # If there is another level to be done....
      if (n + 1 < ubiquitinCount):
        # Recurse
        rec(lvN1List, n + 1)
      else:
        # Otherwise, there is no more levels to be done. Thus, calculate the fragment!
        appendFragment(lvN1List)

# Begin the cycle
rec(temp0FileName)






