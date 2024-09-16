# Fixes the...
#       substrate, (1, ub), (2, ub)
#       substrate, (2, ub), (1, ub)
# problem!

# Currently the fastest script

# Uses the reference csv at https://github.com/specht/proteomics-knowledge-base/blob/master/amino-acids.csv for amino acid information
# Uses pandas for csv manipulation and reading

print("Program started")

# Imports
from csv_reading import loadRefSheet
loadRefSheet()
from protein import *
from fragmentation import *
from parsing import unparse
import pandas as pd


print("================= Starting generation! =================")


# Get the input for number of ubiquitins
while True:
    try:
        u = int(input("How many ubiquitins? "))
        break
    except:
        print("That is not a number!")


print("Substrate: GFP")
print("Generating varied forms...")

# Generate all the possible forms for the given substrate and ubiquitin amount
polyChains = genUbiquitinatedFormsNoDupe(GFP_SEQ, u)

print("Done!")

# Update user about how many protein possibilities
structAmount = len(polyChains)
print(str(structAmount) + " chain structures generated")

# Variables for updating the user on progress
progress = 0
visual = ["|", "\\", "-", "/"]
import time
startTime = time.time()
# Change in code to False for a 7% speed increase
updateUser = True


# The data that will be added to in order to form the final csv
data = {"form" : [],
        "cleavedProtein": [],
        "ionType": [],
        "ionIndex": [],
        "mass" : []
        }


# Make a unique file handle for saving the output
from datetime import datetime
now = datetime.now()
FILENAME = now.strftime("output/NoDupe_%d_%m_%Y_%H_%M_%S_out.csv")

# Begin generating the fragments!
print("Beginning fragmentation... Output: " + FILENAME)

# ====================================================

# Writes to the csv on the fly by directly appending the data to the file
def saveToFile():
  
  global data

  # Append to the csv
  df = pd.DataFrame(data)
  df.to_csv(FILENAME, mode="a", index=False, header=True)

  # Reset the variable
  data = {"form" : [],
        "cleavedProtein": [],
        "ionType": [],
        "ionIndex": [],
        "mass" : []
        }


# Save every 10k fragments generated
fragmentationChunkSize = 10000

# Add a value to the dict util it reaches a the fragmentation chunk max size
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

  #Clear space in ram!
  if (progress % fragmentationChunkSize == 0):
    print("===================== Intermediate saving... =====================", end="\r")
    saveToFile()


# Iterate through each form
for strProtein in polyChains:
    
    # Parse the text into fragmentable list
    collection = getProteinCollection(getNestedListFromString(strProtein, Protein(GFP_SEQ)))

    # Fragment the proteins
    calcFragment(collection)


# Save one final time for the last rows not encapsulated by the modulo
print("")
print("Final saving...")
saveToFile()

# Done!
print("")
print("Fragmentation finished! In " + str((time.time() - startTime)) + " seconds")
print("See " + FILENAME + " for data!")



# Save the fragment data to an output csv
#df = pd.DataFrame(data)

#from datetime import datetime
#now = datetime.now()
#fileName = now.strftime("output/NoDupe_%d_%m_%Y_%H_%M_%S_out.csv")

#df.to_csv(fileName, index=False)

