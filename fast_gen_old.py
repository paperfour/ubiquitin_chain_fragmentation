# Experiments to speed up generation time

# Uses the reference csv at https://github.com/specht/proteomics-knowledge-base/blob/master/amino-acids.csv for amino acid information
# Uses pandas for csv manipulation and reading

print("Program started")

from csv_reading import loadRefSheet
loadRefSheet()


from protein import *
from fragmentation import *
from parsing import unparse
import pandas as pd


#testCZFragmentationFast("NEW_FAST_FIXED")
#testCZFragmentation("OLD_VER")

print("================= Starting NON-APPENDING generation! =================")

while True:
  try:
    u = int(input("How many ubiquitins? "))
    break
  except:
    print("That is not a number!")

print("Substrate: GFP")
print("Generating varied forms...")
polyChains = genUbiquitinatedForms(Protein(GFP_SEQ), u)


structAmount = len(polyChains)

print(str(structAmount) + " chain structures generated!")

print(str(countUbiquitinatedForms(Protein(GFP_SEQ), u)))

print("Beginning fragmentation...")


progress = 0

data = {"form" : [],
        "cleavedProtein": [],
        "ionType": [],
        "ionIndex": [],
        "mass" : []
        }


visual = ["|", "\\", "-", "/"]

# ================================= OPTIMIZATION BEGIN

import time

startTime = time.time()

# Change in code to False for a 7% speed increase
updateUser = True

for collection in polyChains:

  progress += 1

  if (updateUser):

    percent = progress / structAmount

    timeLeft = (time.time() - startTime) * (1 - percent) / percent

    print("Generating... %" + "%.2f" % round(percent * 100, 2) + " Estimated time: " + "%.2f" % round(timeLeft/60, 2) + " minutes remaining " + visual[progress % len(visual)], end="\r")

  dict = genCZFragmentsFast(collection, form = unparse(collection))

  
  for key in dict.keys():
    data[key] += dict[key]





print("")
print("Fragmentation finished! In " + str((time.time() - startTime)) + " seconds")
print("Converting to csv...")


# ================================= END OF OPTIMIZATION

df = pd.DataFrame(data)


from datetime import datetime

# datetime object containing current date and time
now = datetime.now()
dt_string = now.strftime("output/OLD_OPTIMIZED_%d_%m_%Y_%H_%M_%S_out.csv")

df.to_csv(dt_string, index=False)

print("Done! See " + dt_string + " for data")

