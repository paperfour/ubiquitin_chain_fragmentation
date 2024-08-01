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

g = 0

for i in range(0, 0):
    
  import time    

  t = time.time()

  testCZFragmentation("OLD_SLOW")

  o = time.time() - t

  print("Old method done in " + str(o) + " seconds")


  t = time.time()

  # Get the formatted data
  text = "substrate, 41(Ub, 6(Ub), 48(Ub)), 113(Ub, 63(Ub, 6(Ub)))"

  # Establish the substrate for the PTMs
  substrate = Protein(GFP_SEQ)

  # Parse the text into nested list data
  # target_form = [substrate, [41, [ubiquitin1, [6, [ubiquitin2]], [48, [ubiquitin3]]]], [113, [ubiquitin4, [63, [ubiquitin5, [6, [ubiquitin6]]]]]]]
  target_form = getNestedListFromString(text, substrate)


  data = genCZFragmentsFast(getProteinCollection(target_form), form = text)

  # Pass to csv!
  df = pd.DataFrame(data)
  df.to_csv("output/FAST_NEW_FINAL.csv", index=False)

  n = time.time() - t

  print("Old method done in " + str(n) + " seconds")


  print(str(o/n) + "times faster!")
  g += o/n

  print("================================")

print("Average: " + str(g/100) + " times faster!")

print("================= Starting generation! =================")

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

for collection in polyChains:

  progress += 1


  percent = progress / structAmount

  timeLeft = (time.time() - startTime) * (1 - percent) / percent

  print("Generating... %" + "%.2f" % round(percent * 100, 2) + " Estimated time: " + "%.2f" % round(timeLeft/60, 2) + " minutes remaining " + visual[progress % len(visual)], end="\r")


  dict = genCZFragmentsFast(collection, form = unparse(collection))

  for key in dict.keys():
    data[key] += dict[key]

print("")
print("Converting to csv...")


# ================================= END OF OPTIMIZATION

df = pd.DataFrame(data)


from datetime import datetime

# datetime object containing current date and time
now = datetime.now()
dt_string = now.strftime("output/OPTIMIZED_%d_%m_%Y_%H_%M_%S_out.csv")

df.to_csv(dt_string, index=False)

print("Done! See " + dt_string + " for data")
