# Uses the reference csv at https://github.com/specht/proteomics-knowledge-base/blob/master/amino-acids.csv for amino acid information
# Uses pandas for csv manipulation and reading

print("Program started")

from csv_reading import loadRefSheet
loadRefSheet()


from protein import *
from fragmentation import *
from parsing import unparse
import pandas as pd


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
        "chemicalFormula" : [],
        "mass" : []
        }


visual = ["|", "\\", "-", "/"]

for collection in polyChains:

  progress += 1


  print("Generating... %" + "%.2f" % round(progress / structAmount * 100, 2) + " " + visual[progress % len(visual)], end="\r")


  dict = genCZFragments(collection, form = unparse(collection))

  for key in dict.keys():
    data[key] += dict[key]

print("")
print("Converting to csv...")

df = pd.DataFrame(data)


from datetime import datetime

# datetime object containing current date and time
now = datetime.now()
dt_string = now.strftime("output/%d_%m_%Y_%H_%M_%S_out.csv")

df.to_csv(dt_string, index=False)

print("Done! See " + dt_string + " for data")
