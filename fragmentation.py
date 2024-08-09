# Iteration of cleaving

from protein import *
from measurement_utils import *
from cleaving import *
from parsing import *

import pandas as pd


# Returns a dictionary csv form reference of all the possible CZ-ion cleavages in a system of proteins
# Input should be in the form of [substrate, [41, [ub]], [48, [ub]]] or similar
# ======================    IMPORTANT FUNCTION #1   ================================
def genCZFragments(originalCollection, form = ""):

  # Iterate through every combination of cleaved protein and cleaved site

  dataDict = {"form" : [],
              "cleavedProtein": [],
              "ionType": [],
              "ionIndex": [],
              "mass" : [],
              "chemicalFormula" : []
              }



  # Iterate through every protein in the collection
  for index in range(0, len(originalCollection)):

    # DON'T EDIT THIS OBJECT because it references the original collection that the copies are made of... don't mess with it
    ORIGINAL_PROTEIN = originalCollection[index]

    # Iterate through each cleavage site in the protein
    for n in range(1, ORIGINAL_PROTEIN.sequenceLength()):

      # Cleave the collection on target protein
      #TODO: Abstract the cleavage type for more flexible usage
      tempCollection = cleave(originalCollection, originalCollection[index], n, IonType.C)

      # Save the masses of the pieces
      # Filters the collection by only the cleaved ions
      for pIon in filter(lambda p : p.ionType != IonType.NONE,  tempCollection):

        #print(p.name + " : " + p.sequence + " : " + str(getMass(proteinChemicalFormula(p))) + " : " + str(p.ionType))

        pParent = getBaseProtein(pIon)
        pFormula = proteinChemicalFormula(pParent)
        pMass = round(getMass(pFormula), 3)

        dataDict["form"].append(form)
        dataDict["cleavedProtein"].append(getSequenceName(ORIGINAL_PROTEIN.sequence) + ": " + getProteinName(index))
        dataDict["ionIndex"].append(n if pIon.ionType in [IonType.A, IonType.B, IonType.C] else ORIGINAL_PROTEIN.sequenceLength() - n)
        dataDict["ionType"].append(str(pIon.ionType)[8:])
        dataDict["mass"].append(pMass)
        dataDict["chemicalFormula"].append(chemFormulaToString(pFormula))

  return dataDict



# A copy of genCZFragments() but lmited for optimization
def genCZFragmentsFast(originalCollection, form = ""):

  ORIGINAL_MASS = getMass(proteinChemicalFormula(getBaseProtein(originalCollection[0])))

  # Iterate through every combination of cleaved protein and cleaved site

  dataDict = {"form" : [],
              "cleavedProtein": [],
              "ionType": [],
              "ionIndex": [],
              "mass" : []
              }


  # Iterate through every protein in the collection
  for index in range(0, len(originalCollection)):

    # DON'T EDIT THIS OBJECT because it references the original collection that the copies are made of... don't mess with it
    ORIGINAL_PROTEIN = originalCollection[index]
    
    # Iterate through each cleavage site in the protein
    for n in range(1, originalCollection[index].sequenceLength()):

      # Cleave the collection on target protein
      cIonCollection = cIonCleave(originalCollection, originalCollection[index], n)

      cIonParent = getBaseProtein(cIonCollection[-1])

  
      cIonMass = getCFragmentMassFast(cIonParent)

      # Calculate the C fragment
      dataDict["form"].append(form)
      dataDict["cleavedProtein"].append(getSequenceName(ORIGINAL_PROTEIN.sequence) + ": " + getProteinName(index))
      dataDict["ionIndex"].append(n)
      dataDict["ionType"].append("C")
      dataDict["mass"].append(round(cIonMass, 3))

      # Calculate the Z fragment
      dataDict["form"].append(form)
      dataDict["cleavedProtein"].append(getSequenceName(ORIGINAL_PROTEIN.sequence) + ": " + getProteinName(index))
      dataDict["ionIndex"].append(ORIGINAL_PROTEIN.sequenceLength() - n)
      dataDict["ionType"].append("Z")
      dataDict["mass"].append(round(ORIGINAL_MASS - cIonMass, 3))


  return dataDict





# A sample usage of the genCZFragments
def testCZFragmentation(name):
  # Get the formatted data
  text = "substrate, 41(Ub, 6(Ub), 48(Ub)), 113(Ub, 63(Ub, 6(Ub)))"

  # Establish the substrate for the PTMs
  substrate = Protein(GFP_SEQ)

  # Parse the text into nested list data
  # target_form = [substrate, [41, [ubiquitin1, [6, [ubiquitin2]], [48, [ubiquitin3]]]], [113, [ubiquitin4, [63, [ubiquitin5, [6, [ubiquitin6]]]]]]]
  target_form = getNestedListFromString(text, substrate)


  # !TODO! Edit the parsing and child creation to support other proteins


  # Strip away the site attachment/parent/child info (already stored inside the protein object)
  formCollection = getProteinCollection(target_form)

  # Generate fragments off of the proteins formed by the parsing
  data = genCZFragments(formCollection, form = text)

  # Pass to csv!
  df = pd.DataFrame(data)
  df.to_csv("output/" + name + ".csv", index=False)


def testCZFragmentationFast(name):

  # Get the formatted data
  text = "substrate, 41(Ub, 6(Ub), 48(Ub)), 113(Ub, 63(Ub, 6(Ub)))"

  # Establish the substrate for the PTMs
  substrate = Protein(GFP_SEQ)

  # Parse the text into nested list data
  # target_form = [substrate, [41, [ubiquitin1, [6, [ubiquitin2]], [48, [ubiquitin3]]]], [113, [ubiquitin4, [63, [ubiquitin5, [6, [ubiquitin6]]]]]]]
  target_form = getNestedListFromString(text, substrate)

  # !TODO! Edit the parsing and child creation to support other proteins

  # Strip away the site attachment/parent/child info (already stored inside the protein object)
  formCollection = getProteinCollection(target_form)


  # Generate fragments off of the proteins formed by the parsing
  data = genCZFragmentsFast(formCollection, form = text)

  # Pass to csv!
  df = pd.DataFrame(data)
  df.to_csv("output/" + name + ".csv", index=False)



# The following functions are used to find the CZ fragments of every possible chain arrangement for a ubiquitinated substrate

# Returns a list of indexes in a protein that a Ubiquitin can be attached to
# Check if only lysine should be considered or all
def ubSites(protein):
  out = []

  for a in range(1, protein.sequenceLength()):

    if protein.sequence[a-1] == "K":

      out.append(a)

  # Add the N terminus as a peptide bond if not already there
  if out[-1] != protein.sequenceLength():
    out.append(protein.sequenceLength())

  # Remove all the sites where there is already something bound
  for pairedA in protein.children.keys():
    if pairedA in out:
      out.remove(pairedA)

  return out

# Returns a list of ubiquitin-bindable indexes in a sequence without care for children
def simpleUbSites(proteinSequence):

  out = []

  # Append all lysines except the last one (it will always be added)
  for a in range(1, len(proteinSequence)):

    # Subtract one because the first site should correspond to the number 1 instead of the string index 0
    if proteinSequence[a-1] == "K":
      out.append(a)

  # Add the N terminus as a peptide bond
  out.append(len(proteinSequence))

  return out

# Returns a list of lists, the sub-lists each containing proteins linked to eachother in a unique way so that all chain forms are taken into account
# This is great, but is problematic because it creates duplucates
# ======================    IMPORTANT FUNCTION #2    ================================
def genUbiquitinatedForms(substrate, ubiquitinCount):

  # This is the returned list of lists
  out = []

  # Level 0 contains only the substrate
  lv0List = []
  lv0List.append(Protein(substrate.sequence))

  # For level N... starting with 0
  def rec(substList, n = 0):

    # For each protein in level N...
    for lvNProtein in substList:
      # And for each available site in said protein...
      for lvNSite in ubSites(lvNProtein):

        # There is a variation!
        # AS LONG AS the SHAPE is not already present in N1

        # Create a copy of the list - each one will be given a unique branching structure
        lvN1List = copyProteinCollection(substList)

        # Add a different protein to each one
        lvN1List.append(Protein(UBIQUITIN_SEQ, parent = lvN1List[substList.index(lvNProtein)], site = lvNSite))

        # If there is another level to be done....
        if (n + 1 < ubiquitinCount):
          # Recurse
          rec(lvN1List, n + 1)
        else:
          # Otherwise, add the built list to the output!
          out.append(lvN1List)

  rec(lv0List)

  return out

from datetime import datetime
import pandas as pd
import os

# Returns a file path to a list of parsable forms
# TODO: Instead of recursing on UBIQUITIN COUNT recurse on MAX CHAIN DEPTH
def genUbiquitinatedFile(substrateSequence, ubiquitinCount):

  # This is the unique date and time code that will be used to identify the temp files
  dtCode = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")


  # Initiate the first level of recursion (only one form possibility for just the substrate)
  temp0df = pd.DataFrame(["form", "substrate"])
  temp0FileName = "temp/LV0_" + dtCode + "_temp.csv"
  temp0df.to_csv(temp0FileName, index=False, header=False)

  # Create the empty base temp files for all further levels of recursion
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
    #os.remove(prevFileName)
    


  return "temp/LV" + str(ubiquitinCount) + "_" + dtCode + "_temp.csv"


n = 2


#genUbiquitinatedFile(GFP_SEQ, n)


# Returns the amount of forms possible
def countUbiquitinatedForms(substrate, ubiquitinCount):

  # This is the returned value
  # With the amount of 
  out = 1

  totalSites = len(ubSites(substrate))

  for i in range(0, ubiquitinCount):
    
    # Every additional PTM can choose from any available site... so there are [sites available] times more forms
    out *= totalSites

    # One site is occupied by the protein and eight sites are found on the new protein
    # 8 - 1 = 7
    totalSites += 7

  return out
import pprint
pprint.pprint([unparse(getProteinCollection(i)) for i in genUbiquitinatedForms(Protein(GFP_SEQ), n)])
print(countUbiquitinatedForms(Protein(GFP_SEQ), n), " lines in the final form!")