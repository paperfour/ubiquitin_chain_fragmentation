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



# The following functions are used to find the CZ fragments of every possible chain arrangement for a ubiquitinated substrate

# Returns a list of indexes in a protein that a Ubiquitin can be attached to
# Check if only lysine should be considered or all
def ubSites(protein):
  out = []

  for a in range(1, protein.sequenceLength()):

    if protein.sequence[a-1] == "K":

      out.append(a)

  if out[-1] != protein.sequenceLength:
    out.append(protein.sequenceLength())

  # Remove all the sites where there is already something bound
  for pairedA in protein.children.keys():
    if pairedA in out:
      out.remove(pairedA)

  return out


# Returns a list of lists, the sub-lists each containing proteins linked to eachother in a unique way so that all chain forms are taken into account
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
