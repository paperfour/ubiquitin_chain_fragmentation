from iontype import IonType
from csv_reading import SINGLE_LOOKUP_COLUMN, COMPOSITION_LOOKUP_COLUMN, MASS_LOOKUP_COLUMN
from protein import Protein

# Contains functions for retrieving the mass and chemical formula of a given protein


# Generate chemical formula
def chemicalFormula(protein: Protein):

  ionType = protein.ionType

  match ionType:
    case IonType.NONE:
    # H and OH for the N and C terminus
         element_dict = { "C": 0, "H": 2, "N": 0, "O": 1, "S": 0}

    case IonType.A:
    # H on the N terminus but a C and O lost from the cleavage
         element_dict = { "C": -1, "H": 1, "N": 0, "O": -1, "S": 0 }

    case IonType.B:
    # H on the N terminus and no peptide parts lost from the cleavage
         element_dict = { "C": 0, "H": 1, "N": 0, "O": 0, "S": 0 }

    case IonType.C:
    # H on the N terminus and an extra N and H from cleaving off a piece of the following amino acid
         element_dict = { "C": 0, "H": 2, "N": 1, "O": 0, "S": 0 }

    case IonType.X:
    # OH on the C terminus and receives the extra C and O from the complementary A ion
         element_dict = { "C": 1, "H": 1, "N": 0, "O": 2, "S": 0 }

    case IonType.Y:
    # OH on the C terminus and no peptide parts lost from the cleavage
         element_dict = { "C": 0, "H": 2, "N": 1, "O": 1, "S": 0 }

    case IonType.Z:
    # OH on the C terminus but loses an N and H to the complementary C ion
         element_dict = { "C": 0, "H": 0, "N": -1, "O": 1, "S": 0 }

    case _:
       raise Exception("Unexpected ion type!")


  # Remove an OH from the C terminus if the C terminus is bound to another residue (has a parent)

  if (protein.parentAttachmentSite != -1):
    element_dict["O"] -= 1
    element_dict["H"] -= 1


  # Remove an H for every child (the peptide bond discards an H)

  element_dict["H"] -= len(protein.children)

  # Convert the list of amino acids to a list of row indexes in the reference csv so that lookup is easier

  indexes = []

  for n in range(0, protein.sequenceLength()):

    try:
      # Add the row number of each amino acid to the
      indexes.append(SINGLE_LOOKUP_COLUMN.index(protein.sequence[n]))
    except:
      # Unknown amino acid index
      print("Unknown amino acid: " + protein.sequence[n])
      protein.append(21)

  # "Indexes" now is a list of row number

  # For each amino acid...
  for amino_acid_index in indexes:

    # Retrieve the chemical formula of the amino acid
    string_comp = COMPOSITION_LOOKUP_COLUMN[amino_acid_index]

    # Take the string format and covert it to a list of [element symbols followed by their amount]
    # The following SHOULD be done with regex but I don't know regex

    for key in list(element_dict.keys()):
      string_comp = string_comp.replace(key, " " + key)
      string_comp = string_comp.strip()

    list_comp = string_comp.split(" ")

    # Add each element of the amino acid to the running total of the protein
    for element in list_comp:

      element_symbol = element[0]
      # Add 1 if the symbol is alone
      element_frequency = 1
      if (len(element) > 1):
        element_frequency = int(element[1:])

      element_dict[element_symbol] += element_frequency

  return element_dict


# Finds the chemical formula of a protein but takes into account all child proteins
def proteinChemicalFormula(protein: Protein):

  # Establish the parent composition
  element_dict = chemicalFormula(protein)

  for child in protein.children.values():
    # Add the children's compositions to a running total

    formula = proteinChemicalFormula(child)
    for element in formula.keys():
      element_dict[element] += formula[element]

  return element_dict



# Generate mass from chemical formula
def getMass(formula: dict):

  carbon = 12.011 * formula["C"]
  hydrogen = 1.008 * formula["H"]
  nitrogen = 14.007 * formula["N"]
  oxygen =  15.999 * formula["O"]
  sulfur = 32.065 * formula["S"]

  return carbon + hydrogen + nitrogen + oxygen + sulfur



# Mass generation directly from a protein parent, optimized for CZ fragments
def getCFragmentMassFast(protein: Protein):
   
   # Start with a "stolen" NH on the cleaved edge (stolen from the Z fragment)
  mass = 15.015

  # A list of all proteins in the composite protein structure

  l = []

  def addProt(p: Protein, l : list):

     l.append(p)

     for child in p.children.values():
        addProt(child, l)

  addProt(protein, l)

  # l now contains every protein whose mass should be considered


  # Add the masses of each composite protein

  for piece in (l):
      
    for n in range(0, piece.sequenceLength()):

      # Add an H on the N terminus
      mass += 1.008

      try:
        # Add the row number of each amino acid to the
        mass += MASS_LOOKUP_COLUMN[SINGLE_LOOKUP_COLUMN.index(piece.sequence[n])]
      except:
          print("Unknown amino acid: " + protein.sequence[n])

    #lose an H and OH for each isopeptide bond
    mass -= 18.01528 * len(piece.children)
  
  return mass
        
   






# used to print chemical formula in a readable manner
def chemFormulaToString(dict: dict):

  textFormula = ""

  for key in list(dict.keys()):

    element_frequency = dict[key]

    if (element_frequency > 0):

      textFormula = textFormula + str(key) + str(element_frequency) + " "

  return textFormula

