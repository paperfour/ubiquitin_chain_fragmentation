# -*- coding: utf-8 -*-
'''

!!! IMPORTANT

This is the OLD code from a GOOGLE COLAB workspace... JUST FOR REFERENCE THIS SHOULD NOT BE USED

'''


# Uses the reference csv at https://github.com/specht/proteomics-knowledge-base/blob/master/amino-acids.csv for amino acid information
# Uses pandas for csv manipulation and reading

raise Exception("Hey, you're using the wrong file!")


print("Program started")


# CSV file reading

print("Importing...")
import pandas as pd
import math

print("Loading reference URL")
url = "https://raw.githubusercontent.com/specht/proteomics-knowledge-base/master/amino-acids.csv"

print("Getting the csv...")
from os.path import exists

if (exists("amino_acid_lookup.csv")):
  print("Reading amino acids from local...")
  aminoAcidDF = pd.read_csv("amino_acid_lookup.csv")
else:
  print("No local lookup found; retreiving from online. This may take a while...")
  aminoAcidDF = pd.read_csv(url)
  aminoAcidDF.to_csv("amino_acid_lookup.csv", index = False)
  print("File downloaded and ready for reference")

single_lookup_column = aminoAcidDF.loc[:, "single letter code"].tolist()

print("Imports finished!")

#--------------------------------------------------------------------------



# A silly dictionary of names used for debugging
NAMES = ["Breeanna","Fiona","Jailene","Jacquelin","Marla","Terri","Chance","Lilianna","Hamza","Harley","Erin","Zack","Meru","Colin","Sandra","Olivia","Laci","Sebastien","Joel","Mya","Long","Mae","Ayesha","Michael","Mia","Stephanie","Wilson","Charlotte","Shakira","Maleah","Keegan","Hayley","Johanna","Greta","Broderick","Joan","Magaly","Jeff","Keelan","Estefania","Yesica","Reilly","Morris","Haylie","Brent","Jacquelyn","Gregory","Deshaun","Gregorio","Tayla","Lucy","Adrian","Lisbeth","Hailie","Clarissa","Micheal","Clay","Katharine","Keyla","Karyme","Taylor","Samuel","Priscila","Justus","Brieanna","Naomi","Zoie","Haley","Gisselle","Chelsey","Santana","Brycen","Dominique","Ayla","Destinee","Amie","Reese","Keanu","Caylee","Ryley","Alani","Maureen","Graciela","Lyle","Amaris","Tatianna","Deasia","Desmond","Bobby","Cornelius","Asha","Madisen","Ramiro","Niko","Shayne","Jairo","Dajuan","Arjun","Bobbie","Amir","Maddie","Camden","Adrienne","Roberto","Jamir","Gissell","Keyon","Johnna","Peyton","Milan","Angie","Campbell","Amelia","Ronan","Jadon","Alaina","Philip","Christianna","Carleigh","Kianna","Zakary","Raina","Payton","Tracey","Janell","Kristyn","Tyra","Kole","London","Jovanny","Breonna","Annalisa","Brooks","Jerome","Kelvin","Mckinley","Brenna","Dora","Brandt","Kaylynn","Joseph","Dorien","Makenzi","Korey","Breann","Carmen","Crystal","Deonte","Travon","Cory","Parker","Jack","Jessika","Christy","Haily","Mariam","Zechariah","Jennifer","Brody","Leigha","Daylon","Emerson","Jonas","Quincy","Kobie","Alexandra","Maribel","Tess","Rodney","Shianne","Tania","Davon","Jaqueline","Rhiannon","Shae","Brian","Arnold","Josue","Kyara","Elias","Malka","Marquise","Chris","Tyreese","Dylon","Jermaine","Katherine","Dallin","Romeo","Areli","Bayleigh","Christiana","Zane","Jeffrey","Halle","Kurtis","Nora","Ashley","Rachelle","Gunnar","Kalia","Simone","Breanne","Byron","Laney","Marcelo","Tucker","Micayla","Jamaal","Desirae","Eliana","Mckenzie","Dalila","Elliott","Travis","Taliyah","Miah","Dejuan","Sonny","Jania","Dakoda","Vivian","Ximena","Darian","Tristian","Freddie","Allissa","Nick","Karlee","Javier","Forrest","Monet","Kya","Brennan","Abram","Litzy","Tabatha","Demetrius","Koby","Natalie","Kerry","Deandre","Jason","Darien","Dale","Obed","Juan","Joanne","Charlize","Valentin","Rosalinda","Yehuda","Jamari","Valencia","Hudson","Madilyn","Sara","Jorden","Annemarie","Jonatan","Lexus","Yusuf","Tariq","Oswaldo","Marlee","Tavon","Kody","Charlie","Zachery","Jaidyn","Mario","Aliya","Justice","Ammon","Ariel","Felipe","Myra","Armando","Melody","Ivan","Grayson","Amia","Kiley","Lina","Tyrone","Cristopher","Chaya","Rebecca","Cassandra","Ira","Alfredo","Lilian","Anjali","Herman","Keyonna","Rohan","Sammy","Loren","Jazmine","Demond","Karime","Luc","Michaela","Antonio","Harper","Lizet","Remington","Gabriel","Cristina","Chloe","Dario","Leann","Miya","Rhett","Annabel","Willie","Ladarius","Dahlia","Carl","Addie","Corina","Beth","Rafael","Tamya","Helena","Alanis","Delia","Maximillian","Arely","Aisha","Gissel","Colton","Rebeca","Aditya","Amari","Deon","Yolanda","Brisa","Esther","Tyquan","Darrian","Whitney","Reece","Jensen","Jimmie","Katarina","Alexus","Giovani","Holden","Erica","Danae","Melissa","Reginald","Francisco","Wilfredo","Sheridan","Giovanni","Antwan","Harold","Marcus","Jarret","Fabiola","Cruz","Treasure","Miguelangel","Cameron","Larry","Adela","Johana","Aylin","Billy","Collin","Coy","Annalee","Notnamed","Lara","Holland","Teresa","Kassandra","Kendrick","Amanda","Candice","Ray","Diane","Carrie","Ian","Kyron","Mara","Sylvia","Tea","Vladimir","Ayanna","Alisha","Kristofer","Devin","Trinity","Jocelynn","Jailyn","Gwyneth","Rayven","Keira","Kendell","Kentrell","Luis","Jolene","Mark","Melinda","Brylee","Tabitha","Tierra","Harrison","Kaylyn","Nia","Randall","Breana","Samaria","Braydon","Sariah","Zack","Brook","Danika","Devon","Jaylon","Aubrie","Fernando","Jennie","Antony","Kamren","Prince","Jesse","Rachel","Shanna","Fletcher","George","Lidia","Walter","Kane","Hayleigh","Sheldon","Avery","Leonel","Monique","Citlalli","Kelton","Alexis","Tyrell","Abdul","Clint","Beyonce","Roy","Catherine","Marion","Kobe","Susan","Tyrique","Pilar","Allyssa","Gino","Heidi","Yulisa","Angela","Alonzo","Presley","Riley","Janessa","Savana","Keshon","Neha","Mason","Jamar","Louie","Chantel","Carlee","Lawson","Elana","Brooklynn","Carter","Malia","Courtney","Steven","Irvin","Conrad","Rayshawn","Francesca","Chynna","Jacie","Mohammad","Toni","Gaige","Janaya","Abriana","Barbara","Jude","Dashawn","Ciarra","Braiden","Ken","Addison","Nikki","Austen","Lance","Juanita","Marguerite","Jaron","Nathalie","Katelyn","Jeremiah","Carolina","Destinie","Gabriella","Sarina","Shyla","Alissa","Daron","Carlo","Dylan","Vincenzo","Gabriela","Veronica","River","Adolfo","Geovanni","Seth","Alisa","Sarahi","Earl","Josephine","Anissa","Jamal","Kelsey","Zaira","Christion","Robert","Francis","Garret","Joaquin","Sophia","Brannon","Malcolm","Mona","Liberty","Bryce","Tyshawn","Jacob","Alaysia","Trevion","Dajah","Jerod","Codey","Bryn","Jerrod","Ramon","Tatiana","Mathew","Ciana","Analise","Elisa","Wesley","Iman","Jaylynn","Starr","Jaila","Savanah","Ezra","Ana","Dustyn","Cori","Noe","Tobias","Marc","Taj","Malorie","Angelina","Guillermo","Marcanthony","Kaylen","Kirstyn","Kaylea","Larissa","Frank","Jevon","Dianna","Ross","Kaylee","Everett","Janelle","Jean","Hailee","Lydia","Anna","Allysa","Patrick","Lana","Edmund","Noel","Markell","Kaleigh","Liliana","Issac","Theresa","Jocelyne","Moses","Makena","Matteo","Hailey","Jala","Juancarlos","Zoe","Maryam","Analisa","Kaya","Helen","Jaime","Ireland","Madisyn","Rahul","Gideon","Magali","Shaelyn","Amirah","Kyla","Miguel","Alannah","Frankie","Jajuan","Keona","Anjelica","Eileen","Britney","Ernest","Shania","Bella","Sahil","Shirley","Jaleel","Kearra","Paul","Denisse","Keven","Jesenia","Cole","Raphael","Katie","Rogelio","Jalisa","Giana","Willow","Josie","Sandy","Judy","Jedidiah","Annabelle","Stefani","Cielo","Alek","Hernan","Logan","Tyreek","Brionna","Ahmed","Nicklaus","Raelynn","Katia","Margaret","Sabina","Marie","Lindsay","Giancarlo","Kaitlynn","Fredrick","Dequan","Randolph","Marlene","Myah","Destany","Esmeralda","Kira","Tyasia","Mina","Ezekiel","Hunter","Alexandro","Jaret","Bruno","Alessandro","Natalee","Warren","Cade","Ashanti","Joann","Isai","Cameryn","Tahj","Christin","Kerri","Nathanial","Drew","Kent","Gemma","Amos","Martin","Aldo","Juwan","Calvin","Diego","Salma","Kaci","Luisa","Janet","Cullen","Brennon","Tyla","Jaelyn","Jena","Stone","Jaiden","Reid","Ethan","Evan","Jamarcus","Shivani","Caley","Belinda","Marques","Lilly","Cristobal","Sahara","Varun","Junior","Annie","Marcellus","Daveon","Branson","Mika","Lena","Cristofer","Davion","Cecil","Colby","Sienna","Wendy","Ashlynn","Jamison","Kristi","Santino","Siena","Margarita","Nehemiah","Annette","Leonardo","Leo","Carlos","Anaya","Darrion","Lars","Rianna","Shane","Ricky","Josiah","Chandler","Shelby","Everardo","Azaria","Magdalena","Riya","Kieran","Klarissa","Abagail","Shea","Keshaun","Ashtyn","Nickolas","Devan","Wade","Kristina","Markel","Kyree","Lela","Dejon","Louis","Alexandria","Erik","Abigale","Jakobe","Johnathon","Jon","Rene","Kyra","Felix","Adonis","Anita","Sydnee","Ahmad","Maryann","Danny","Ishmael","Jessie","Skylar","Darrien","Jacinda","Theron","Darin","Anthony","Kayli","Rolando","Mikaila","Karly","Bilal","Leila","Alfonso","Andy","Marissa","Brenton","Sana","Carla","Corey","Reuben","Imani","Kellen","Fatima","Willis","Perry","Sarah","Lacey","Jabari","Phoenix","Iyana","Shamar","Shawn","Yamileth","Bailey","Serina","Kolton","Anika","Harry","Breanna","Oscar","Treyvon","Garrison","Ali","Angeles","Jayda","Kellie","Griffin","Ainsley","Samara","Priya","Mireya","Bernard","Rosalie","Vanesa","Phoebe","Charles","Bridget","Travion","Ashlee","Lynette","Araceli","Aislinn","Jared","Laken","Dionte","Kori","Grant","Louise","Beatrice","Athena","Madelyne","Diamond","Yahaira","Edgar","Keith","Jackie","Martha","Adia","Richard","Nathaniel","Trevon","Kyndall","Noa","Marin","Jericho","Jackeline","Arianna","Odalis","Lane","Marco","Lianna","Shekinah","Ean","Sebastian","Kassidy","Jasmine","Lee","Shreya","Carli","Mackenzie","Terrence","Kayleen","Arlene","Maegan","Tanner","Randy","Christen","Gilberto","Jelani","Ameer","Kalista","Lauren","Rivka","Arielle","Simon","Dorothy","Lexis","Milton","Jakayla","Trace","Axel","Jaxson","Ignacio","Infant","Shyann","Mustafa","Hugh","Janine","Gina","Augustus","Chaz","Emilio","Dallas","Noah","Kobi","Elaina","Rory","Jovany","Stefan","Rubi","Theodore","Derick","Christopher","Francesco","Keon","Jamia","Kaiya","Jalynn","Devontae","Delanie","Abigayle","Andrea","Erich","Carson","Lynn","Emmanuel","Frances","Terrell","Lacy","Armand","Jessalyn","Joe","Harmony","Lazaro","Krysta","Mekhi","Mariano","Linda","Yesenia","Patience","Celeste","Ivette","Sincere","Ava","Divya","Catrina","Morgan","Cristal","Seamus","Farrah","Leeann","Shyanne","Damion","Jadyn","Jaeden","Jazmyne","Allen","Luiz","Anton","Maia","Maurice","Dan","Jacklyn","Mykayla","Alexzander","Julian","Daija","Kate","Alice","Sterling","Yosef","Eliezer","Zackary"]

from enum import Enum

IonType = Enum('IonType', ['NONE', 'A', 'B', 'C', 'X', 'Y', 'Z'])


def getProteinName(index):
  try:
    return NAMES[index]
  except:
    return "Bob #" + str(index - len(NAMES))
  


# Protein class, stores helpful information related to proteins

class Protein:


  def __init__(self, singleLetterSequence, parent = None, site = -1, ionType = IonType.NONE, updateParent = True):

    # Key-value pairs are [attachment_site : protein_child]
    self.children = {}

    self.sequence = singleLetterSequence
    self.parent = parent
    self.ionType = ionType

    self.setParent(parent, site)

    # Assign a name for debug purposes
    self.name = str(id(self))

  def setParent(self, proteinParent, site):

    self.parent = proteinParent
    self.parentAttachmentSite = site

    # Update the parent's children
    if (proteinParent != None):

      # Attach self to proteinParent at site

      if (site == -1):
        raise Exception("Cannot specify parent without an attachment amino acid!")
      elif(site > proteinParent.sequenceLength()):
        raise Exception("Amino acid doesn't exist!")
      else:
        # Assign itself as the parent's child
        proteinParent.children[site] = self


  # For debugging
  def __str__(self):
    return("I am " + self.name + ". Attached to " + ("nobody; I am on top of the world!" if (self.parent == None) else self.parent.name) + " at " + str(self.parentAttachmentSite))

  # Returns the number of amino acids in the peptide (not counting parent/children)
  def sequenceLength(self):
    return len(self.sequence)


# Constant protein sequences
UBIQUITIN_SEQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
GFP_SEQ = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

SEQUENCE_DICT = {UBIQUITIN_SEQ:"Ubiquitin", GFP_SEQ:"GFP"}

# Returns the common name for a protein sequence if it is included in the SEQUENCE_DICT
# Otherwise returns "Unknown protein"
def getSequenceName(sequence):
  if sequence in SEQUENCE_DICT.keys():
    return SEQUENCE_DICT[sequence]
  else:
    return "Unknown protein"


# Generate chemical formula
def chemicalFormula(protein):

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
      indexes.append(single_lookup_column.index(protein.sequence[n]))
    except:
      # Unknown amino acid index
      protein.append(21)

  # "Indexes" now is a list of row numbers

  comp_lookup_column = aminoAcidDF.loc[:, "composition"].tolist()

  # For each amino acid...
  for amino_acid_index in indexes:

    # Retrieve the chemical formula of the amino acid
    string_comp = comp_lookup_column[amino_acid_index]

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
def proteinChemicalFormula(protein):

  # Establish the parent composition
  element_dict = chemicalFormula(protein)

  for child in protein.children.values():
    # Add the children's compositions to a running total

    formula = proteinChemicalFormula(child)
    for element in formula.keys():
      element_dict[element] += formula[element]

  return element_dict



# Generate mass from chemical formula
def getMass(formula):

  carbon = 12.011 * formula["C"]
  hydrogen = 1.008 * formula["H"]
  nitrogen = 14.007 * formula["N"]
  oxygen =  15.999 * formula["O"]
  sulfur = 32.065 * formula["S"]

  return carbon + hydrogen + nitrogen + oxygen + sulfur


# used to print chemical formula in a readable manner
def chemFormulaToString(dict):

  textFormula = ""

  for key in list(dict.keys()):

    element_frequency = dict[key]

    if (element_frequency > 0):

      textFormula = textFormula + str(key) + str(element_frequency) + " "

  return textFormula


# Parses and prepares a nested list from given text
def getNestedListFromString(text, substrate):


  # Parses the creation of nested proteins
  # Returns the text parsed into nested lists
  # Text must be in the following form: " substrate, site#( ---- ),  site#( ---- ), ...  " where ---- can be replaced with a similar structure
  def parse(string):

    s = string.replace(" ", "")

    nested_list_protein = []

    # Count is used to keep track of when parenthesis are "neutral"... this indicates that the initial level has been reached
    count = 0

    # Pretending that there was a comma right before the whole string so that it can be trimmed neatly
    lastIndex = -1


    # Divide the string by commas on a shared parenthesis level
    for index in range(0, len(s)):

      # s[index] is the character that is being examined

      if s[index] == '(':
        count += 1
      elif s[index] == ')':
        count -= 1
      elif (s[index] == ',') and count == 0:

        nested_list_protein.append(s[(lastIndex + 1):index])

        lastIndex = index

    # Add the last piece (tail end)
    nested_list_protein.append(s[lastIndex + 1:])


    # Split up children (if any) into [where they bind] and [what they are]
    for n in range(1, len(nested_list_protein)):

      tempSiteInfo = nested_list_protein[n]

      bindingSite = int(tempSiteInfo[:tempSiteInfo.index("(")])

      whatIsBound = tempSiteInfo[tempSiteInfo.index("("):]

      # Recurse on each child protein, but chop off the surrounding parenthesis for the text
      nested_list_protein[n] = [bindingSite, parse(whatIsBound[1:(len(whatIsBound) - 1)])]


    return(nested_list_protein)

  # Create child proteins for X in "X, n(Y), n(Z)..." from the nested lists using recursion
  # parentLevel is a nested list that begins "[X, ..."
  # Updates the parentLevel list with protein objects replacing the strings
  def createChildProteins(parentLevel, substrate):

      if substrate != None:
        parentLevel[0] = substrate

      for childInfo in parentLevel[1:]:

        # Hardcoded UBIQUITIN_SEQ can be replaced by reading childInfo[1][0] and finding the right protein
        newProtein = Protein(UBIQUITIN_SEQ, parent=parentLevel[0], site=childInfo[0])

        # Updates the nested list representation to include the unique protein model
        childInfo[1][0] = newProtein

        # Can input None because the 'substrate' is was just set to the newProtein and does not need to be set
        createChildProteins(childInfo[1], None)


  # The right structure, but with strings instead of protein objects
  stringList = parse(text)

  # Injects the protein objects into the structure
  createChildProteins(stringList, substrate)

  # Return
  return stringList


# Returns a flat list of every protein in a nested list... makes it easier to iterate
def getProteinCollection(nestedList, final = True):

  proteinCollection = []

  for entry in nestedList:

    if isinstance(entry, list):
      proteinCollection += getProteinCollection(entry)
    elif isinstance(entry, Protein):
      proteinCollection.append(entry)


  return proteinCollection


# Returns identical proteins but different objects... for the copy() method isn't shallow enough because the parent and child pointers are still to the uncopied objects
def copyProteinCollection(proteinCollection):

  out = []

  for p in proteinCollection:

    # Make copies of each protein, but with no
    out.append(Protein(p.sequence, parent = None, site = p.parentAttachmentSite, ionType=p.ionType))

  for i in range(0, len(out)):

    # Finds where the corresponding parent is and assigns the recently created version of it to each protein

    counterpart = proteinCollection[i]

    if counterpart.parentAttachmentSite != -1:

      counterpartParentIndex = proteinCollection.index(counterpart.parent)

      out[i].setParent(out[counterpartParentIndex], counterpart.parentAttachmentSite)

  return out


# Returns a new collection of proteins that have been cleaved at a specified site of a specified protein
# Cleavage site is which amino acid to be cleaved AFTER (i.e. the last amino acid in the N-terminus fragment)
# For instance: AAGG cleaved at 3 would yield AAG + G
#               AAGG cleaved at 1 would yield A + AGG
def cleave(referenceCollection, referenceCleavedProtein, cleavageSite: int, nTerminusIonType: IonType):

  proteinCollection = copyProteinCollection(referenceCollection)
  cleavedProtein = proteinCollection[referenceCollection.index(referenceCleavedProtein)]

  # Ensure that the protein is long enough to be cleaved at the requested site
  if ((cleavageSite >= cleavedProtein.sequenceLength()) or (cleavageSite < 1)):
    raise Exception("Cannot cleave at " + str(cleavageSite))

  # Assign the correct complimentary ion to the cleavage type
  match nTerminusIonType:

    case IonType.A:
      cTerminusIonType = IonType.X

    case IonType.B:
      cTerminusIonType = IonType.Y

    case IonType.C:
      cTerminusIonType = IonType.Z

    case _:
       raise Exception("Unsupported cleavage type")



  # Create the C-terminus (x, y, or z ion) fragment; the parent is retained for this fragment as are any children right of the cleavage site

  cFragment = Protein(cleavedProtein.sequence[(cleavageSite):], ionType = cTerminusIonType, parent = cleavedProtein.parent, site = cleavedProtein.parentAttachmentSite)

  # Create the N-terminus (a, b, or c ion) fragment; this fragment will be parentless by nature but still retains children left of the cleavage site

  nFragment = Protein(cleavedProtein.sequence[:(cleavageSite)], ionType = nTerminusIonType)

  # Reassign the old protein's children to the correct ions

  for attachmentSite in cleavedProtein.children.keys():

    # If left of cleavage site
    if (attachmentSite <= cleavageSite):
      # Attach to N fragment
      cleavedProtein.children[attachmentSite].setParent(nFragment, attachmentSite)
    # The rest must be right of cleavage site
    else:
      # Attach to C fragment and set the new attachment site to be based of the new leftmost amino acid
      cleavedProtein.children[attachmentSite].setParent(cFragment, attachmentSite - cleavageSite)

  # Replace the old protein with the two cleaved ones in the protein collection
  # print(cleavedProtein.name + " has split into " + nFragment.name + " and " + cFragment.name)
  proteinCollection.remove(cleavedProtein)
  proteinCollection.append(nFragment)
  proteinCollection.append(cFragment)

  return proteinCollection

# Returns the parent of the parent of... a given protein
# Goes down the tree until the selected node has no parent itself
def getBaseProtein(protein):

  if(protein.parentAttachmentSite == -1):
    return protein
  else:
    return getBaseProtein(protein.parent)


# Prints the structure of the proteins in a nested list
def printStruct(protein, enclosingCollection, indent = 0):

  try:
    print(("  " * indent) + getProteinName(enclosingCollection.index(protein)) + " (" + getSequenceName(protein.sequence) + "): " + str(protein.parentAttachmentSite))
  except:
    print("Protein is not part of the enclosing collection. This is normally caused by a child not part of the parent's collection")


  for child in protein.children.values():
    printStruct(child, enclosingCollection, indent = indent + 1)


# Returns a text version of a given collection
# The list will have proteins not strings
# Uses a similar method to printStruct but replaces indentations with parenthesis
output = ""
def unparse(enclosingCollection):


  global output
  output = ""

  maxParent = enclosingCollection[0]

  def rec(protein):

    #try:
    string = str(protein.parentAttachmentSite) + ("(") + getSequenceName(protein.sequence)  + ": " + getProteinName(enclosingCollection.index(protein)) + ", "

    global output
    output += string

   # except:
    #  print("Protein is not part of the enclosing collection. This is normally caused by a child not part of the parent's collection")

    for child in protein.children.values():
      rec(child)

    output += "), "


  # Recurse for each child
  rec(maxParent)

  # Format correctly
  output = output[3:(len(output) - 5)]

  output = output.replace(", )", ")")


  return output


# Returns a dictionary csv form reference of all the possible CZ-ion cleavages in a system of proteins
# Input should be in the form of [substrate, [41, [ub]], [48, [ub]]] or similar
# ======================    IMPORTANT FUNCTION #1   ================================
def genCZFragments(originalCollection, form = ""):

  # Iterate through every combination of cleaved protein and cleaved site

  dataDict = {"form" : [],
              "cleavedProtein": [],
              "ionType": [],
              "ionIndex": [],
              "chemicalFormula" : [],
              "mass" : []
              }



  # Iterate through every protein in the collection
  for index in range(0, len(originalCollection)):


    # Iterate through each cleavage site in the protein
    for n in range(1, originalCollection[index].sequenceLength()):

      # Cleave the collection on target protein
      tempCollection = cleave(originalCollection, originalCollection[index], n, IonType.C)

      # Save the masses of the pieces
      # Filters the collection by only the cleaved ions
      for pIon in filter(lambda p : p.ionType != IonType.NONE,  tempCollection):

        #print(p.name + " : " + p.sequence + " : " + str(getMass(proteinChemicalFormula(p))) + " : " + str(p.ionType))

        pParent = getBaseProtein(pIon)
        pFormula = proteinChemicalFormula(pParent)
        pMass = round(getMass(pFormula), 3)

        # DON'T EDIT THIS OBJECT because it references the original collection that the copies are made of... don't mess with it
        ORIGINAL_PROTEIN = originalCollection[index]

        dataDict["form"].append(form)
        dataDict["cleavedProtein"].append(getSequenceName(ORIGINAL_PROTEIN.sequence) + ": " + getProteinName(index))
        dataDict["ionIndex"].append(n if pIon.ionType in [IonType.A, IonType.B, IonType.C] else ORIGINAL_PROTEIN.sequenceLength() - n)
        dataDict["ionType"].append(str(pIon.ionType)[8:])
        dataDict["chemicalFormula"].append(chemFormulaToString(pFormula))
        dataDict["mass"].append(pMass)

  return dataDict



# A sample usage of the genCZFragments
def testCZFragmentation():
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
  data = genCZFragments(formCollection)

  # Pass to csv!
  df = pd.DataFrame(data)
  df.to_csv('out.csv', index=False)


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
df.to_csv('out.csv', index=False)



print("Done! See out.csv for data")
