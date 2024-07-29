from protein import *
from names import getProteinName

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




# Returns a text version of a given collection
# The list argument should contain only proteins
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
