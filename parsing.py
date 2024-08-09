from protein import *
from names import getProteinName


# Parses the creation of nested proteins
# Returns the text parsed into nested lists
# Text must be in the following form: " substrate, site#( ---- ),  site#( ---- ), ...  " where ---- can be replaced with a similar structure
def parse(string):

  s = string.replace(" ", "")

  nestedList = []

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

      nestedList.append(s[(lastIndex + 1):index])

      lastIndex = index

  # Add the last piece (tail end)
  nestedList.append(s[lastIndex + 1:])


  # Split up children (if any) into [where they bind] and [what they are]
  for n in range(1, len(nestedList)):

    tempSiteInfo = nestedList[n]

    bindingSite = int(tempSiteInfo[:tempSiteInfo.index("(")])

    whatIsBound = tempSiteInfo[tempSiteInfo.index("("):]

    # Recurse on each child protein, but chop off the surrounding parenthesis for the text
    nestedList[n] = [bindingSite, parse(whatIsBound[1:(len(whatIsBound) - 1)])]


  return(nestedList)



# Parses and prepares a nested list from given text
def getNestedListFromString(text, substrate):

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

        # Can input None because the 'substrate' was just set to the newProtein and does not need to be set
        createChildProteins(childInfo[1], None)


  # The right structure, but with strings instead of protein objects
  nestedList = parse(text)

  # Injects the protein objects into the structure
  createChildProteins(nestedList, substrate)

  # Return
  return nestedList


# Returns a collection of simple proteins off of a text
def simpleParse(text, substrate):

  collection = [substrate]

  # Create child proteins for X in "X, n(Y), n(Z)..." from the nested lists using recursion
  # parentLevel is a nested list that begins "[X, ..."
  # Updates the parentLevel list with protein objects replacing the strings
  def createSimpleChildProteins(parentLevel):
      
      for childInfo in parentLevel[1:]:

        # We are... "HERE, n(Y), n(Z)..."

        # Hardcoded UBIQUITIN_SEQ can be replaced by reading childInfo[1][0] and finding the right protein
        childProtein = simpleUbiquitin(parent=parentLevel[0], parentSite=childInfo[0])

        # Now that the child is paired to the parent, its site becomes off-limits
        try:
          parentLevel[0].attachmentSites.remove(childInfo[0])
        except:
          raise Exception("Problem with removing " + str(childInfo[0]) + " from " + str(parentLevel[0]) + " ! This is probably due to passing a substrate object not usable in the text's format")

        # Updates the nested list representation to include the unique protein model
        childInfo[1][0] = childProtein
        collection.append(childProtein)

        # Now we are... "X, n(HERE), n(Z)..."
        createSimpleChildProteins(childInfo[1])


  # The right structure, but with strings instead of protein objects
  nestedList = parse(text)

  # Give the substrate as a starting point
  nestedList[0] = substrate

  # Injects the SIMPLE protein objects into the structure assuming that nestedList[0] = substrate
  createSimpleChildProteins(nestedList)

  # Return
  return collection


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



# Returns a text version of a given collection
# The list argument should contain only proteins
simpleOutput = ""
def simpleUnparse(enclosingCollection, names = False):


  global simpleOutput
  simpleOutput = ""

  maxParent = enclosingCollection[0]

  def rec(simpleProtein):


    #try:

    string = str(simpleProtein.parentSite) + "(" + (getProteinName(enclosingCollection.index(simpleProtein)) if names else "x") + ","

    global simpleOutput
    simpleOutput += string

   # except:
    #  print("Protein is not part of the enclosing collection. This is normally caused by a child not part of the parent's collection")
    
    # Recurse for each child
    for child in filter(lambda p: p.parent == simpleProtein, enclosingCollection):
      rec(child)

    # Cap off with a )
    simpleOutput += "),"

  # Begin recursion
  rec(maxParent)

  # Format correctly... the output now looks something like -1(x,3(x,),4(x,),),
  simpleOutput = simpleOutput[3:-3]

  # The output now looks  like x,3(x,),4(x,)
  simpleOutput = simpleOutput.replace(",)", ")")

  # The output now looks  like x,3(x),4(x)... which was the goal!
  return simpleOutput


