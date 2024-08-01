from protein import *

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


# Returns the C fragment as a collection when cleaved at a site
def cIonCleave(referenceCollection, referenceCleavedProtein, cleavageSite: int):

  proteinCollection = copyProteinCollection(referenceCollection)
  cleavedProtein = proteinCollection[referenceCollection.index(referenceCleavedProtein)]

  # Ensure that the protein is long enough to be cleaved at the requested site
  if ((cleavageSite >= cleavedProtein.sequenceLength()) or (cleavageSite < 1)):
    raise Exception("Cannot cleave at " + str(cleavageSite))

  # Create the C-n fragment; this fragment will be parentless by nature but should still retain children left of the cleavage site
  nFragment = Protein(cleavedProtein.sequence[:(cleavageSite)], ionType = IonType.C)

  # Replace the old protein with the cleaved one in the protein collection
  proteinCollection.append(nFragment)

  # Reassign the old protein's children to the correct ions
  for attachmentSite in cleavedProtein.children.keys():

    child = cleavedProtein.children[attachmentSite]

    # If left of cleavage site
    if (attachmentSite <= cleavageSite):
      # Attach to N fragment
      child.setParent(nFragment, attachmentSite)
    
    
  # Remove any protein that is not connected to the nFragment
  for p in proteinCollection:
    if getBaseProtein(p) != nFragment:
        proteinCollection.remove(p)


  return proteinCollection
