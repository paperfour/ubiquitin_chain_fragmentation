from iontype import IonType

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

# ============================= Utils for protein identification =================================

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


# Returns the parent of the parent of... a given protein
# Goes down the tree until the selected node has no parent itself
def getBaseProtein(protein):

  if(protein.parentAttachmentSite == -1):
    return protein
  else:
    return getBaseProtein(protein.parent)
