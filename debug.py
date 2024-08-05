def getMass(formula: dict):

  carbon = 12.0107 * formula["C"]
  hydrogen = 1.00794 * formula["H"]
  nitrogen = 14.0067 * formula["N"]
  oxygen =  15.9994 * formula["O"]
  sulfur = 32.065 * formula["S"]

  return carbon + hydrogen + nitrogen + oxygen + sulfur


t = { "C": 6, "H": 12, "N": 4, "O": 1, "S": 0}


print(getMass(t))