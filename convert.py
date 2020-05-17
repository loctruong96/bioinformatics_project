def letterToThree(oneLetter):
    oneLetter = oneLetter.upper()
    letters = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "B": "ASX",
               "C": "CYS", "E": "GLU", "Q": "GLN", "Z": "GLX", "G": "GLY",
               "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET",
               "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
               "Y": "TYR", "V": "VAL"}
    retVal = letters.get(oneLetter, 'N/A')
    return retVal

def threeToLetter(threeLetter):
    threeLetter = threeLetter.upper()
    letters = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASX": "B",
               "CYS": "C", "GLU": "E", "GLN": "Q", "GLX": "Z", "GLY": "G",
               "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M",
               "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
               "TYR": "Y", "VAL": "V"}
    retVal = letters.get(oneLetter, 'N/A')
    return retVal


