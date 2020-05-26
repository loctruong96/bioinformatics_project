# bioinformatics_project

#multiMultiMutant.sh

To use this script, place it in the multiMutant repo, in the same directory as the original multiMutant.sh
It takes the following parameters: pdbID, chainID startResidue:endResidue numMutations, -flags
Example call: ./multiMultiMutant.sh 1HHP A 40:45 2 -em

This will mutate residues 40-45 of 1HHP chain A, and perform energy mimimization. It will do this for all possible 2 point mutations within the given residue range.

In addition, this script makes parellel calls to promute.

