package require psfgen

resetpsf

topology minimization/scripts/top_all27_prot_lipid_na.inp

pdbalias residue HIS HSD
pdbalias residue CYX CYS
pdbalias atom ILE CD1 CD

segment A {

first GLYP
last CTER
pdb chainA.pdb

}
coordpdb chainA.pdb A
guesscoord

writepsf MBDLG2HB010.B2842V_autopsf.psf

writepdb MBDLG2HB010.B2842V_autopsf.pdb

