#ifndef _FixAtomLabels_h_
#define _FixAtomLabels_h_ 1

void FixAtomLabels(PDB *pdb);
BOOL NeedToSwapSP2Atoms(REAL tor1, REAL tor2);
REAL AngleDistanceFromZero(REAL angle);
void SwapAtomCoords(PDB *atom1, PDB *atom2);
REAL CalcAngleDiff(REAL tor1, REAL tor2);
REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians);
#define FAL_ERROR_VALUE 9999.0

#endif
