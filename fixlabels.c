#include <stdio.h>
#include <math.h>

#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"
#include "FixAtomLabels.h"

int main(int argc, char **argv);
void DoAnalysis(PDB *pdb);



int main(int argc, char **argv)
{
   FILE *in;
   
   if((in=fopen(argv[1], "r"))!=NULL)
   {
      PDB *pdb;
      int natoms;
      
      if((pdb = blReadPDB(in, &natoms))!=NULL)
      {
         printf("Before...\n");
         DoAnalysis(pdb);
         FixAtomLabels(pdb);
         printf("After...\n");
         DoAnalysis(pdb);
         blWritePDB(stdout, pdb);
         FREELIST(pdb, PDB);
      }
      else
      {
         fprintf(stderr, "No atoms read from file %s\n", argv[1]);
         return(1);
      }
      fclose(in);
   }
   else
   {
      fprintf(stderr, "Can't read file %s\n", argv[1]);
      return(1);
   }
   return(0);
}

void DoAnalysis(PDB *pdb)
{
   PDB *res     = NULL,
       *nextres = NULL;

   for(res=pdb; res!=NULL; res=nextres)
   {
      nextres = blFindNextResidue(res);
      if(!strncmp(res->resnam, "LEU", 3))
      {
         PDB *CA, *CB, *CG, *CD1, *CD2;
         REAL tor1, tor2, diff;
         char label[16];
            
         CA  = blFindAtomInRes(res, "CA  ");
         CB  = blFindAtomInRes(res, "CB  ");
         CG  = blFindAtomInRes(res, "CG  ");
         CD1 = blFindAtomInRes(res, "CD1 ");
         CD2 = blFindAtomInRes(res, "CD2 ");
         
         tor1 = CalcTorsion(CA, CB, CG, CD1, FALSE);
         tor2 = CalcTorsion(CA, CB, CG, CD2, FALSE);

         diff = CalcAngleDiff(tor1, tor2);
         strcpy(label, "OK");
         
         if((diff < 90) || (diff > 180))
            strcpy(label, "SWAPPED!");
         
         printf("%s Tor1: %8.3f Tor2: %8.3f Diff: %8.3f %s\n",
                res->resnam, tor1, tor2, diff, label);

      }
      else if(!strncmp(res->resnam, "VAL", 3))
      {
         PDB *N, *CA, *CB, *CG1, *CG2;
         REAL tor1, tor2, diff;
         char label[16];
            
         N   = blFindAtomInRes(res, "N   ");
         CA  = blFindAtomInRes(res, "CA  ");
         CB  = blFindAtomInRes(res, "CB  ");
         CG1 = blFindAtomInRes(res, "CG1 ");
         CG2 = blFindAtomInRes(res, "CG2 ");
         
         tor1 = CalcTorsion(N, CA, CB, CG1, FALSE);
         tor2 = CalcTorsion(N, CA, CB, CG2, FALSE);

         diff = CalcAngleDiff(tor1, tor2);
         strcpy(label, "OK");
         
         if((diff < 90) || (diff > 180))
            strcpy(label, "SWAPPED!");
         
         printf("VAL  Tor1: %8.3f Tor2: %8.3f Diff: %8.3f %s\n",
                tor1, tor2, diff, label);

      }
      else if(!strncmp(res->resnam, "ILE", 3))
      {
         PDB  *N, *CA, *CB, *CG1, *CG2;
         REAL tor1, tor2, diff;
         char label[16];
            
         N   = blFindAtomInRes(res, "N   ");
         CA  = blFindAtomInRes(res, "CA  ");
         CB  = blFindAtomInRes(res, "CB  ");
         CG1 = blFindAtomInRes(res, "CG1 ");
         CG2 = blFindAtomInRes(res, "CG2 ");

         /* Note this is the other way round! */
         tor1 = CalcTorsion(N, CA, CB, CG2, FALSE);
         tor2 = CalcTorsion(N, CA, CB, CG1, FALSE);

         diff = CalcAngleDiff(tor1, tor2);
         strcpy(label, "OK");
         
         if((diff < 90) || (diff > 180))
            strcpy(label, "SWAPPED!");
         
         printf("ILE  Tor1: %8.3f Tor2: %8.3f Diff: %8.3f %s\n",
                tor1, tor2, diff, label);

      }
      else if(!strncmp(res->resnam, "PHE", 3) ||
              !strncmp(res->resnam, "TYR", 3))
      {
         PDB *CA, *CB, *CG, *CD1, *CD2;
         REAL tor1, tor2;
         char label[16];
            
         CA  = blFindAtomInRes(res, "CA  ");
         CB  = blFindAtomInRes(res, "CB  ");
         CG  = blFindAtomInRes(res, "CG  ");
         CD1 = blFindAtomInRes(res, "CD1 ");
         CD2 = blFindAtomInRes(res, "CD2 ");
         
         tor1 = CalcTorsion(CA, CB, CG, CD1, FALSE);
         tor2 = CalcTorsion(CA, CB, CG, CD2, FALSE);

         strcpy(label, "OK");
         
         if(NeedToSwapSP2Atoms(tor1, tor2))
            strcpy(label, "SWAPPED!");
         
         printf("%s Tor1: %8.3f Tor2: %8.3f %s\n",
                res->resnam, tor1, tor2, label);

      }
      else if(!strncmp(res->resnam, "ASP", 3))
      {
         PDB *CA, *CB, *CG, *OD1, *OD2;
         REAL tor1, tor2;
         char label[16];
            
         CA  = blFindAtomInRes(res, "CA  ");
         CB  = blFindAtomInRes(res, "CB  ");
         CG  = blFindAtomInRes(res, "CG  ");
         OD1 = blFindAtomInRes(res, "OD1 ");
         OD2 = blFindAtomInRes(res, "OD2 ");
         
         tor1 = CalcTorsion(CA, CB, CG, OD1, FALSE);
         tor2 = CalcTorsion(CA, CB, CG, OD2, FALSE);

         strcpy(label, "OK");
         if(NeedToSwapSP2Atoms(tor1, tor2))
         {
            strcpy(label, "SWAPPED!");
         }
         printf("%s Tor1: %8.3f Tor2: %8.3f %s\n",
                res->resnam, tor1, tor2, label);
      }
      else if(!strncmp(res->resnam, "GLU", 3))
      {
         PDB  *CB, *CG, *CD, *OE1, *OE2;
         REAL tor1, tor2;
         char label[16];
            
         CB  = blFindAtomInRes(res, "CB  ");
         CG  = blFindAtomInRes(res, "CG  ");
         CD  = blFindAtomInRes(res, "CD  ");
         OE1 = blFindAtomInRes(res, "OE1 ");
         OE2 = blFindAtomInRes(res, "OE2 ");
         
         tor1 = CalcTorsion(CB, CG, CD, OE1, FALSE);
         tor2 = CalcTorsion(CB, CG, CD, OE2, FALSE);

         strcpy(label, "OK");
         if(NeedToSwapSP2Atoms(tor1, tor2))
         {
            strcpy(label, "SWAPPED!");
         }
         
         printf("%s Tor1: %8.3f Tor2: %8.3f %s\n",
                res->resnam, tor1, tor2, label);
      }
      else if(!strncmp(res->resnam, "ARG", 3))
      {
         PDB  *CD, *NE, *CZ, *NH1, *NH2;
         REAL tor1, tor2;
         char label[16];
            
         CD  = blFindAtomInRes(res, "CD  ");
         NE  = blFindAtomInRes(res, "NE  ");
         CZ  = blFindAtomInRes(res, "CZ  ");
         NH1 = blFindAtomInRes(res, "NH1 ");
         NH2 = blFindAtomInRes(res, "NH2 ");

         tor1 = CalcTorsion(CD, NE, CZ, NH1, FALSE);
         tor2 = CalcTorsion(CD, NE, CZ, NH2, FALSE);

         strcpy(label, "OK");
         if(NeedToSwapSP2Atoms(tor1, tor2))
         {
            strcpy(label, "SWAPPED!");
         }
         printf("%s Tor1: %8.3f Tor2: %8.3f %s\n",
                res->resnam, tor1, tor2, label);
      }
   }
}

