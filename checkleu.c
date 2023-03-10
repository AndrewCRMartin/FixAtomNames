#include <stdio.h>
#include <math.h>

#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"

#define ERROR_VALUE 9999.0

int main(int argc, char **argv);
void DoAnalysis(PDB *pdb);
REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians);
REAL CalcDiff(REAL tor1, REAL tor2);
void FixLeu(PDB *pdb);
void SwapAtomCoords(PDB *atom1, PDB *atom2);






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
         FixLeu(pdb);
         printf("After...\n");
         DoAnalysis(pdb);
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

         diff = CalcDiff(tor1, tor2);
         strcpy(label, "OK");
         
         if((diff < 90) || (diff > 180))
            strcpy(label, "SWAPPED!");
         
         printf("LEU Tor1: %8.3f Tor2: %8.3f Diff: %8.3f %s\n",
                tor1, tor2, diff, label);

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

         diff = CalcDiff(tor1, tor2);
         strcpy(label, "OK");
         
         if((diff < 90) || (diff > 180))
            strcpy(label, "SWAPPED!");
         
         printf("VAL Tor1: %8.3f Tor2: %8.3f Diff: %8.3f %s\n",
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

         diff = CalcDiff(tor1, tor2);
         strcpy(label, "OK");
         
         if((diff < 90) || (diff > 180))
            strcpy(label, "SWAPPED!");
         
         printf("ILE Tor1: %8.3f Tor2: %8.3f Diff: %8.3f %s\n",
                tor1, tor2, diff, label);

      }
   }
}

void FixLeu(PDB *pdb)
{
   PDB *res     = NULL,
       *nextres = NULL;

   for(res=pdb; res!=NULL; res=nextres)
   {
      PDB *atom[8];
      REAL tor1, tor2, diff;
      int i;
      for(i=0; i<8; i++)
         atom[i] = NULL;
            
      nextres = blFindNextResidue(res);
      if(!strncmp(res->resnam, "LEU", 3))
      {
         atom[0] = blFindAtomInRes(res, "CA  ");
         atom[1] = blFindAtomInRes(res, "CB  ");
         atom[2] = blFindAtomInRes(res, "CG  ");
         atom[3] = blFindAtomInRes(res, "CD1 ");
         atom[4] = blFindAtomInRes(res, "CD2 ");
      }
      else if(!strncmp(res->resnam, "VAL", 3))
      {
         atom[0] = blFindAtomInRes(res, "N   ");
         atom[1] = blFindAtomInRes(res, "CA  ");
         atom[2] = blFindAtomInRes(res, "CB  ");
         atom[3] = blFindAtomInRes(res, "CG1 ");
         atom[4] = blFindAtomInRes(res, "CG2 ");
      }
      else if(!strncmp(res->resnam, "ILE", 3))
      {
         atom[0] = blFindAtomInRes(res, "N   ");
         atom[1] = blFindAtomInRes(res, "CA  ");
         atom[2] = blFindAtomInRes(res, "CB  ");
         /* Note this is the other way round! */
         atom[3] = blFindAtomInRes(res, "CG2 ");
         atom[4] = blFindAtomInRes(res, "CG1 ");
      }
      
         
      if(atom[0] != NULL)
      {
         tor1 = CalcTorsion(atom[0],atom[1],atom[2],atom[3],FALSE);
         tor2 = CalcTorsion(atom[0],atom[1],atom[2],atom[4],FALSE);

         if((tor1 < 9998.0) && (tor2 < 9998.0))
         {
            diff = CalcDiff(tor1, tor2);
         
            if((diff < 90) || (diff > 180))
            {
               SwapAtomCoords(atom[3], atom[4]);
            }
         }
         
      }
   }
}


void SwapAtomCoords(PDB *atom1, PDB *atom2)
{
   PDB tmp;
   REAL x,y,z;
   
   fprintf(stderr,"Before: atom1 %.3f %.3f %.3f\n", atom1->x, atom1->y, atom1->z);
   fprintf(stderr,"        atom2 %.3f %.3f %.3f\n", atom2->x, atom2->y, atom2->z);
   
   
   x = atom1->x;
   y = atom1->y;
   z = atom1->z;
   
   atom1->x = atom2->x;
   atom1->y = atom2->y;
   atom1->z = atom2->z;
   
   atom2->x = x;
   atom2->y = y;
   atom2->x = z;

   fprintf(stderr,"After:  atom1 %.3f %.3f %.3f\n",   atom1->x, atom1->y, atom1->z);
   fprintf(stderr,"        atom2 %.3f %.3f %.3f\n\n", atom2->x, atom2->y, atom2->z);
}



REAL CalcDiffOld(REAL tor1, REAL tor2)
{
   REAL diff;
   diff = tor2 - tor1;
   while(diff < 0)
      diff += 180.0;
   return(diff);
}

REAL CalcDiffOld2(REAL tor1, REAL tor2)
{
   REAL diff;
   if(tor1 < 0)
      tor1 += 360;
   if(tor2 < 0)
      tor2 += 360;
   
   diff = tor2 - tor1;
/*
   if(diff < -180.0)
      diff += 360.0;
*/ return(diff);
}


REAL CalcDiff(REAL tor1, REAL tor2)
{
   REAL diff;
   diff = tor2 - tor1;
   while(diff < 0)
      diff += 360;
   while(diff > 360)                   
      diff -= 360;
   return(diff);
}



REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians)
{
   REAL tor;
   
   if((p1==NULL)||(p2==NULL)||(p3==NULL)||(p4==NULL))
   {
      return(ERROR_VALUE);
   }

   tor = blPhi(p1->x, p1->y, p1->z,
               p2->x, p2->y, p2->z,
               p3->x, p3->y, p3->z,
               p4->x, p4->y, p4->z);

   if(!Radians)
      tor *= 180/PI;

   return(tor);
}

