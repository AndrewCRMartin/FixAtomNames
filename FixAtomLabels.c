#include <math.h>

#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"
#include "FixAtomLabels.h"
void FixAtomLabels(PDB *pdb)
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
      if(!strncmp(res->resnam, "LEU", 3) ||
         !strncmp(res->resnam, "PHE", 3) ||
         !strncmp(res->resnam, "TYR", 3))
      {
         atom[0] = blFindAtomInRes(res, "CA  ");
         atom[1] = blFindAtomInRes(res, "CB  ");
         atom[2] = blFindAtomInRes(res, "CG  ");
         atom[3] = blFindAtomInRes(res, "CD1 ");
         atom[4] = blFindAtomInRes(res, "CD2 ");
      }
      else if(!strncmp(res->resnam, "ASP", 3))
      {
         atom[0] = blFindAtomInRes(res, "CA  ");
         atom[1] = blFindAtomInRes(res, "CB  ");
         atom[2] = blFindAtomInRes(res, "CG  ");
         atom[3] = blFindAtomInRes(res, "OD1 ");
         atom[4] = blFindAtomInRes(res, "OD2 ");
      }
      else if(!strncmp(res->resnam, "GLU", 3))
      {
         atom[0] = blFindAtomInRes(res, "CB  ");
         atom[1] = blFindAtomInRes(res, "CG  ");
         atom[2] = blFindAtomInRes(res, "CD  ");
         atom[3] = blFindAtomInRes(res, "OE1 ");
         atom[4] = blFindAtomInRes(res, "OE2 ");
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
      else if(!strncmp(res->resnam, "ARG", 3))
      {
         atom[0] = blFindAtomInRes(res, "CD  ");
         atom[1] = blFindAtomInRes(res, "NE  ");
         atom[2] = blFindAtomInRes(res, "CZ  ");
         atom[3] = blFindAtomInRes(res, "NH1 ");
         atom[4] = blFindAtomInRes(res, "NH2 ");
      }
      
         
      if(atom[0] != NULL)
      {
         tor1 = CalcTorsion(atom[0],atom[1],atom[2],atom[3],FALSE);
         tor2 = CalcTorsion(atom[0],atom[1],atom[2],atom[4],FALSE);

         if((tor1 < FAL_ERROR_VALUE-1) && (tor2 < FAL_ERROR_VALUE-1))
         {
            if(!strncmp(res->resnam, "LEU", 3) ||
               !strncmp(res->resnam, "ILE", 3) ||
               !strncmp(res->resnam, "VAL", 3))
            {
            
               diff = CalcAngleDiff(tor1, tor2);
         
               if((diff < 90) || (diff > 180))
               {
                  SwapAtomCoords(atom[3], atom[4]);
                  
                  if(!strncmp(res->resnam, "PHE", 3) ||
                     !strncmp(res->resnam, "TYR", 3))
                  {
                     atom[5] = blFindAtomInRes(res, "CE1 ");
                     atom[6] = blFindAtomInRes(res, "CE2 ");
                     
                     SwapAtomCoords(atom[5], atom[6]);
                  }
               }
            }
            else
            {
               if(NeedToSwapSP2Atoms(tor1, tor2))
               {
                  SwapAtomCoords(atom[3], atom[4]);
                  
                  if(!strncmp(res->resnam, "PHE", 3) ||
                     !strncmp(res->resnam, "TYR", 3))
                  {
                     atom[5] = blFindAtomInRes(res, "CE1 ");
                     atom[6] = blFindAtomInRes(res, "CE2 ");
                     
                     SwapAtomCoords(atom[5], atom[6]);
                  }
               }
            }
         }
      }
   }
}

BOOL NeedToSwapSP2Atoms(REAL tor1, REAL tor2)
{
   if(AngleDistanceFromZero(tor1) > AngleDistanceFromZero(tor2))
      return(TRUE);
   return(FALSE);
}

REAL AngleDistanceFromZero(REAL angle)
{
   /* Make sure angles are -180...180 */
   if(angle > 180.0)
      angle -= 360.0;

   return(fabs(angle));
}


void SwapAtomCoords(PDB *atom1, PDB *atom2)
{
   PDB tmp;

   tmp = *atom1;
   
   atom1->x = atom2->x;
   atom1->y = atom2->y;
   atom1->z = atom2->z;

   atom2->x = tmp.x;
   atom2->y = tmp.y;
   atom2->z = tmp.z;
}



REAL CalcAngleDiff(REAL tor1, REAL tor2)
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
      return(FAL_ERROR_VALUE);
   }

   tor = blPhi(p1->x, p1->y, p1->z,
               p2->x, p2->y, p2->z,
               p3->x, p3->y, p3->z,
               p4->x, p4->y, p4->z);

   if(!Radians)
      tor *= 180/PI;

   return(tor);
}

