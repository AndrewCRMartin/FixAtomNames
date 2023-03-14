/************************************************************************/
/**

   Program:    
   \file       FixAtomLabels.c
   
   \version    V1.0
   \date       13.03.23   
   \brief      Routines to fix symmetrical atom labels
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2023
   \author     Prof. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   Corrects SP1 hybridized symmetrical atoms (ARG, ASP, GLU, PHE, and TYR)
   such that the Cx1 has the torsion angle closer to zero and corrects
   SP2 hybridized symmetrical atoms (LEU, VAL) such that the smaller 
   angle going from Cx1 to Cx2 (which is ~120 degrees) is positive.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0    13.03.23   Original   By: ACRM

*************************************************************************/
/* Includes
*/
#include <math.h>

#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"
#include "FixAtomLabels.h"

/************************************************************************/
/* Defines and macros
*/
#define FAL_ERROR_VALUE 9999.0

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static BOOL NeedToSwapSP2Atoms(REAL tor1, REAL tor2);
static REAL AngleDistanceFromZero(REAL angle);
static void SwapAtomCoords(PDB *atom1, PDB *atom2);
static REAL CalcAngleDiff(REAL tor1, REAL tor2);
static REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians);

/************************************************************************/
void blFixAtomLabels(PDB *pdb, int verbose)
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
      else if(!strncmp(res->resnam, "VAL", 3))
      {
         atom[0] = blFindAtomInRes(res, "N   ");
         atom[1] = blFindAtomInRes(res, "CA  ");
         atom[2] = blFindAtomInRes(res, "CB  ");
         atom[3] = blFindAtomInRes(res, "CG1 ");
         atom[4] = blFindAtomInRes(res, "CG2 ");
      }
#ifdef ILE
      else if(!strncmp(res->resnam, "ILE", 3))
      {
         atom[0] = blFindAtomInRes(res, "N   ");
         atom[1] = blFindAtomInRes(res, "CA  ");
         atom[2] = blFindAtomInRes(res, "CB  ");
         /* Note this is the other way round! */
         atom[3] = blFindAtomInRes(res, "CG2 ");
         atom[4] = blFindAtomInRes(res, "CG1 ");
      }
#endif
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
#ifdef ILE
               !strncmp(res->resnam, "ILE", 3) ||
#endif
               !strncmp(res->resnam, "VAL", 3))
            {
            
               diff = CalcAngleDiff(tor1, tor2);
         
               if((diff < 90) || (diff > 180))
               {
                  if(verbose >= 1)
                  {
                     fprintf(stderr,"Swapped atom labels for %s %s%d%s\n",
                             atom[0]->resnam,
                             atom[0]->chain,
                             atom[0]->resnum,
                             atom[0]->insert);
                  }
                  SwapAtomCoords(atom[3], atom[4]);
               }
               else if(verbose >= 2)
               {
                  fprintf(stderr,"Atom labels for %s %s%d%s are OK\n",
                          atom[0]->resnam,
                          atom[0]->chain,
                          atom[0]->resnum,
                          atom[0]->insert);
               }
            }
            else /* PHE, TYR, ASP, GLU, ARG                             */
            {
               if(NeedToSwapSP2Atoms(tor1, tor2))
               {
                  if(verbose >= 1)
                  {
                     fprintf(stderr,"Swapped atom labels for %s %s%d%s\n",
                             atom[0]->resnam,
                             atom[0]->chain,
                             atom[0]->resnum,
                             atom[0]->insert);
                  }

                  SwapAtomCoords(atom[3], atom[4]);
                  
                  if(!strncmp(res->resnam, "PHE", 3) ||
                     !strncmp(res->resnam, "TYR", 3))
                  {
                     atom[5] = blFindAtomInRes(res, "CE1 ");
                     atom[6] = blFindAtomInRes(res, "CE2 ");
                     
                     SwapAtomCoords(atom[5], atom[6]);
                  }
               }
               else if(verbose >= 2)
               {
                  fprintf(stderr,"Atom labels for %s %s%d%s are OK\n",
                          atom[0]->resnam,
                          atom[0]->chain,
                          atom[0]->resnum,
                          atom[0]->insert);
               }
            }
         }
      }
   }
}

/************************************************************************/
void blPrintTorsionAtomLabels(FILE *out, PDB *pdb)
{
   PDB  *res     = NULL,
        *nextres = NULL;
   char resspec[16];
   
   
   for(res=pdb; res!=NULL; res=nextres)
   {
      nextres = blFindNextResidue(res);
      blBuildResSpec(res, resspec);
      
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
         
         fprintf(out, "%s %6s Tor1: %8.3f Tor2: %8.3f %s \
(Diff: %8.3f)\n",
                 res->resnam, resspec, tor1, tor2, label, diff);
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
         
         fprintf(out, "%s %6s Tor1: %8.3f Tor2: %8.3f %s \
(Diff: %8.3f)\n",
                 res->resnam, resspec, tor1, tor2, label, diff);

      }
#ifdef ILE
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
         
         fprintf(out, "%s %6s Tor1: %8.3f Tor2: %8.3f %s \
(Diff: %8.3f)\n",
                 res->resnam, resspec, tor1, tor2, label, diff);

      }
#endif
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
         
         fprintf(out, "%s %6s Tor1: %8.3f Tor2: %8.3f %s \n",
                 res->resnam, resspec, tor1, tor2, label);

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
         fprintf(out, "%s %6s Tor1: %8.3f Tor2: %8.3f %s \n",
                 res->resnam, resspec, tor1, tor2, label);
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
         
         fprintf(out, "%s %6s Tor1: %8.3f Tor2: %8.3f %s \n",
                 res->resnam, resspec, tor1, tor2, label);
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
         fprintf(out, "%s %6s Tor1: %8.3f Tor2: %8.3f %s \n",
                 res->resnam, resspec, tor1, tor2, label);
      }
   }
}

/************************************************************************/
static BOOL NeedToSwapSP2Atoms(REAL tor1, REAL tor2)
{
   if(AngleDistanceFromZero(tor1) > AngleDistanceFromZero(tor2))
      return(TRUE);
   return(FALSE);
}

/************************************************************************/
static REAL AngleDistanceFromZero(REAL angle)
{
   /* Make sure angles are -180...180 */
   if(angle > 180.0)
      angle -= 360.0;

   return(fabs(angle));
}


/************************************************************************/
static void SwapAtomCoords(PDB *atom1, PDB *atom2)
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



/************************************************************************/
static REAL CalcAngleDiff(REAL tor1, REAL tor2)
{
   REAL diff;
   diff = tor2 - tor1;
   while(diff < 0)
      diff += 360;
   while(diff > 360)                   
      diff -= 360;
   return(diff);
}



/************************************************************************/
static REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians)
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

