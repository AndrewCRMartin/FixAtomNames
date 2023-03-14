/************************************************************************/
/**

   \file       pdbflip.c
   
   \version    V2.0
   \date       13.03.23
   \brief      Standardise equivalent atom labelling
   
   \copyright  (c) UCL, Prof. Andrew C. R. Martin 1996-2023
   \author     Prof. Andrew C. R. Martin
   \par
               Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This code is NOT IN THE PUBLIC DOMAIN, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC.

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   
-  V1.1   22.07.14 Renamed deprecated functions with bl prefix.
                   Added doxygen annotation. By: CTP
-  V1.2   19.08.14 Removed unused variable in DoFlipping() By: CTP
-  V1.3   06.11.14 Renamed from flip
-  V1.4   13.02.15 Added whole PDB support
-  V1.5   12.03.15 Changed to allow multi-character chain names
-  V2.0   13.03.23 Complete rewrite for new flipping code which is now
                   in BiopLib

*************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "bioplib/SysDefs.h"
#include "bioplib/general.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"
#include "FixAtomLabels.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define FAL_ERROR_VALUE 9999.0

/************************************************************************/
/* Globals
*/


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *verbosity, BOOL *reportOnly);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program to flip sidechain equivalent atom naming

-  08.11.96 Original   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  13.02.15 Added whole PDB support.  By: ACRM
-  13.03.23 Complete rewrite to use new flip routines
*/
int main(int argc, char **argv)
{
   FILE     *in        = stdin,
            *out       = stdout;
   int      verbosity  = 0;
   BOOL     reportOnly = FALSE;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];
   WHOLEPDB *wpdb;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &verbosity, &reportOnly))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((wpdb = blReadWholePDB(in)) != NULL)
         {
            PDB *pdb;
            pdb = wpdb->pdb;
            if(reportOnly)
            {
               blPrintTorsionAtomLabels(out, pdb);
            }
            else
            {
               blFixAtomLabels(pdb, verbosity);
               blWriteWholePDB(out, wpdb);
            }
            FREELIST(pdb, PDB);
         }
         else
         {
            fprintf(stderr,"No atoms read from PDB file\n");
         }
      }
   }
   else
   {
      Usage();
   }

   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     int *verbosity, BOOL *reportOnly)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *verbosity   Information level
   \param[out]     *reportOnly  Report wrong residues rather than fixing
   \return                      Success?

   Parse the command line
   
-  08.11.96 Original    By: ACRM
-  13.02.23 Updated for V2.0
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *verbosity, BOOL *reportOnly)
{
   argc--;
   argv++;

   infile[0]   = outfile[0] = '\0';
   *verbosity  = 0;
   *reportOnly = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'v':
            *verbosity = 1;
            if(argv[0][2] == 'v')
               *verbosity=2;
            break;
         case 'r':
            *reportOnly = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  08.11.96 Original   By: ACRM
-  22.07.14 V1.1 By: CTP
-  06.11.14 V1.2 By: ACRM
-  12.03.15 V1.5
-  13.03.23 V2.0
*/
void Usage(void)
{
   fprintf(stderr,"\npdbflip V2.0 (c) 2014-2023 Prof. Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"\nUsage: pdbflip [-v[v]] [-r] [in.pdb [out.pdb]]\n");
   fprintf(stderr,"               -v   Report fixed atoms\n");
   fprintf(stderr,"               -vv  Report unfixed atoms as well\n");
   fprintf(stderr,"               -r   Only report atoms rather than \
fixing\n");

   fprintf(stderr,"\npdbflip V2 is a much-improved program for fixing \
the names of\n");
   printf("symmetrical atoms in ARG, ASP, GLU, PHE, TYR, LEU and VAL.\n");
   printf("This new version simply swaps the coordinates rather than \
trying to worry\n");
   printf("explicitly about connectivity, etc., making the code much \
simpler. It \n");
   printf("still assumes that the connectivity is correct (e.g. in PHE, \
CE1 is\n");
   printf("connected to CD1 and CE2 is connected to CD2).\n");
   printf("\nLEU and VAL were not handled by the old version.\n\n");
}


