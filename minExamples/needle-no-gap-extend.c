/*This is a minimal example of how to run alnSeq's Needlman
` alignment without gap extension penalties
*/

#include "../needleman/needleNoGap.h"
#include "../general/alnStruct.h"

int main(
){ /*main*/
   uchar errUC = 0; /*Reports errors from functions*/

   char *outRefSeqStr = 0; /*Holds the aligned sequence*/
   char *outQrySeqStr = 0; /*Holds the aligned sequence*/
   long scoreL = 0;        /*Holds the final score*/

   struct seqStruct refST;
   struct seqStruct qryST;

   /*alignment settings*/   
   struct alnSet settings;

   /*For holding alignments*/
   struct alnStruct *alnST = 0;
   struct alnMatrix *matrixST = 0;

   FILE *refFILE =
      fopen("../analysis/genomes/Small-ref.fasta", "r");
   FILE *qryFILE = 
      fopen("../analysis/genomes/Small-query.fasta", "r");

   /*Initialize structures*/
   initSeqST(&refST);
   initSeqST(&qryST);
   initAlnSet(&settings);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Read in the reference and query sequences
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
   /*Check If I opened the files*/
   if(refFILE == 0 || qryFILE == 0)
   { /*If: I could not open the reference file*/
      if(refFILE != 0) fclose(refFILE);
      if(qryFILE != 0) fclose(qryFILE);
      return -1;
   } /*If: I could not open the reference file*/

   errUC = readFaSeq(refFILE, &refST);
   fclose(refFILE);
   refFILE = 0;

   if(errUC != 1)
   { /*If: I had an error*/
      fclose(qryFILE);
      return -1;
   } /*If: I had an error*/

   errUC = readFaSeq(qryFILE, &qryST);
   fclose(qryFILE);
   qryFILE = 0;

   if(errUC != 1)
   { /*If: I had an error reading the query sequence*/
      freeSeqSTStack(&refST);
      return -1;
   } /*If: I had an error reading the query sequence*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Set up for alignment and do alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Set last base to align (-1 to convert index 1 to 0).
   ` The default is 0, so it needs to be set
   */
   refST.endAlnUL = refST.lenSeqUL - 1;
   qryST.endAlnUL = qryST.lenSeqUL - 1;

   /*Set first base to align (as index 0). This can be
   ` skipped if you are just doing 0 (is default).
   */
   refST.offsetUL = 0;
   qryST.offsetUL = 0;

   /*Convert sequences to lookup index's. This step
   ` speeds up the alignment step. It does nothing when
   ` -DNOSEQCNVT is used during the compile step and so,
   ` can be skipped with -DNOSEQCNVT
   */
   seqToLookupIndex(refST.seqCStr);
   seqToLookupIndex(qryST.seqCStr);

   /*Do alignments. This is for the non-two-bit array
   ` Needlemans and Watermans
   */
   matrixST = NeedleAlnNoGap(&qryST, &refST, &settings);

   /*Convert sequences back to bases*/
   lookupIndexToSeq(refST.seqCStr);
   lookupIndexToSeq(qryST.seqCStr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Get alignment and print out
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*These next steps seem to add extra effort, however,
   ` they are here so I can use the same print functions
   ` on my programs. It also alows me to extract multiple
   ` alignments when doing query reference scans. I just
   ` need to provide a separate index.
   */
   /*Convert the directional matrix to an alignment*/
   alnST =
      dirMatrixToAln(
         &refST,
         &qryST,
         matrixST->bestEndIndexUL,
         &settings,
         matrixST
   );

   scoreL = matrixST->bestScoreL;
   freeAlnMatrix(matrixST);
   matrixST = 0;

   /*Get the aligned sequences. alnSTToSeq and fprint can
   ` be replaced with prntAln(). This also removes the
   ` need for the outRefSeqStr and outQrySeqStr c-strings.
   */
   errUC =
      alnSTToSeq(
         &refST,
         &qryST,
         alnST,
         &settings,
         &outRefSeqStr,
         &outQrySeqStr
   );

   printf(
      "Score: %li\n%s%s\n%s%s\n",
      scoreL,
      refST.idCStr,
      outRefSeqStr,
      qryST.idCStr,
      outQrySeqStr
   );
   /*The refST.idCStr and qryST.idCStr still have the
   ` starting ">" and ending newline ("\n")
   */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Clean up and exit
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   freeAlnST(alnST);
   freeSeqSTStack(&refST);
   freeSeqSTStack(&qryST);
   free(outRefSeqStr);
   free(outQrySeqStr);

   return 0;
} /*main*/