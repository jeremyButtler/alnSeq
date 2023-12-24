/*#########################################################
# Name: alnSeqDefaults.h
# Use:
#  o Holds the default settings and global definitions for
#    find alignSeq
# Libraries:
# C Standard Libraries:
#########################################################*/

#ifndef ALNSEQDEFAULTS_H
#define ALNSEQDEFAULTS_H

/*Gereral find co-infections settings*/
#define defVersion 20231224  /*Version number for alnSeq*/

/*Aligment method being used*/
#define defUseNeedle 1  /*do Needleman-Wunsch alignment*/
#define defUseWater 0   /*Do a waterman smith alignment*/
#define defUseHirsch 0  /*Use a Hirschberg alignment*/
#define defUseMemWater 0 /*Use a memory efficent water*/

#define defUseTwoBit 0 /*1:use two bit arrays; 0 do not*/
   /*This only applies to Needleman and Waterman*/

#define defJustScoresBl 0 /*memory efficent water only*/
   /*Only print out coordiantes and scores for memWater*/

/*Output file formating (only change defFormat)*/
#define defExpandCig 0 /*Expanded cigar format flag*/
#define defEMBOSS 1    /*EMBOSS format flag*/
#define defClustal 2   /*Clustal format flag*/
#define defFasta 3     /*Fasta format flag*/
#define defFormat defExpandCig
   /*Change defFormat to change the default output
   ` The current default format is expand cigar
   */

#define defPAln 1
   /* 1: Include softmasked regions in output alignment
   `  0: Only print out the aligned regions
   */
#define defPPos 0
  /* 1: Print out the starting and ending position for
  `     each sequence line in the alignment file
  `  0: Do not print positions
  */
#define defLineWrap 59 /*Max chars per line in out file*/
   /*Minimun is 40 characters*/


/*Query reference scan settings*/
#define defQueryRefScan 0
   /* 1: Recored the best score for each reference and
   `     query base in in a Waterman Smith alnigment
   */
#define defMinScore 1000
    /*Min score to keep an non-best alignment when doing
    ` an mutli-entry Waterman Smith alignment
    */

#define defFilterByQuery 0
  /*1: Remove alternative alignments with overlaping query
  `    scores
  ` 0: Do not remove alternative alignments
  */
#define defFilterByRef 0
  /*1: Remove alternative alignments with overlaping
  `    reference scores
  ` 0: Do not remove alternative alignments
  */
#define defFilterByQueryRef 1
  /*1: Remove alternative alignments with overlaps in both
  ` the query and reference scores
  ` 0: Do not remove alternative alignments
  */
#define defPAltAln 0
  /*1: Print out the alignment for each alternative
  `    alignment (soft masked regions are not printed)
  `0: Print the scores and positions for the alternative
  `   alignments
  */

/*Alignment matrix movements*/
/*Do not change these values*/
#define defMvStop 0    /*Stop*/
#define defMvDel 1     /*Move left (deletion)*/
#define defMvIns 2     /*Move up (insertion)*/
#define defMvSnp 3     /*Move on a diagnol (snp/match)*/

/*Matrix filling settings*/

/* Settings for which direction to keep when their are
`  multiple equal directions.
   Order is AlwaysKeep->secondBest->lastOption
*/
#define defSnpInsDel 0
#define defSnpDelIns 1
#define defInsSnpDel 2
#define defInsDelSnp 3
#define defDelSnpIns 4
#define defDelInsSnp 5

#define defBestDir defInsDelSnp

/*Scoring variables*/
#define defNoGapExtend 0 /*1: do not use gap extension*/
#define defGapOpen -10  /*Penalty for starting indel*/
#define defGapExtend -1 /*Penalty for extading an indel*/
/* Scoring matrix is EDNAFULL or something close to it. I
`   think I have the anonymous bases with the correct
`  scores.
*/

/*Scoring matrix name*/
static char defMatrixNameStr[] = "EDNAFULL";

#define defWordMatchScore 5
#define defWordSNPScore -4 

/*Alignment scoring matrix (non-anonymous)*/
#define defAToA 5     /*Score for query A to ref A*/
#define defAToT -4    /*Score for query A to ref T*/
#define defAToU -4    /*Score for query A to ref T*/
#define defAToG -4    /*Score for query A to ref G*/
#define defAToC -4    /*Score for query A to ref C*/

#define defTToA -4    /*Score for query T to ref A*/
#define defTToT 5     /*Score for query T to ref T*/
#define defTToU 5     /*Score for query T to ref T*/
#define defTToG -4    /*Score for query T to ref G*/
#define defTToC -4    /*Score for query T to ref C*/

#define defGToA -4    /*Score for query G to ref A*/
#define defGToT -4    /*Score for query G to ref T*/
#define defGToU -4    /*Score for query G to ref T*/
#define defGToG 5     /*Score for query G to ref G*/
#define defGToC -4    /*Score for query G to ref C*/

#define defCToA -4    /*Score for query C to ref A*/
#define defCToT -4    /*Score for query C to ref T*/
#define defCToU -4    /*Score for query C to ref T*/
#define defCToG -4    /*Score for query C to ref G*/
#define defCToC 5     /*Score for query C to ref C*/

/*Alignment scoring matrix anonymous (A)*/
#define defAToW 1     /*Score for query A to ref W (AT)*/
#define defAToS -4    /*Score for query A to ref S (CG)*/
#define defAToM 1     /*Score for query A to ref M (AC)*/
#define defAToK -4    /*Score for query A to ref K (GT)*/
#define defAToR 1     /*Score for query A to ref R (AG)*/
#define defAToY -4    /*Score for query A to ref Y (CT)*/
#define defAToB -4    /*Score for query A to ref B (CGT)*/
#define defAToD -1    /*Score for query A to ref D (AGT)*/
#define defAToH -1    /*Score for query A to ref H (ACT)*/
#define defAToV -1    /*Score for query A to ref V (ACG)*/
#define defAToN -2    /*Score for query A to ref N (ACGT)*/
#define defAToX -2    /*Score for query A to ref X (ACGT)*/

#define defWToA 1     /*Score for query W to ref A (AT)*/
#define defSToA -4    /*Score for query S to ref A (CG)*/
#define defMToA 1     /*Score for query M to ref A (AC)*/
#define defKToA -4    /*Score for query K to ref A (GT)*/
#define defRToA 1     /*Score for query R to ref A (AG)*/
#define defYToA -4    /*Score for query Y to ref A (CT)*/
#define defBToA -4    /*Score for query B to ref A (CGT)*/
#define defDToA -1    /*Score for query D to ref A (AGT)*/
#define defHToA -1    /*Score for query H to ref A (ACT)*/
#define defVToA -1    /*Score for query V to ref A (ACG)*/
#define defNToA -2    /*Score for query N to ref A (ACGT)*/
#define defXToA -2    /*Score for query X to ref A (ACGT)*/

/*Alignment scoring matrix anonymous (C)*/
#define defCToW -4    /*Score for query C to ref W (AT)*/
#define defCToS 1     /*Score for query C to ref S (CG)*/
#define defCToM 1     /*Score for query C to ref M (AC)*/
#define defCToK -4    /*Score for query C to ref K (GT)*/
#define defCToR -4    /*Score for query C to ref R (AG)*/
#define defCToY 1     /*Score for query C to ref Y (CT)*/
#define defCToB -1    /*Score for query C to ref B (CGT)*/
#define defCToD -4    /*Score for query C to ref D (AGT)*/
#define defCToH -1    /*Score for query C to ref H (ACT)*/
#define defCToV -1    /*Score for query C to ref V (ACG)*/
#define defCToN -2    /*Score for query C to ref N (ACGT)*/
#define defCToX -2    /*Score for query C to ref X (ACGT)*/

#define defWToC -4    /*Score for query W to ref C (AT)*/
#define defSToC 1     /*Score for query S to ref C (CG)*/
#define defMToC 1     /*Score for query M to ref C (AC)*/
#define defKToC -4    /*Score for query K to ref C (GT)*/
#define defRToC -4    /*Score for query R to ref C (AG)*/
#define defYToC 1     /*Score for query Y to ref C (CT)*/
#define defBToC -1    /*Score for query B to ref C (CGT)*/
#define defDToC -4    /*Score for query D to ref C (AGT)*/
#define defHToC -1    /*Score for query H to ref C (ACT)*/
#define defVToC -1    /*Score for query V to ref C (ACG)*/
#define defNToC -2    /*Score for query N to ref C (ACGT)*/
#define defXToC -2    /*Score for query X to ref C (ACGT)*/

/*Alignment scoring matrix anonymous (G)*/
#define defGToW -4    /*Score for query G to ref W (AT)*/
#define defGToS 1     /*Score for query G to ref S (CG)*/
#define defGToM -4    /*Score for query G to ref M (AC)*/
#define defGToK 1     /*Score for query G to ref K (GT)*/
#define defGToR 1     /*Score for query G to ref R (AG)*/
#define defGToY -4    /*Score for query G to ref Y (CT)*/
#define defGToB -1    /*Score for query G to ref B (CGT)*/
#define defGToD -1    /*Score for query G to ref D (AGT)*/
#define defGToH -4    /*Score for query G to ref H (ACT)*/
#define defGToV -1    /*Score for query G to ref V (ACG)*/
#define defGToN -2    /*Score for query G to ref N (ACGT)*/
#define defGToX -2    /*Score for query G to ref X (ACGT)*/

#define defWToG -4    /*Score for query G to ref W (AT)*/
#define defSToG 1     /*Score for query G to ref S (CG)*/
#define defMToG -4    /*Score for query G to ref M (AC)*/
#define defKToG 1     /*Score for query G to ref K (GT)*/
#define defRToG 1     /*Score for query G to ref R (AG)*/
#define defYToG -4    /*Score for query G to ref Y (CT)*/
#define defBToG -1    /*Score for query G to ref B (CGT)*/
#define defDToG -1    /*Score for query G to ref D (AGT)*/
#define defHToG -4    /*Score for query G to ref H (ACT)*/
#define defVToG -1    /*Score for query G to ref V (ACG)*/
#define defNToG -2    /*Score for query G to ref N (ACGT)*/
#define defXToG -2    /*Score for query G to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (T)*/
#define defTToW 1     /*Score for query T to ref W (AT)*/
#define defTToS -4    /*Score for query T to ref S (CG)*/
#define defTToM -4    /*Score for query T to ref M (AC)*/
#define defTToK 1     /*Score for query T to ref K (GT)*/
#define defTToR -4    /*Score for query T to ref R (AG)*/
#define defTToY 1     /*Score for query T to ref Y (CT)*/
#define defTToB -1    /*Score for query T to ref B (CGT)*/
#define defTToD -1    /*Score for query T to ref D (AGT)*/
#define defTToH -1    /*Score for query T to ref H (ACT)*/
#define defTToV -4    /*Score for query T to ref V (ACG)*/
#define defTToN -2    /*Score for query T to ref N (ACGT)*/
#define defTToX -2    /*Score for query T to ref X (ACGT)*/

#define defWToT 1     /*Score for query T to ref W (AT)*/
#define defSToT -4    /*Score for query T to ref S (CG)*/
#define defMToT -4    /*Score for query T to ref M (AC)*/
#define defKToT 1     /*Score for query T to ref K (GT)*/
#define defRToT -4    /*Score for query T to ref R (AG)*/
#define defYToT -4    /*Score for query T to ref Y (CT)*/
#define defBToT -1    /*Score for query T to ref B (CGT)*/
#define defDToT -1    /*Score for query T to ref D (AGT)*/
#define defHToT -1    /*Score for query T to ref H (ACT)*/
#define defVToT -4    /*Score for query T to ref V (ACG)*/
#define defNToT -2    /*Score for query T to ref N (ACGT)*/
#define defXToT -2    /*Score for query T to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (U)*/
#define defUToW 1     /*Score for query T to ref W (AT)*/
#define defUToS -4    /*Score for query T to ref S (CG)*/
#define defUToM -4    /*Score for query T to ref M (AC)*/
#define defUToK 1     /*Score for query T to ref K (GT)*/
#define defUToR -4    /*Score for query T to ref R (AG)*/
#define defUToY 1     /*Score for query T to ref Y (CT)*/
#define defUToB -1    /*Score for query T to ref B (CGT)*/
#define defUToD -1    /*Score for query T to ref D (AGT)*/
#define defUToH -1    /*Score for query T to ref H (ACT)*/
#define defUToV -4    /*Score for query T to ref V (ACG)*/
#define defUToN -2    /*Score for query T to ref N (ACGT)*/
#define defUToX -2    /*Score for query T to ref X (ACGT)*/

#define defWToU 1     /*Score for query T to ref W (AT)*/
#define defSToU -4    /*Score for query T to ref S (CG)*/
#define defMToU -4    /*Score for query T to ref M (AC)*/
#define defKToU 1     /*Score for query T to ref K (GT)*/
#define defRToU -4    /*Score for query T to ref R (AG)*/
#define defYToU -4    /*Score for query T to ref Y (CT)*/
#define defBToU -1    /*Score for query T to ref B (CGT)*/
#define defDToU -1    /*Score for query T to ref D (AGT)*/
#define defHToU -1    /*Score for query T to ref H (ACT)*/
#define defVToU -4    /*Score for query T to ref V (ACG)*/
#define defNToU -2    /*Score for query T to ref N (ACGT)*/
#define defXToU -2    /*Score for query T to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (W)*/
#define defWToW -1    /*Score for query W to ref W (AT)*/
#define defWToS -4    /*Score for query W to ref S (CG)*/
#define defWToM -1    /*Score for query W to ref M (AC)*/
#define defWToK -1    /*Score for query W to ref K (GT)*/
#define defWToR -1    /*Score for query W to ref R (AG)*/
#define defWToY -1    /*Score for query W to ref Y (CT)*/
#define defWToB -1    /*Score for query W to ref B (CGT)*/
#define defWToD -1    /*Score for query W to ref D (AGT)*/
#define defWToH -1    /*Score for query W to ref H (ACT)*/
#define defWToV -1    /*Score for query W to ref V (ACG)*/
#define defWToN -1    /*Score for query W to ref N (ACGT)*/
#define defWToX -1    /*Score for query W to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (S)*/
#define defSToW -4    /*Score for query S to ref W (AT)*/
#define defSToS -1    /*Score for query S to ref S (CG)*/
#define defSToM -2    /*Score for query S to ref M (AC)*/
#define defSToK -2    /*Score for query S to ref K (GT)*/
#define defSToR -2    /*Score for query S to ref R (AG)*/
#define defSToY -2    /*Score for query S to ref Y (CT)*/
#define defSToB -1    /*Score for query S to ref B (CGT)*/
#define defSToD -3    /*Score for query S to ref D (AGT)*/
#define defSToH -3    /*Score for query S to ref H (ACT)*/
#define defSToV -1    /*Score for query S to ref V (ACG)*/
#define defSToN -1    /*Score for query S to ref N (ACGT)*/
#define defSToX -1    /*Score for query S to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (M)*/
#define defMToW -2    /*Score for query M to ref W (AT)*/
#define defMToS -2    /*Score for query M to ref S (CG)*/
#define defMToM -1    /*Score for query M to ref M (AC)*/
#define defMToK -4    /*Score for query M to ref K (GT)*/
#define defMToR -2    /*Score for query M to ref R (AG)*/
#define defMToY -2    /*Score for query M to ref Y (CT)*/
#define defMToB -3    /*Score for query M to ref B (CGT)*/
#define defMToD -3    /*Score for query M to ref D (AGT)*/
#define defMToH -1    /*Score for query M to ref H (ACT)*/
#define defMToV -1    /*Score for query M to ref V (ACG)*/
#define defMToN -1    /*Score for query M to ref N (ACGT)*/
#define defMToX -1    /*Score for query M to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (K)*/
#define defKToW -2    /*Score for query K to ref W (AT)*/
#define defKToS -2    /*Score for query K to ref S (CG)*/
#define defKToM -4    /*Score for query K to ref M (AC)*/
#define defKToK -1    /*Score for query K to ref K (GT)*/
#define defKToR -2    /*Score for query K to ref R (AG)*/
#define defKToY -2    /*Score for query K to ref Y (CT)*/
#define defKToB -1    /*Score for query K to ref B (CGT)*/
#define defKToD -1    /*Score for query K to ref D (AGT)*/
#define defKToH -3    /*Score for query K to ref H (ACT)*/
#define defKToV -3    /*Score for query K to ref V (ACG)*/
#define defKToN -1    /*Score for query K to ref N (ACGT)*/
#define defKToX -1    /*Score for query K to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (R)*/
#define defRToW -1    /*Score for query R to ref W (AT)*/
#define defRToS -1    /*Score for query R to ref S (CG)*/
#define defRToM -1    /*Score for query R to ref M (AC)*/
#define defRToK -2    /*Score for query R to ref K (GT)*/
#define defRToR -1    /*Score for query R to ref R (AG)*/
#define defRToY -4    /*Score for query R to ref Y (CT)*/
#define defRToB -3    /*Score for query R to ref B (CGT)*/
#define defRToD -1    /*Score for query R to ref D (AGT)*/
#define defRToH -1    /*Score for query R to ref H (ACT)*/
#define defRToV -3    /*Score for query R to ref V (ACG)*/
#define defRToN -1    /*Score for query R to ref N (ACGT)*/
#define defRToX -1    /*Score for query R to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (Y)*/
#define defYToW -2    /*Penalty for query Y to ref W (AT)*/
#define defYToS -2    /*Score for query Y to ref S (CG)*/
#define defYToM -2    /*Score for query Y to ref M (AC)*/
#define defYToK -2    /*Score for query Y to ref K (GT)*/
#define defYToR -4    /*Score for query Y to ref R (AG)*/
#define defYToY -1    /*Score for query Y to ref Y (CT)*/
#define defYToB -2    /*Score for query Y to ref B (CGT)*/
#define defYToD -3    /*Score for query Y to ref D (AGT)*/
#define defYToH -1    /*Score for query Y to ref H (ACT)*/
#define defYToV -3    /*Score for query Y to ref V (ACG)*/
#define defYToN -1    /*Score for query Y to ref N (ACGT)*/
#define defYToX -1    /*Score for query Y to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (B)*/
#define defBToW -3    /*Score for query B to ref W (AT)*/
#define defBToS -1    /*Score for query B to ref S (CG)*/
#define defBToM -3    /*Score for query B to ref M (AC)*/
#define defBToK -1    /*Score for query B to ref K (GT)*/
#define defBToR -3    /*Score for query B to ref R (AG)*/
#define defBToY -1    /*Score for query B to ref Y (CT)*/
#define defBToB -1    /*Score for query B to ref B (CGT)*/
#define defBToD -2    /*Score for query B to ref D (AGT)*/
#define defBToH -2    /*Score for query B to ref H (ACT)*/
#define defBToV -2    /*Score for query B to ref V (ACG)*/
#define defBToN -1    /*Score for query B to ref N (ACGT)*/
#define defBToX -1    /*Score for query B to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (D)*/
#define defDToW -1    /*Score for query D to ref W (AT)*/
#define defDToS -3    /*Score for query D to ref S (CG)*/
#define defDToM -3    /*Score for query D to ref M (AC)*/
#define defDToK -1    /*Score for query D to ref K (GT)*/
#define defDToR -1    /*Score for query D to ref R (AG)*/
#define defDToY -3    /*Score for query D to ref Y (CT)*/
#define defDToB -2    /*Score for query D to ref B (CGT)*/
#define defDToD -1    /*Score for query D to ref D (AGT)*/
#define defDToH -2    /*Score for query D to ref H (ACT)*/
#define defDToV -2    /*Score for query D to ref V (ACG)*/
#define defDToN -1    /*Score for query D to ref N (ACGT)*/
#define defDToX -1    /*Score for query D to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (H)*/
#define defHToW -1    /*Score for query H to ref W (AT)*/
#define defHToS -3    /*Score for query H to ref S (CG)*/
#define defHToM -1    /*Score for query H to ref M (AC)*/
#define defHToK -3    /*Score for query H to ref K (GT)*/
#define defHToR -3    /*Score for query H to ref R (AG)*/
#define defHToY -1    /*Score for query H to ref Y (CT)*/
#define defHToB -2    /*Score for query H to ref B (CGT)*/
#define defHToD -2    /*Score for query H to ref D (AGT)*/
#define defHToH -1    /*Score for query H to ref H (ACT)*/
#define defHToV -2    /*Score for query H to ref V (ACG)*/
#define defHToN -1    /*Score for query H to ref N (ACGT)*/
#define defHToX -1    /*Score for query H to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (V)*/
#define defVToW -3    /*Score for query V to ref W (AT)*/
#define defVToS -1    /*Score for query V to ref S (CG)*/
#define defVToM -1    /*Score for query V to ref M (AC)*/
#define defVToK -3    /*Score for query V to ref K (GT)*/
#define defVToR -1    /*Score for query V to ref R (AG)*/
#define defVToY -3    /*Score for query V to ref Y (CT)*/
#define defVToB -2    /*Score for query V to ref B (CGT)*/
#define defVToD -2    /*Score for query V to ref D (AGT)*/
#define defVToH -2    /*Score for query V to ref H (ACT)*/
#define defVToV -1    /*Score for query V to ref V (ACG)*/
#define defVToN -1    /*Score for query V to ref N (ACGT)*/
#define defVToX -1    /*Score for query V to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (N)*/
#define defNToW -2    /*Score for query N to ref W (AT)*/
#define defNToS -2    /*Score for query N to ref S (CG)*/
#define defNToM -2    /*Score for query N to ref M (AC)*/
#define defNToK -2    /*Score for query N to ref K (GT)*/
#define defNToR -2    /*Score for query N to ref R (AG)*/
#define defNToY -2    /*Score for query N to ref Y (CT)*/
#define defNToB -2    /*Score for query N to ref B (CGT)*/
#define defNToD -2    /*Score for query N to ref D (AGT)*/
#define defNToH -2    /*Score for query N to ref H (ACT)*/
#define defNToV -2    /*Score for query N to ref V (ACG)*/
#define defNToN -1    /*Score for query N to ref N (ACGT)*/
#define defNToX -1    /*Score for query N to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (X)*/
#define defXToW -1    /*Score for query X to ref W (AT)*/
#define defXToS -1    /*Score for query X to ref S (CG)*/
#define defXToM -1    /*Score for query X to ref M (AC)*/
#define defXToK -1    /*Score for query X to ref K (GT)*/
#define defXToR -1    /*Score for query X to ref R (AG)*/
#define defXToY -1    /*Score for query X to ref Y (CT)*/
#define defXToB -1    /*Score for query X to ref B (CGT)*/
#define defXToD -1    /*Score for query X to ref D (AGT)*/
#define defXToH -1    /*Score for query X to ref H (ACT)*/
#define defXToV -1    /*Score for query X to ref V (ACG)*/
#define defXToN -1    /*Score for query X to ref N (ACGT)*/
#define defXToX -1    /*Score for query X to ref X (ACGT)*/

/*Anonymous bases for matching*/

#define defUMatchT 1     /*query U to ref T*/
#define defTMatchU 1     /*query T to ref U*/

/*A as query*/
#define defAMatchW 1     /*query A to ref W (AT)*/
#define defAMatchS 0    /*query A to ref S (CG)*/
#define defAMatchM 1     /*query A to ref M (AC)*/
#define defAMatchK 0    /*query A to ref K (GT)*/
#define defAMatchR 1     /*query A to ref R (AG)*/
#define defAMatchY 0    /*query A to ref Y (CT)*/
#define defAMatchB 0    /*query A to ref B (CGT)*/
#define defAMatchD 1    /*query A to ref D (AGT)*/
#define defAMatchH 1    /*query A to ref H (ACT)*/
#define defAMatchV 1    /*query A to ref V (ACG)*/
#define defAMatchN 1    /*query A to ref N (ACGT)*/
#define defAMatchX 1    /*query A to ref X (ACGT)*/

/*A as reference*/
#define defWMatchA 1     /*query W to ref A (AT)*/
#define defSMatchA 0    /*query S to ref A (CG)*/
#define defMMatchA 1     /*query M to ref A (AC)*/
#define defKMatchA 0    /*query K to ref A (GT)*/
#define defRMatchA 1     /*query R to ref A (AG)*/
#define defYMatchA 0    /*query Y to ref A (CT)*/
#define defBMatchA 0    /*query B to ref A (CGT)*/
#define defDMatchA 1    /*query D to ref A (AGT)*/
#define defHMatchA 1    /*query H to ref A (ACT)*/
#define defVMatchA 1    /*query V to ref A (ACG)*/
#define defNMatchA 1    /*query N to ref A (ACGT)*/
#define defXMatchA 1    /*query X to ref A (ACGT)*/

/*C as query*/
#define defCMatchW 0    /*query C to ref W (AT)*/
#define defCMatchS 1     /*query C to ref S (CG)*/
#define defCMatchM 1     /*query C to ref M (AC)*/
#define defCMatchK 0    /*query C to ref K (GT)*/
#define defCMatchR 0    /*query C to ref R (AG)*/
#define defCMatchY 1     /*query C to ref Y (CT)*/
#define defCMatchB 1    /*query C to ref B (CGT)*/
#define defCMatchD 0    /*query C to ref D (AGT)*/
#define defCMatchH 1    /*query C to ref H (ACT)*/
#define defCMatchV 1    /*query C to ref V (ACG)*/
#define defCMatchN 1    /*query C to ref N (ACGT)*/
#define defCMatchX 1    /*query C to ref X (ACGT)*/

/*C as reference*/
#define defWMatchC 0    /*query W to ref C (AT)*/
#define defSMatchC 1     /*query S to ref C (CG)*/
#define defMMatchC 1     /*query M to ref C (AC)*/
#define defKMatchC 0    /*query K to ref C (GT)*/
#define defRMatchC 0    /*query R to ref C (AG)*/
#define defYMatchC 1     /*query Y to ref C (CT)*/
#define defBMatchC 1    /*query B to ref C (CGT)*/
#define defDMatchC 0    /*query D to ref C (AGT)*/
#define defHMatchC 1    /*query H to ref C (ACT)*/
#define defVMatchC 1    /*query V to ref C (ACG)*/
#define defNMatchC 1    /*query N to ref C (ACGT)*/
#define defXMatchC 1    /*query X to ref C (ACGT)*/

/*G as query*/
#define defGMatchW 0    /*query G to ref W (AT)*/
#define defGMatchS 1     /*query G to ref S (CG)*/
#define defGMatchM 0    /*query G to ref M (AC)*/
#define defGMatchK 1     /*query G to ref K (GT)*/
#define defGMatchR 1     /*query G to ref R (AG)*/
#define defGMatchY 0    /*query G to ref Y (CT)*/
#define defGMatchB 1    /*query G to ref B (CGT)*/
#define defGMatchD 1    /*query G to ref D (AGT)*/
#define defGMatchH 0    /*query G to ref H (ACT)*/
#define defGMatchV 1    /*query G to ref V (ACG)*/
#define defGMatchN 1    /*query G to ref N (ACGT)*/
#define defGMatchX 1    /*query G to ref X (ACGT)*/

/*G as reference*/
#define defWMatchG 0    /*query G to ref W (AT)*/
#define defSMatchG 1     /*query G to ref S (CG)*/
#define defMMatchG 0    /*query G to ref M (AC)*/
#define defKMatchG 1     /*query G to ref K (GT)*/
#define defRMatchG 1     /*query G to ref R (AG)*/
#define defYMatchG 0    /*query G to ref Y (CT)*/
#define defBMatchG 1    /*query G to ref B (CGT)*/
#define defDMatchG 1    /*query G to ref D (AGT)*/
#define defHMatchG 0    /*query G to ref H (ACT)*/
#define defVMatchG 1    /*query G to ref V (ACG)*/
#define defNMatchG 1    /*query G to ref N (ACGT)*/
#define defXMatchG 1    /*query G to ref X (ACGT)*/

/*T as query*/
#define defTMatchW 1     /*query T to ref W (AT)*/
#define defTMatchS 0    /*query T to ref S (CG)*/
#define defTMatchM 0    /*query T to ref M (AC)*/
#define defTMatchK 1     /*query T to ref K (GT)*/
#define defTMatchR 0    /*query T to ref R (AG)*/
#define defTMatchY 1     /*query T to ref Y (CT)*/
#define defTMatchB 1    /*query T to ref B (CGT)*/
#define defTMatchD 1    /*query T to ref D (AGT)*/
#define defTMatchH 1    /*query T to ref H (ACT)*/
#define defTMatchV 0    /*query T to ref V (ACG)*/
#define defTMatchN 1    /*query T to ref N (ACGT)*/
#define defTMatchX 1    /*query T to ref X (ACGT)*/

/*T as reference*/
#define defWMatchT 1    /*query T to ref W (AT)*/
#define defSMatchT 0    /*query T to ref S (CG)*/
#define defMMatchT 0    /*query T to ref M (AC)*/
#define defKMatchT 1    /*query T to ref K (GT)*/
#define defRMatchT 0    /*query T to ref R (AG)*/
#define defYMatchT 1    /*query T to ref Y (CT)*/
#define defBMatchT 1    /*query T to ref B (CGT)*/
#define defDMatchT 1    /*query T to ref D (AGT)*/
#define defHMatchT 1    /*query T to ref H (ACT)*/
#define defVMatchT 0    /*query T to ref V (ACG)*/
#define defNMatchT 1    /*query T to ref N (ACGT)*/
#define defXMatchT 1    /*query T to ref X (ACGT)*/

/*U as query*/
#define defUMatchW 1     /*query T to ref W (AT)*/
#define defUMatchS 0    /*query T to ref S (CG)*/
#define defUMatchM 0    /*query T to ref M (AC)*/
#define defUMatchK 1     /*query T to ref K (GT)*/
#define defUMatchR 0    /*query T to ref R (AG)*/
#define defUMatchY 1     /*query T to ref Y (CT)*/
#define defUMatchB 1    /*query T to ref B (CGT)*/
#define defUMatchD 1    /*query T to ref D (AGT)*/
#define defUMatchH 1    /*query T to ref H (ACT)*/
#define defUMatchV 0    /*query T to ref V (ACG)*/
#define defUMatchN 1    /*query T to ref N (ACGT)*/
#define defUMatchX 1    /*query T to ref X (ACGT)*/

/*U as reference*/
#define defWMatchU 1     /*query U to ref T*/
#define defSMatchU 0    /*query T to ref S (CG)*/
#define defMMatchU 0    /*query T to ref M (AC)*/
#define defKMatchU 1     /*query T to ref K (GT)*/
#define defRMatchU 0    /*query T to ref R (AG)*/
#define defYMatchU 1    /*query T to ref Y (CT)*/
#define defBMatchU 1    /*query T to ref B (CGT)*/
#define defDMatchU 1    /*query T to ref D (AGT)*/
#define defHMatchU 1    /*query T to ref H (ACT)*/
#define defVMatchU 0    /*query T to ref V (ACG)*/
#define defNMatchU 1    /*query T to ref N (ACGT)*/
#define defXMatchU 1    /*query T to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (W)*/
#define defWMatchW 1    /*query W to ref W (AT)*/
#define defWMatchS 0    /*query W to ref S (CG)*/
#define defWMatchM 1    /*query W to ref M (AC)*/
#define defWMatchK 1    /*query W to ref K (GT)*/
#define defWMatchR 1    /*query W to ref R (AG)*/
#define defWMatchY 1    /*query W to ref Y (CT)*/
#define defWMatchB 1    /*query W to ref B (CGT)*/
#define defWMatchD 1    /*query W to ref D (AGT)*/
#define defWMatchH 1    /*query W to ref H (ACT)*/
#define defWMatchV 1    /*query W to ref V (ACG)*/
#define defWMatchN 1    /*query W to ref N (ACGT)*/
#define defWMatchX 1    /*query W to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (S)*/
#define defSMatchW 0    /*query S to ref W (AT)*/
#define defSMatchS 1    /*query S to ref S (CG)*/
#define defSMatchM 1    /*query S to ref M (AC)*/
#define defSMatchK 1    /*query S to ref K (GT)*/
#define defSMatchR 1    /*query S to ref R (AG)*/
#define defSMatchY 1    /*query S to ref Y (CT)*/
#define defSMatchB 1    /*query S to ref B (CGT)*/
#define defSMatchD 1    /*query S to ref D (AGT)*/
#define defSMatchH 1    /*query S to ref H (ACT)*/
#define defSMatchV 1    /*query S to ref V (ACG)*/
#define defSMatchN 1    /*query S to ref N (ACGT)*/
#define defSMatchX 1    /*query S to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (M)*/
#define defMMatchW 1    /*query M to ref W (AT)*/
#define defMMatchS 1    /*query M to ref S (CG)*/
#define defMMatchM 1    /*query M to ref M (AC)*/
#define defMMatchK 0    /*query M to ref K (GT)*/
#define defMMatchR 1    /*query M to ref R (AG)*/
#define defMMatchY 1    /*query M to ref Y (CT)*/
#define defMMatchB 1    /*query M to ref B (CGT)*/
#define defMMatchD 1    /*query M to ref D (AGT)*/
#define defMMatchH 1    /*query M to ref H (ACT)*/
#define defMMatchV 1    /*query M to ref V (ACG)*/
#define defMMatchN 1    /*query M to ref N (ACGT)*/
#define defMMatchX 1    /*query M to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (K)*/
#define defKMatchW 1    /*query K to ref W (AT)*/
#define defKMatchS 1    /*query K to ref S (CG)*/
#define defKMatchM 0    /*query K to ref M (AC)*/
#define defKMatchK 1    /*query K to ref K (GT)*/
#define defKMatchR 1    /*query K to ref R (AG)*/
#define defKMatchY 1    /*query K to ref Y (CT)*/
#define defKMatchB 1    /*query K to ref B (CGT)*/
#define defKMatchD 1    /*query K to ref D (AGT)*/
#define defKMatchH 1    /*query K to ref H (ACT)*/
#define defKMatchV 1    /*query K to ref V (ACG)*/
#define defKMatchN 1    /*query K to ref N (ACGT)*/
#define defKMatchX 1    /*query K to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (R)*/
#define defRMatchW 1    /*query R to ref W (AT)*/
#define defRMatchS 1    /*query R to ref S (CG)*/
#define defRMatchM 1    /*query R to ref M (AC)*/
#define defRMatchK 1    /*query R to ref K (GT)*/
#define defRMatchR 1    /*query R to ref R (AG)*/
#define defRMatchY 0    /*query R to ref Y (CT)*/
#define defRMatchB 1    /*query R to ref B (CGT)*/
#define defRMatchD 1    /*query R to ref D (AGT)*/
#define defRMatchH 1    /*query R to ref H (ACT)*/
#define defRMatchV 1    /*query R to ref V (ACG)*/
#define defRMatchN 1    /*query R to ref N (ACGT)*/
#define defRMatchX 1    /*query R to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (Y)*/
#define defYMatchW 1    /*r query Y to ref W (AT)*/
#define defYMatchS 1    /*query Y to ref S (CG)*/
#define defYMatchM 1    /*query Y to ref M (AC)*/
#define defYMatchK 1    /*query Y to ref K (GT)*/
#define defYMatchR 0    /*query Y to ref R (AG)*/
#define defYMatchY 1    /*query Y to ref Y (CT)*/
#define defYMatchB 1    /*query Y to ref B (CGT)*/
#define defYMatchD 1    /*query Y to ref D (AGT)*/
#define defYMatchH 1    /*query Y to ref H (ACT)*/
#define defYMatchV 1    /*query Y to ref V (ACG)*/
#define defYMatchN 1    /*query Y to ref N (ACGT)*/
#define defYMatchX 1    /*query Y to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (B)*/
#define defBMatchW 1    /*query B to ref W (AT)*/
#define defBMatchS 1    /*query B to ref S (CG)*/
#define defBMatchM 1    /*query B to ref M (AC)*/
#define defBMatchK 1    /*query B to ref K (GT)*/
#define defBMatchR 1    /*query B to ref R (AG)*/
#define defBMatchY 1    /*query B to ref Y (CT)*/
#define defBMatchB 1    /*query B to ref B (CGT)*/
#define defBMatchD 1    /*query B to ref D (AGT)*/
#define defBMatchH 1    /*query B to ref H (ACT)*/
#define defBMatchV 1    /*query B to ref V (ACG)*/
#define defBMatchN 1    /*query B to ref N (ACGT)*/
#define defBMatchX 1    /*query B to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (D)*/
#define defDMatchW 1    /*query D to ref W (AT)*/
#define defDMatchS 1    /*query D to ref S (CG)*/
#define defDMatchM 1    /*query D to ref M (AC)*/
#define defDMatchK 1    /*query D to ref K (GT)*/
#define defDMatchR 1    /*query D to ref R (AG)*/
#define defDMatchY 1    /*query D to ref Y (CT)*/
#define defDMatchB 1    /*query D to ref B (CGT)*/
#define defDMatchD 1    /*query D to ref D (AGT)*/
#define defDMatchH 1    /*query D to ref H (ACT)*/
#define defDMatchV 1    /*query D to ref V (ACG)*/
#define defDMatchN 1    /*query D to ref N (ACGT)*/
#define defDMatchX 1    /*query D to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (H)*/
#define defHMatchW 1    /*query H to ref W (AT)*/
#define defHMatchS 1    /*query H to ref S (CG)*/
#define defHMatchM 1    /*query H to ref M (AC)*/
#define defHMatchK 1    /*query H to ref K (GT)*/
#define defHMatchR 1    /*query H to ref R (AG)*/
#define defHMatchY 1    /*query H to ref Y (CT)*/
#define defHMatchB 1    /*query H to ref B (CGT)*/
#define defHMatchD 1    /*query H to ref D (AGT)*/
#define defHMatchH 1    /*query H to ref H (ACT)*/
#define defHMatchV 1    /*query H to ref V (ACG)*/
#define defHMatchN 1    /*query H to ref N (ACGT)*/
#define defHMatchX 1    /*query H to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (V)*/
#define defVMatchW 1    /*query V to ref W (AT)*/
#define defVMatchS 1    /*query V to ref S (CG)*/
#define defVMatchM 1    /*query V to ref M (AC)*/
#define defVMatchK 1    /*query V to ref K (GT)*/
#define defVMatchR 1    /*query V to ref R (AG)*/
#define defVMatchY 1    /*query V to ref Y (CT)*/
#define defVMatchB 1    /*query V to ref B (CGT)*/
#define defVMatchD 1    /*query V to ref D (AGT)*/
#define defVMatchH 1    /*query V to ref H (ACT)*/
#define defVMatchV 1    /*query V to ref V (ACG)*/
#define defVMatchN 1    /*query V to ref N (ACGT)*/
#define defVMatchX 1    /*query V to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (N)*/
#define defNMatchW 1    /*query N to ref W (AT)*/
#define defNMatchS 1    /*query N to ref S (CG)*/
#define defNMatchM 1    /*query N to ref M (AC)*/
#define defNMatchK 1    /*query N to ref K (GT)*/
#define defNMatchR 1    /*query N to ref R (AG)*/
#define defNMatchY 1    /*query N to ref Y (CT)*/
#define defNMatchB 1    /*query N to ref B (CGT)*/
#define defNMatchD 1    /*query N to ref D (AGT)*/
#define defNMatchH 1    /*query N to ref H (ACT)*/
#define defNMatchV 1    /*query N to ref V (ACG)*/
#define defNMatchN 1    /*query N to ref N (ACGT)*/
#define defNMatchX 1    /*query N to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (X)*/
#define defXMatchW 1    /*query X to ref W (AT)*/
#define defXMatchS 1    /*query X to ref S (CG)*/
#define defXMatchM 1    /*query X to ref M (AC)*/
#define defXMatchK 1    /*query X to ref K (GT)*/
#define defXMatchR 1    /*query X to ref R (AG)*/
#define defXMatchY 1    /*query X to ref Y (CT)*/
#define defXMatchB 1    /*query X to ref B (CGT)*/
#define defXMatchD 1    /*query X to ref D (AGT)*/
#define defXMatchH 1    /*query X to ref H (ACT)*/
#define defXMatchV 1    /*query X to ref V (ACG)*/
#define defXMatchN 1    /*query X to ref N (ACGT)*/
#define defXMatchX 1    /*query X to ref X (ACGT)*/


#endif
