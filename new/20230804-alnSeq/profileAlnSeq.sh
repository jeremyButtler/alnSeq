algoStr="$1";
algoParmStr="";

smallTestQueryStr="small-test-query.fasta";
smallTestRefStr="small-test-ref.fasta";

midTestQueryStr="TickBornEncephalitis-query.fasta";
midTestRefStr="TickBornEncephalitis-reference.fasta";

largeTestQueryStr="PEDV_Austria_2015.fasta";
largeTestRefStr="PEDV_Germany_2015.fasta";

#oldProgramStr="20230804-alnSeq/alnSeq";
oldProgramStr="old-alnSeq/alnSeq";
   # Currently old changed version
newProgramStr="new-alnSeq/alnSeq";

numRepI=5;
statsFileStr="alnSeq-bench-changes.tsv"
iCnt=0;
genomeSizeStr="";

if [[ "$algoStr" == "needle" ]]; then
# If using the needleman alignment
    algoStr="needle";
    algoParmStr="-use-needle";
# If using the needleman alignment

elif [[ "$algoStr" == "water" ]]; then
# else If using the waterman alignment
    algoStr="water";
    algoParmStr="-use-water";
# else If using the waterman alignment

elif [[ "$algoStr" == "hirsch" ]]; then
# else If using the waterman alignment
    algoStr="hirsch";
    algoParmStr="-use-hirschberg";
# else If using the waterman alignment

else 
# else using the default (needle) alignment
    algoStr="needle";
    algoParmStr="-use-needle";
# else using the default (needle) alignment
fi

# Make sure using the correct programs
cd "$(dirname "$newProgramStr")" && make && cd ..;
cd "$(dirname "$oldProgramStr")" && make && cd ..;

{
    printf "Program\ttest\talgorithm\telapsedTime\t";
    printf "userTime\tsystemTime\tmaxResidentMemory\tCPU";
    printf "\n";
} > "$statsFileStr";

for strRef in "$smallTestRefStr" "$midTestRefStr" "$largeTestRefStr"; do
# For all references
    iCnt=0;

    while [[ "$iCnt" -lt "$numRepI" ]]; do
    # While I have replicates to run
    
       iCnt="$((iCnt + 1))";

       if [[ "$strRef" == "$smallTestRefStr" ]]; then
         strQuery="$smallTestQueryStr";
         genomeSizeStr="small";
       elif [[ "$strRef" == "$midTestRefStr" ]]; then
         strQuery="$midTestQueryStr";
         genomeSizeStr="mid";
       elif [[ "$strRef" == "$largeTestRefStr" ]]; then
         strQuery="$largeTestQueryStr";
         genomeSizeStr="large";
       else
         continue; # null entry??
       fi # done

       printf "Genome: %5s\tRep: %s\n" "$genomeSizeStr" "$iCnt";

       /usr/bin/time \
           -f "old\t$genomeSizeStr\t$algoStr\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$oldProgramStr" \
           "$algoParmStr" \
           -query "$strQuery" \
           -ref "$strRef" \
           -out "tmp.aln";

       /usr/bin/time \
           -f "new\t$genomeSizeStr\t$algoStr\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$newProgramStr" \
           "$algoParmStr" \
           -query "$strQuery" \
           -ref "$strRef" \
           -out "tmp.aln";
    done # While I have replicates to run
done # For all references
