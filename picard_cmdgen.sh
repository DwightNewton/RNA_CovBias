for file in *.bam;
do
Base_Name=${file%.bam}
echo "java -jar /external/mgmt3/quarantine/Picard/2.90/picard.jar CollectRnaSeqMetrics \
REFERENCE_SEQUENCE=\"/external/rprshnas01/eslab/dnewton/Ram_MDD/Reference/fasta/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa\" \
REF_FLAT=\"Grch38-r98-refFlat-rearranged.txt\" \
INPUT=\"/external/rprshnas01/eslab/dnewton/PITT_tetrad_bam_5_3__3_15/all_bam_files/bam_sorted/${Base_Name}.bam\" \
OUTPUT=\"/external/rprshnas01/eslab/dnewton/PITT_tetrad_bam_5_3__3_15/all_bam_files/picard_output/${Base_Name}_picardoutput.txt\" \
STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND \
CHART_OUTPUT=\"/external/rprshnas01/eslab/dnewton/PITT_tetrad_bam_5_3__3_15/all_bam_files/picard_output/${Base_Name}_picardoutput.pdf\""
done