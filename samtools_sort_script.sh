for file in *.bam;
do
Base_Name=${file%.bam}
echo "samtools sort -@ 8 ${Base_Name}.bam > bam_sorted/${Base_Name}_sorted.bam"
done