for i in {25..40};
do
echo "Processing sample ${i}"
# Assuming salmon is one directory above the reads directory, ie
# in the directory where salmon resides there is a directory reads.
salmon quant -i athal_index -l A \
         -r ./reads/${i}_1.fastq.gz \
            ./reads/${i}_2.fastq.gz \
         -p 8 --gcBias -o quants/${i}_quant
done 
