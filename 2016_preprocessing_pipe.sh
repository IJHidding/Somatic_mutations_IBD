#!/bin/bash

#SBATCH --job-name=2016_Preprocessing
#SBATCH --output=jobName.out
#SBATCH --error=jobName.err
#SBATCH --time=145:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

uniq_run_id=$(echo $RANDOM)
biopt_folder=$1
wes_1=$2
wes_file=$( echo $2 | sed 's/C1933_//g')
fastq_name_path=$(echo $1) # | cut -d"." -f1)
wes_name_path=$(echo $wes_file | cut -d"." -f1)
fastq_name=$(echo ${fastq_name_path##*/})
wes_name=$(echo ${wes_name_path##*/})

echo $wes_name
echo $fastq_name

if [ -f "$wes_file" ]; then
	echo "$wes_file exists."
else 
	echo "$wes_file does not exist, exiting..." >> /groups/umcg-weersma/tmp01/Iwan/pipeline/log.txt
	exit 1
fi


#Star aligner twice 

module load STAR
COUNTER=1
for file in $(grep -v "md5" $1 | sort | uniq)  
do
	echo $file
	if [ -f "$file" ]; then
        	echo "$file exists."
	else
    		echo "$file from folder ${biopt_folder} does not exist, exiting..." >> /groups/umcg-weersma/tmp01/Iwan/pipeline/log.txt
        	exit 1
	fi

	STAR --genomeDir /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38 --outSAMmapqUnique 60 --readFilesIn <(gunzip -c ${file}) --outFileNamePrefix ${fastq_name_path}

	mv ${fastq_name_path}Aligned.out.sam ${uniq_run_id}${fastq_name}fastq_${COUNTER}.sam
	(( COUNTER++ ))
done

#sleep 1m
# merge

module load SAMtools

samtools view -S -b ${uniq_run_id}${fastq_name}fastq_1.sam > ${uniq_run_id}${fastq_name}1.bam
samtools view -S -b ${uniq_run_id}${fastq_name}fastq_2.sam > ${uniq_run_id}${fastq_name}2.bam
samtools view -S -b ${uniq_run_id}${fastq_name}fastq_3.sam > ${uniq_run_id}${fastq_name}3.bam
samtools view -S -b ${uniq_run_id}${fastq_name}fastq_4.sam > ${uniq_run_id}${fastq_name}4.bam
samtools view -S -b ${uniq_run_id}${fastq_name}fastq_5.sam > ${uniq_run_id}${fastq_name}5.bam
samtools view -S -b ${uniq_run_id}${fastq_name}fastq_6.sam > ${uniq_run_id}${fastq_name}6.bam
samtools view -S -b ${uniq_run_id}${fastq_name}fastq_7.sam > ${uniq_run_id}${fastq_name}7.bam
samtools view -S -b ${uniq_run_id}${fastq_name}fastq_8.sam > ${uniq_run_id}${fastq_name}8.bam

samtools merge ${uniq_run_id}${fastq_name}.bam ${uniq_run_id}${fastq_name}1.bam ${uniq_run_id}${fastq_name}2.bam ${uniq_run_id}${fastq_name}3.bam ${uniq_run_id}${fastq_name}4.bam ${uniq_run_id}${fastq_name}5.bam ${uniq_run_id}${fastq_name}6.bam ${uniq_run_id}${fastq_name}7.bam ${uniq_run_id}${fastq_name}8.bam

rm ${uniq_run_id}${fastq_name}1.bam
rm ${uniq_run_id}${fastq_name}2.bam
rm ${uniq_run_id}${fastq_name}3.bam
rm ${uniq_run_id}${fastq_name}4.bam
rm ${uniq_run_id}${fastq_name}5.bam
rm ${uniq_run_id}${fastq_name}6.bam
rm ${uniq_run_id}${fastq_name}7.bam
rm ${uniq_run_id}${fastq_name}8.bam
rm ${uniq_run_id}${fastq_name}fastq_1.sam
rm ${uniq_run_id}${fastq_name}fastq_2.sam
rm ${uniq_run_id}${fastq_name}fastq_3.sam
rm ${uniq_run_id}${fastq_name}fastq_4.sam
rm ${uniq_run_id}${fastq_name}fastq_5.sam
rm ${uniq_run_id}${fastq_name}fastq_6.sam
rm ${uniq_run_id}${fastq_name}fastq_7.sam
rm ${uniq_run_id}${fastq_name}fastq_8.sam


#rm ${fastq_name_path}Aligned.out.sam

# convert cram to bam
samtools view -b -T /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/Homo_sapiens_assembly38.fasta -o ${uniq_run_id}${wes_name}.bam ${wes_file}


# readgroups
module load picard

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I=${uniq_run_id}${fastq_name}.bam \
O=${uniq_run_id}output_${fastq_name}.bam \
RGID=4 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=${fastq_name}

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I=${uniq_run_id}${wes_name}.bam \
O=${uniq_run_id}output_${wes_name}.bam \
RGID=4 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=${wes_name}

rm ${uniq_run_id}${wes_name}.bam
rm ${uniq_run_id}${fastq_name}.bam

# sort and index

# add filter steps here  
module load GATK

gatk BaseRecalibrator \
-I ${uniq_run_id}output_${fastq_name}.bam \
-R /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/Homo_sapiens_assembly38.fasta \
--known-sites /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/somatic-hg38_af-only-gnomad.hg38.vcf \
--known-sites /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/somatic-hg38_1000g_pon.hg38.vcf \
-O ${uniq_run_id}${fastq_name}_recal_data.table

gatk ApplyBQSR \
-R /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/Homo_sapiens_assembly38.fasta \
-I ${uniq_run_id}output_${fastq_name}.bam \
--bqsr-recal-file ${uniq_run_id}${fastq_name}_recal_data.table \
-O ${uniq_run_id}${fastq_name}_output_test.bam

rm ${uniq_run_id}${fastq_name}_recal_data.table
mv ${uniq_run_id}${fastq_name}_output_test.bam ${uniq_run_id}output_${fastq_name}.bam


gatk MarkDuplicatesSpark \
-I ${uniq_run_id}output_${fastq_name}.bam \
-O ${uniq_run_id}${fastq_name}_marked_duplicates.bam \
--remove-sequencing-duplicates \
--tmp-dir /groups/umcg-weersma/tmp01/Iwan/pipeline/tmp

mv ${uniq_run_id}${fastq_name}_marked_duplicated.bam ${uniq_run_id}output_${fastq_name}.bam

gatk BaseRecalibrator \
-I ${uniq_run_id}output_${wes_name}.bam \
-R /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/Homo_sapiens_assembly38.fasta \
--known-sites /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/somatic-hg38_af-only-gnomad.hg38.vcf \
--known-sites /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/somatic-hg38_1000g_pon.hg38.vcf \
-O ${uniq_run_id}${wes_name}_recal_data.table

gatk ApplyBQSR \
-R /groups/umcg-weersma/tmp01/Iwan/datastorage/reference/hg38/Homo_sapiens_assembly38.fasta \
-I ${uniq_run_id}output_${wes_name}.bam \
--bqsr-recal-file ${uniq_run_id}${wes_name}_recal_data.table \
-O ${uniq_run_id}${wes_name}_output_test.bam

mv ${uniq_run_id}${wes_name}_output_test.bam ${uniq_run_id}output_${wes_name}.bam


gatk MarkDuplicatesSpark \
-I ${uniq_run_id}output_${wes_name}.bam \
-O ${uniq_run_id}${wes_name}_marked_duplicates.bam \
--remove-sequencing-duplicates \
--tmp-dir /groups/umcg-weersma/tmp01/Iwan/pipeline/tmp

rm ${uniq_run_id}${wes_name}_recal_data.table
mv ${uniq_run_id}${wes_name}_marked_duplicated.bam ${uniq_run_id}output_${wes_name}.bam

module load SAMtools

samtools sort ${uniq_run_id}output_${fastq_name}.bam -o ${uniq_run_id}${fastq_name}sorted.bam
mv ${uniq_run_id}${fastq_name}sorted.bam /groups/umcg-weersma/tmp01/Iwan/pipeline/output_preprocessing/${uniq_run_id}output_${fastq_name}.bam
samtools index /groups/umcg-weersma/tmp01/Iwan/pipeline/output_preprocessing/${uniq_run_id}output_${fastq_name}.bam


samtools sort ${uniq_run_id}output_${wes_name}.bam -o ${uniq_run_id}${wes_name}sorted.bam
mv ${uniq_run_id}${wes_name}sorted.bam /groups/umcg-weersma/tmp01/Iwan/pipeline/output_preprocessing/${uniq_run_id}output_${wes_name}.bam
samtools index /groups/umcg-weersma/tmp01/Iwan/pipeline/output_preprocessing/${uniq_run_id}output_${wes_name}.bam

sbatch /groups/umcg-weersma/tmp01/Iwan/run_mutect2.sh /groups/umcg-weersma/tmp01/Iwan/pipeline/output_preprocessing/${uniq_run_id}output_${fastq_name}.bam /groups/umcg-weersma/tmp01/Iwan/pipeline/output_preprocessing/${uniq_run_id}output_${wes_name}.bam ${uniq_run_id}
