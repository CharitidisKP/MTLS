#!/usr/bin/env bash

set -euo pipefail

## Directories and Variables ##
SRV=admin4
Exercise_Dir="$HOME/Exercise"
Raw_Dir="/home/${SRV}/SNP_calling/first_part/raw_data"
Data_Dir="/home/${SRV}/2024_yeast_NGS"
Ref_Src_Dir="/home/${SRV}/SNP_calling/first_part/reference"
Picard_Jar="/home/${SRV}/SNP_calling/first_part/software/picard.jar"
Ref_Fasta="S288C_reference_plus_2u.fa"

## Create working directories ##
#mkdir -p "${Exercise_Dir}/Variant_Calling"
#mkdir -p "${Exercise_Dir}/reference_dir"

## Copy in the raw fastq ##
#cp "${Raw_Dir}/control1_R1_sub.fastq"   "${Exercise_Dir}/Variant_Calling/"
#cp "${Data_Dir}/3_S3_R1_001.fastq.gz"   "${Exercise_Dir}/Variant_Calling/"

## Run FastQC ##
cd "${Exercise_Dir}/Variant_Calling"
#gzip -d 3_S3_R1_001.fastq.gz
fastqc 3_S3_R1_001.fastq

## Build HISAT2 index ##
cd "${Exercise_Dir}/reference_dir"
#cp "${Ref_Src_Dir}/${Ref_Fasta}" .
hisat2-build "${Ref_Fasta}" yeast

## Align ##
cd "${Exercise_Dir}/Variant_Calling"
hisat2 -q -x "${Exercise_Dir}/reference_dir/yeast" \
        -U control1_R1.fastq  -S control1_R1.sam
hisat2 -q -x "${Exercise_Dir}/reference_dir/yeast" \
        -U 3_S3_R1_001.fastq    -S 3_S3_R1_001.sam

## Convert to BAM, sort and filter ##
for SAMPLE in control1_R1 3_S3_R1_001; do
  samtools view -bS -h ${SAMPLE}.sam \
    | samtools view -b -F 4 - \
    > ${SAMPLE}_mapped.bam

  # Sort
  gatk SortSam \
      -I ${SAMPLE}_mapped.bam \
      -O ${SAMPLE}_sorted.bam \
      -SORT_ORDER coordinate

  # Mark (and remove) duplicates
  gatk MarkDuplicates \
      -INPUT ${SAMPLE}_sorted.bam \
      -OUTPUT ${SAMPLE}_sorted_marked.bam \
      -METRICS_FILE ${SAMPLE}_metrics.txt \
      -VALIDATION_STRINGENCY LENIENT \
      -REMOVE_DUPLICATES true \
      -CREATE_INDEX true

  # Add read groups
  gatk AddOrReplaceReadGroups \
      -I ${SAMPLE}_sorted_marked.bam \
      -O ${SAMPLE}_final.bam \
      -RGID 4 -RGLB lib1 -RGPL ILLUMINA \
      -RGPU unit1 -RGSM 20

  samtools index ${SAMPLE}_final.bam
done

## Prepare reference ##
cd "${Exercise_Dir}/reference_dir"
java -jar "${Picard_Jar}" CreateSequenceDictionary \
     R=${Ref_Fasta} O=${Ref_Fasta%.fa}.dict
samtools faidx ${Ref_Fasta}

## Variant calling ##
cd "${Exercise_Dir}/Variant_Calling"
for SAMPLE in control1_R1 3_S3_R1_001; do
  gatk HaplotypeCaller \
      -R "${Exercise_Dir}/reference_dir/${Ref_Fasta}" \
      -I ${SAMPLE}_final.bam \
      -O ${SAMPLE}_raw_variants.vcf

  # Split SNPs / INDELs
  gatk SelectVariants \
      -R "${Exercise_Dir}/reference_dir/${Ref_Fasta}" \
      -V ${SAMPLE}_raw_variants.vcf \
      --select-type-to-include SNP \
      -O ${SAMPLE}_raw_snps.vcf

  gatk SelectVariants \
      -R "${Exercise_Dir}/reference_dir/${Ref_Fasta}" \
      -V ${SAMPLE}_raw_variants.vcf \
      --select-type-to-include INDEL \
      -O ${SAMPLE}_raw_indels.vcf

  # Filter
  gatk VariantFiltration \
      -R "${Exercise_Dir}/reference_dir/${Ref_Fasta}" \
      -V ${SAMPLE}_raw_snps.vcf \
      --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' \
      --filter-name basic_snp_filter \
      -O ${SAMPLE}_filtered_snps.vcf

  gatk VariantFiltration \
      -R "${Exercise_Dir}/reference_dir/${Ref_Fasta}" \
      -V ${SAMPLE}_raw_indels.vcf \
      --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' \
      --filter-name basic_indel_filter \
      -O ${SAMPLE}_filtered_indels.vcf

  # Base recalibration
  gatk BaseRecalibrator \
      -R "${Exercise_Dir}/reference_dir/${Ref_Fasta}" \
      -I ${SAMPLE}_final.bam \
      --known-sites ${SAMPLE}_filtered_snps.vcf \
      --known-sites ${SAMPLE}_filtered_indels.vcf \
      -O ${SAMPLE}_recal_data.table

  gatk ApplyBQSR \
      -R "${Exercise_Dir}/reference_dir/${Ref_Fasta}" \
      -I ${SAMPLE}_final.bam \
      --bqsr-recal-file ${SAMPLE}_recal_data.table \
      -O ${SAMPLE}_recal_reads.bam

  # Final variant calling
  gatk HaplotypeCaller \
      -R "${Exercise_Dir}/reference_dir/${Ref_Fasta}" \
      -I ${SAMPLE}_recal_reads.bam \
      -O ${SAMPLE}_variants_recal.vcf
done

## Mental sanity ##
echo "All steps completed."

