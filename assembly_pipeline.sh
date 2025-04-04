#!/bin/bash
# Script: assembly_pipeline_with_repeats.sh
# Description:
# This script performs genome assembly, scaffolding, and repeat analysis using PacBio reads and Pore-C data.
# Author: Lissa Cruz-Saavedra
# Date: 04-04-2025
# Version: 1.2

# Exit immediately if a command exits with a non-zero status
set -e

# Load necessary modules
module load hifiasm/0.15.2
module load bedtools/2.30.0
module load samtools/1.11
module load repeatmodeler/2.0.1
module load repeatmasker/4.1.2
module load python/3.8

# Set paths (update these as needed)
RAW_DATA_DIR="/path/to/raw_data"                # Directory containing raw data
ASSEMBLY_DIR="${RAW_DATA_DIR}/Assembly"         # Directory for assembly results
POREC_DIR="${ASSEMBLY_DIR}/PoreC"               # Directory for Pore-C data
OUTPUT_DIR="/path/to/output"                    # Output directory
SAMPLE_NAME="<SAMPLE_NAME>"                     # Replace <SAMPLE_NAME> with your sample name

# Create necessary directories
mkdir -p "$ASSEMBLY_DIR" "$POREC_DIR" "$OUTPUT_DIR"

### Step 1: PacBio Assembly ###
echo "Step 1: Performing PacBio reads assembly..."
hifiasm -o ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm -t 48 ${RAW_DATA_DIR}/${SAMPLE_NAME}_PacBio_data.fa

### Step 2: Assembly Improvement ###
echo "Step 2: Improving assembly with Pore-C data..."
minimap2 -ay ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg.fa ${RAW_DATA_DIR}/${SAMPLE_NAME}_PoreC.fastq.gz > ${POREC_DIR}/${SAMPLE_NAME}_aln.sam
samtools view -bS ${POREC_DIR}/${SAMPLE_NAME}_aln.sam > ${POREC_DIR}/${SAMPLE_NAME}_poreC_mappedreads.bam
bamToBed -i ${POREC_DIR}/${SAMPLE_NAME}_poreC_mappedreads.bam > ${POREC_DIR}/${SAMPLE_NAME}_poreC_mappedreads.bed
sort -k 4 ${POREC_DIR}/${SAMPLE_NAME}_poreC_mappedreads.bed > ${POREC_DIR}/${SAMPLE_NAME}_poreC_sorted.bed
awk '{$4 = $4 "/1"; print}' ${POREC_DIR}/${SAMPLE_NAME}_poreC_sorted.bed > ${POREC_DIR}/${SAMPLE_NAME}_poreC_final.bed

# Generate contig lengths file
samtools faidx ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg.fa

# Run SALSA for scaffolding
python2.7 /path/to/SALSA/run_pipeline.py \
  -a ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg.fa \
  -l ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg.fa.fai \
  -b ${POREC_DIR}/${SAMPLE_NAME}_poreC_final.bed \
  -e AAGCTT \
  -o ${ASSEMBLY_DIR}/${SAMPLE_NAME}_scaffolds \
  -m yes \
  -g ${ASSEMBLY_DIR}/${SAMPLE_NAME}_scaffolds/${SAMPLE_NAME}_contigs_graph.gfa

### Step 3: Chromosome Sorting ###
echo "Step 3: Sorting chromosomes based on reference genome..."
java -jar /path/to/NGSEPcore.jar AssemblyReferenceSorter \
  -i ${ASSEMBLY_DIR}/${SAMPLE_NAME}_genome.final.fasta \
  -r ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg.fa \
  -o ${ASSEMBLY_DIR}/${SAMPLE_NAME}_genome_sorted.fasta -rcp 2

### Step 4: Repeat Analysis ###
echo "Step 4: Performing repeat analysis..."
BuildDatabase -name ${OUTPUT_DIR}/${SAMPLE_NAME}_repeat_db ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg.fa
RepeatModeler -database ${OUTPUT_DIR}/${SAMPLE_NAME}_repeat_db -pa 36 -LTRStruct > ${OUTPUT_DIR}/${SAMPLE_NAME}_repeat_modeler.log
RepeatMasker -pa 36 -gff -lib consensi-refined.fa -dir ${OUTPUT_DIR}/${SAMPLE_NAME}_repeat_masker ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg.fa

echo "Assembly, scaffolding, and repeat analysis completed successfully."
