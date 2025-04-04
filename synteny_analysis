#!/bin/bash
# Script: synteny_analysis.sh
# Description:
# This script performs synteny analysis and visualization using fastANI, pyGenomeViz, and minimap2.
# Author: Lissa Cruz-Saavedra
# Date: 04-04-2025
# Version: 1.1

# Exit immediately if a command exits with a non-zero status
set -e

# Load necessary modules
module load fastani/1.33
module load python/3.8
module load minimap2/2.24
module load R/4.2.0

# Set paths (update these as needed)
ASSEMBLY_DIR="/path/to/Assembly"  # Directory containing assembly results
GENOME_FA="/path/to/genome.fa"    # Reference genome file
SAMPLE_NAME="<SAMPLE_NAME>"       # Replace <SAMPLE_NAME> with your sample name

# Create necessary directories if needed
mkdir -p "$ASSEMBLY_DIR"

### Step 6: Synteny Analysis ###
echo "Step 6: Performing synteny analysis..."

# Run fastANI for synteny analysis
fastANI -q ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg.1000000.fa \
  -r ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg_order.fa \
  --visualize -o ${ASSEMBLY_DIR}/${SAMPLE_NAME}_fastANI_visualized.out -t 8

# Visualize synteny using pyGenomeViz
python /path/to/pyGenomeViz/notebooks/fastANI/visualize.py \
  ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg.1000000.fa ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_ctg.less1000000.fa \
  ${ASSEMBLY_DIR}/${SAMPLE_NAME}_fastANI_visualized.out.visual \
  ${ASSEMBLY_DIR}/${SAMPLE_NAME}_fastANI_visualized.png

# Generate PAF file using minimap2
minimap2 ${ASSEMBLY_DIR}/${SAMPLE_NAME}.asm.bp.p_cta.edit.sorted.fasta ${GENOME_FA} > ${ASSEMBLY_DIR}/${SAMPLE_NAME}_vs_reference.paf

### Step 7: Synteny Visualization ###
echo "Step 7: Performing synteny visualization..."
Rscript -e "
library(tidyverse)
library(circlize)
library(viridis)

# Read the PAF file
paf_file <- read_tsv('${ASSEMBLY_DIR}/${SAMPLE_NAME}_vs_reference.paf', col_names = FALSE, show_col_types = FALSE)

# Assign column names
colnames(paf_file) <- c('query', 'qlen', 'qstart', 'qend', 'strand', 'target', 'tlen', 'tstart', 'tend', 'mapq')

# Select necessary columns for synteny
synteny_data <- paf_file %>%
  select(query, qstart, qend, target, tstart, tend, strand)

# Convert to a format for Circlize
links <- synteny_data %>%
  mutate(qstart = as.numeric(qstart), qend = as.numeric(qend),
         tstart = as.numeric(tstart), tend = as.numeric(tend))

# Define the order of chromosomes
genome_1_order <- c(paste0('T', 1:22), 'T22b', paste0('T', 23:31))
genome_2_order <- c(paste0('c', 1:100), 43:1)
final_order <- c(genome_2_order, genome_1_order)

# Filter links to include only valid sectors
links <- links %>%
  filter(query %in% final_order & target %in% final_order)

# Generate unique colors for genome_1 chromosomes
genome_1_colors <- setNames(viridis(32), genome_1_order)

# Assign colors to chromosomes
sector_colors <- sapply(final_order, function(sector) {
  if (sector %in% genome_1_order) {
    genome_1_colors[sector]
  } else {
    'gray'
  }
})

# Calculate xlim for each sector
xlim <- data.frame(
  sector = final_order,
  start = 0,
  end = sapply(final_order, function(sector) {
    if (sector %in% links$query | sector %in% links$target) {
      max(c(
        max(links$qend[links$query == sector], na.rm = TRUE),
        max(links$tend[links$target == sector], na.rm = TRUE)
      ), na.rm = TRUE)
    } else {
      10
    }
  })
)

# Remove sectors with zero or negative width
xlim <- xlim[xlim$end > 0, ]

# Initialize Circlize
ordered_sectors <- factor(xlim$sector, levels = final_order)
circos.par(cell.padding = c(0.01, 0, 0.01, 0))
circos.initialize(factors = ordered_sectors, xlim = xlim[, c('start', 'end')])

# Add track regions
circos.trackPlotRegion(factors = ordered_sectors, ylim = c(0, 1), panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data('sector.index')
  circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = sector_colors[sector.name], border = NA)
  circos.text(CELL_META$xcenter, 0.5, sector.name, cex = 0.8, facing = 'outside', niceFacing = TRUE)
})

# Add synteny links
for (i in 1:nrow(links)) {
  link_color <- if (links$query[i] %in% genome_1_order) {
    genome_1_colors[links$query[i]]
  } else if (links$target[i] %in% genome_1_order) {
    genome_1_colors[links$target[i]]
  } else {
    'gray'
  }
  circos.link(links$query[i], c(links$qstart[i], links$qend[i]),
              links$target[i], c(links$tstart[i], links$tend[i]),
              col = link_color)
}

# Clear Circlize
circos.clear()
"

echo "Synteny analysis and visualization completed successfully."
