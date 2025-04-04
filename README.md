# genome-assembly
Code for assembly of genome

# Genome Assembly and Analysis Pipeline
**Overview**
This repository contains a comprehensive pipeline for genome assembly, scaffolding, repeat analysis and synteny analysis. The pipeline is designed for processing PacBio HIFI reads and includes steps for assembly improvement using Pore-C data, repeat annotation, and comparative genomics.

The pipeline is modular, allowing users to execute specific steps independently or as a complete workflow.

# Features
* Genome Assembly:

Assembles PacBio reads using hifiasm.
Improves assembly using Pore-C data for scaffolding with SALSA.
* Chromosome Sorting:

Sorts chromosomes based on a reference genome using NGSEP.
* Repeat Analysis:

Identifies and annotates repetitive elements using RepeatModeler and RepeatMasker.
* Synteny Analysis:

Compares assemblies using fastANI and visualizes synteny with pyGenomeViz.
Generates PAF files using minimap2 for further synteny visualization.
* Synteny Visualization:

Creates circular synteny plots using circlize and viridis in R.

# Pipeline Steps
* Step 1: Genome Assembly
Tool: hifiasm
Description: Assembles PacBio reads into contigs.
Output: Assembled genome files in the Assembly directory.
* Step 2: Assembly Improvement
Tool: minimap2, samtools, SALSA
Description: Improves assembly using Pore-C data for scaffolding.
Output: Scaffolds and contig graph files.
* Step 3: Chromosome Sorting
Tool: NGSEP
Description: Sorts chromosomes based on a reference genome.
Output: Sorted genome files.
* Step 4: Repeat Analysis
Tool: RepeatModeler, RepeatMasker
Description: Identifies and annotates repetitive elements in the genome.
Output: Repeat annotation files.
* Step 5: Synteny Analysis
Tool: fastANI, pyGenomeViz, minimap2
Description: Compares assemblies and visualizes synteny using FastANI and pyGenomeViz.
Output: FastANI results, visualizations, and PAF files.
* Step 6: Synteny Visualization
Tool: circlize, viridis (R)
Description: Generates circular synteny plots from PAF files.
Output: Circular synteny plots.

