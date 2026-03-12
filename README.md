# Geisler2026
All code associated with the Geisler 2026 publication "Differential control of both cell cycle-regulated and quantitative histone mRNA expression by Drosophila Mute"

normalize_chrHis_by100.py 
- This script reads a bigWig file and divides ALL signal on chrHis by 100 to account for the collapsed histone gene array (~100 copies). All other chromosomes remain unchanged.
- Dependencies: pyBigWig (install with: pip install pyBigWig)
- Usage: python normalize_chrHis_by100.py input.bigWig output_normalized.bigWig

visualize_genomic_bins.py
- Genome-wide signal visualization comparing Mute vs MXC CUT&RUN samples, with special emphasis on the chrHis locus. Produces multiple manhanttan and MA plots.
- Input: A TSV file from multiBigwigSummary (binned genome-wide signal across all samples), hardcoded at the top of the script as filepath
- Dependencies: pandas, matplotlib, numpy

volcano_plot_DEseq2,histones.R
- generates volcano plot of -log10(padj) vs Log2(fc) at all genes from DEseq2 output, highlighting histone genes and the top ten up and down DEGs.
- Dependencies: R

extract_matrix_stats.py
- Parses deepTools computeMatrix output and summarizes per-gene CUT&RUN signal statistics into an Excel file.
- Input: A gzipped matrix file (<sample>_gene_matrix.gz) produced by computeMatrix scale-regions
- Output: An Excel file (<sample>_signal_stats.xlsx) with one row per gene containing: Chr, Start, End, Gene (FBgn ID), Score, Strand, and five signal summary stats — Mean, Median, Max, Min, and Sum Signal.
- Dependencies: Python 3, numpy, pandas, openpyxl
- Usage: python extract_matrix_stats.py <matrix.gz> <output.xlsx> 

merge_deseq2+CNR.py
- Merges the per-sample signal stats Excel files with DESeq2 differential expression results, matching genes by FBgn ID.
- Input: (1) Per-sample *_signal_stats.xlsx files from extract_matrix_stats.py; (2) a tab-separated DESeq2 results file (.tsv) containing columns including gene_id, log2FoldChange, padj, gene_symbol, etc.
- Output: Per-sample *_signal_deseq2_merged.xlsx files 
- Dependencies: Python 3, pandas, openpyxl

volcano_plot_CNRvsRNAseq.R
- Generates volcano plot of Log2FC expression vs Sum Mute signal at all singificant genes from DEseq2 output as well as histone genes.
- Dependencies: R

Raw data available on NCBI GEO: GSE324623 (Cut&Run), GSE324624 (RNA-seq) upon publication
