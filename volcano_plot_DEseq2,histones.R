#!/usr/bin/env Rscript

# Volcano Plot for DESeq2 Results
# Labels: Top 10 up/down DEGs (colored) + Histone genes (bold, circled)

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)  # for better label positioning
library(ggrastr)  # for rasterizing points (keeps vectors editable in Illustrator)

# Check if svglite is installed (needed for SVG output)
if (!requireNamespace("svglite", quietly = TRUE)) {
  cat("Installing svglite package for SVG output...\n")
  # Install to user's personal library
  personal_lib <- Sys.getenv("R_LIBS_USER")
  if (!dir.exists(personal_lib)) {
    dir.create(personal_lib, recursive = TRUE)
  }
  install.packages("svglite", repos = "https://cloud.r-project.org", lib = personal_lib)
  # Add the personal library to the search path
  .libPaths(c(personal_lib, .libPaths()))
}
library(svglite)

# Set file path
data_file <- "/path/to/input/deseq2.results.tsv"

# Read the DESeq2 results with more robust settings
deseq_results <- read.table(data_file, 
                            header = TRUE, 
                            sep = "\t", 
                            stringsAsFactors = FALSE,
                            quote = "",  # Don't treat quotes specially
                            comment.char = "",  # Don't treat # as comments
                            fill = TRUE)  # Fill incomplete rows

# Display basic information about the data
cat("Data dimensions:", dim(deseq_results), "\n")
cat("Column names:", colnames(deseq_results), "\n")
cat("\nFirst few rows:\n")
print(head(deseq_results))

# Filter out rows with NA values in padj (these can't be plotted meaningfully)
deseq_filtered <- deseq_results %>%
  filter(!is.na(padj) & !is.na(log2FoldChange))

cat("\nRows after filtering NAs:", nrow(deseq_filtered), "\n")

# Define significance thresholds
padj_threshold <- 0.05
log2fc_threshold <- 1  # 2-fold change cutoff

# Add a column to categorize genes
deseq_filtered <- deseq_filtered %>%
  mutate(
    gene_type = case_when(
      padj < padj_threshold & log2FoldChange > log2fc_threshold ~ "Up-regulated",
      padj < padj_threshold & log2FoldChange < -log2fc_threshold ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

# Count genes in each category
cat("\nGene categories:\n")
print(table(deseq_filtered$gene_type))

# Define histone genes to label (original names in data)
histone_genes <- c("H2a.A", "H2b.A", "H1.A", "H3.A", "H4.A")

# Create display name mapping (remove .A suffix)
histone_display_names <- c(
  "H2a.A" = "H2a",
  "H2b.A" = "H2b", 
  "H1.A" = "H1",
  "H3.A" = "H3",
  "H4.A" = "H4"
)

# Get top 10 up-regulated genes
top_up <- deseq_filtered %>%
  filter(gene_type == "Up-regulated") %>%
  arrange(padj) %>%
  head(10)

# Get top 10 down-regulated genes
top_down <- deseq_filtered %>%
  filter(gene_type == "Down-regulated") %>%
  arrange(padj) %>%
  head(10)

# Combine top DEGs
top_degs <- bind_rows(top_up, top_down)

# Get histone genes
histone_data <- deseq_filtered %>%
  filter(gene_symbol %in% histone_genes)

# Create labeling dataset
# Combine top DEGs and histones, avoiding duplicates
all_labels <- deseq_filtered %>%
  filter(gene_symbol %in% c(top_degs$gene_symbol, histone_genes)) %>%
  mutate(
    is_histone = gene_symbol %in% histone_genes,
    is_top_deg = gene_symbol %in% top_degs$gene_symbol,
    display_name = ifelse(is_histone, 
                         histone_display_names[gene_symbol],
                         gene_symbol),
    label_color = case_when(
      is_histone & !is_top_deg ~ "black",  # Histone only: black
      is_histone & is_top_deg & gene_type == "Up-regulated" ~ "#00594C",  # Histone + up DEG: red
      is_histone & is_top_deg & gene_type == "Down-regulated" ~ "#4B9CD3",  # Histone + down DEG: blue
      gene_type == "Up-regulated" ~ "#00594C",  # Top up DEG: red
      gene_type == "Down-regulated" ~ "#4B9CD3",  # Top down DEG: blue
      TRUE ~ "black"
    ),
    label_face = ifelse(is_histone, "bold", "plain")
  )

cat("\nGenes to be labeled:\n")
print(all_labels %>% select(gene_symbol, display_name, is_histone, is_top_deg, gene_type))

cat("\nHistone genes found:\n")
print(histone_data %>% select(gene_symbol, log2FoldChange, padj, gene_type))

# Create base volcano plot - RASTERIZED points
p <- ggplot(deseq_filtered, aes(x = log2FoldChange, y = -log10(padj), color = gene_type)) +
  geom_point_rast(alpha = 0.6, size = 1.5, raster.dpi = 300) +
  scale_color_manual(values = c("Up-regulated" = "#00594C",
                                "Down-regulated" = "#4B9CD3", 
                                "Not significant" = "grey")) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), 
             linetype = "dashed", 
             color = "black", 
             alpha = 0.5) +
  geom_hline(yintercept = -log10(padj_threshold), 
             linetype = "dashed", 
             color = "black", 
             alpha = 0.5) +
  labs(
    title = "Volcano Plot: Mute1281 vs WT",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Gene Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Add circles around histone genes (keep as vector, only 5 genes)
p_labeled <- p +
  geom_point(data = histone_data,
             aes(x = log2FoldChange, y = -log10(padj)),
             color = "black",
             size = 3,
             shape = 21,
             fill = NA,
             stroke = 1.5,
             inherit.aes = FALSE)

# Add labels with custom colors and fonts
# We need to add labels one group at a time to control formatting

# Add all labels together using geom_text_repel
for(i in 1:nrow(all_labels)) {
  gene_data <- all_labels[i, ]
  
  p_labeled <- p_labeled +
    geom_text_repel(
      data = gene_data,
      aes(x = log2FoldChange, y = -log10(padj), label = display_name),
      size = 4,
      fontface = gene_data$label_face,
      color = gene_data$label_color,
      max.overlaps = 30,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "grey50",
      segment.size = 0.3,
      min.segment.length = 0,
      inherit.aes = FALSE
    )
}

# Save the combined plot as PNG
ggsave("volcano_plot_combined_labels.png", 
       plot = p_labeled, 
       width = 12, 
       height = 9, 
       dpi = 300)

# Save the combined plot as SVG for Illustrator
ggsave("volcano_plot_combined_labels.svg", 
       plot = p_labeled, 
       width = 12, 
       height = 9)

cat("\nCombined volcano plot saved as:\n")
cat("  - volcano_plot_combined_labels.png (for viewing)\n")
cat("  - volcano_plot_combined_labels.svg (for Illustrator)\n")

# Create summary table
summary_table <- all_labels %>%
  select(gene_symbol, display_name, baseMean, log2FoldChange, padj, gene_type, is_histone, is_top_deg) %>%
  arrange(padj)

write.csv(summary_table, "labeled_genes_summary.csv", row.names = FALSE)
cat("Labeled genes summary saved as: labeled_genes_summary.csv\n")

# Print summary
cat("\nTop 10 up-regulated genes:\n")
print(top_up %>% select(gene_symbol, log2FoldChange, padj))

cat("\nTop 10 down-regulated genes:\n")
print(top_down %>% select(gene_symbol, log2FoldChange, padj))

cat("\nHistone genes in dataset:\n")
print(histone_data %>% select(gene_symbol, log2FoldChange, padj, gene_type))

cat("\nScript completed successfully!\n")
