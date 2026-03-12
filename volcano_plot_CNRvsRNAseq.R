#!/usr/bin/env Rscript
# ============================================================================
# Volcano Plot: Alternative version for TSV/CSV input
# ============================================================================
# Use this version if your data is in tab-separated (.tsv) or comma-separated
# (.csv) format instead of Excel (.xlsx)
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)     # For plotting
  library(dplyr)       # For data manipulation
  library(ggrepel)     # For non-overlapping labels
  library(svglite)     # For SVG export
  library(ggrastr)     # For rasterizing points (keeps vectors editable in Illustrator)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input file path - UPDATE THIS
input_file <- "/path/to/input/CNRvsRNAseq.csv"

# File format - specify delimiter
file_delimiter <- ","  # Use "\t" for TSV, "," for CSV

# Significance threshold
padj_threshold <- 0.05

# Histone genes
histone_genes <- c("H1.A", "H2a.A", "H2b.A", "H3.A", "H4.A")

# Display names
display_names <- c(
  "H1.A" = "H1",
  "H2a.A" = "H2a",
  "H2b.A" = "H2b",
  "H3.A" = "H3",
  "H4.A" = "H4"
)

# Colors
histone_color <- "#EF426F"
deg_color <- "#00A5AD"

# Number of top DEGs to label
n_top_degs <- 10

# Output files
output_png <- "mute_volcano_plot_cutandrun.png"
output_svg <- "mute_volcano_plot_cutandrun.svg"
output_csv <- "mute_labeled_genes_summary.csv"

# ============================================================================
# READ DATA
# ============================================================================

cat("Reading data file...\n")
cat(paste("Looking for:", input_file, "\n"))

# Check if file exists
if (!file.exists(input_file)) {
  cat("\nERROR: File not found!\n")
  cat("Please verify the file path is correct.\n")
  cat("Current working directory:", getwd(), "\n")
  stop("Cannot proceed without input file.")
}

cat("File found! Reading...\n")

# Read with robust parameters
data <- read.table(
  input_file,
  header = TRUE,
  sep = file_delimiter,
  stringsAsFactors = FALSE,
  quote = "",           # Don't treat quotes specially
  comment.char = "",    # Don't treat # as comment
  fill = TRUE           # Fill incomplete rows
)

cat(paste("Read", nrow(data), "rows\n"))
cat("Columns:", paste(names(data), collapse = ", "), "\n\n")

# ============================================================================
# FILTER DATA
# ============================================================================

cat("\nPreparing data...\n")

# Handle column names with spaces (R converts spaces to dots)
# Check for various possible column name formats
if ("Sum.Signal" %in% names(data) && !("Sum Signal" %in% names(data))) {
  data$`Sum Signal` <- data$Sum.Signal
  cat("Note: Renamed 'Sum.Signal' column to 'Sum Signal'\n")
}

# Also check for log2FoldChange variants
if ("log2FoldChange" %in% names(data)) {
  # Good, already exists
} else if ("log2.Fold.Change" %in% names(data)) {
  data$log2FoldChange <- data$log2.Fold.Change
  cat("Note: Renamed 'log2.Fold.Change' to 'log2FoldChange'\n")
}

# Print available columns
cat("\nAvailable columns:\n")
print(names(data))

# Check for required columns
required_cols <- c("gene_symbol", "Sum Signal", "log2FoldChange", "padj")
missing_cols <- setdiff(required_cols, names(data))

if (length(missing_cols) > 0) {
  cat("\nERROR: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
  cat("Please ensure your data has these columns:\n")
  cat("  - gene_symbol\n")
  cat("  - Sum Signal (or Sum.Signal)\n")
  cat("  - log2FoldChange\n")
  cat("  - padj\n")
  stop("Missing required columns")
}

# Filter and categorize
data_filtered <- data %>%
  filter(!is.na(padj) & !is.na(log2FoldChange) & !is.na(`Sum Signal`)) %>%
  filter(padj < padj_threshold | gene_symbol %in% histone_genes) %>%
  mutate(
    gene_category = case_when(
      gene_symbol %in% histone_genes ~ "Histone",
      TRUE ~ "Non-histone DEG"
    ),
    display_name = ifelse(
      gene_symbol %in% names(display_names),
      display_names[gene_symbol],
      gene_symbol
    )
  )

cat(paste("Filtered to", nrow(data_filtered), "genes\n"))
print(table(data_filtered$gene_category))

# ============================================================================
# IDENTIFY LABELS
# ============================================================================

genes_to_label <- bind_rows(
  data_filtered %>% filter(gene_category == "Histone"),
  data_filtered %>% filter(gene_category == "Non-histone DEG" & log2FoldChange > 0) %>% 
    arrange(padj) %>% head(n_top_degs),
  data_filtered %>% filter(gene_category == "Non-histone DEG" & log2FoldChange < 0) %>% 
    arrange(padj) %>% head(n_top_degs)
)

cat(paste("\nLabeling", nrow(genes_to_label), "genes\n"))

# ============================================================================
# CREATE PLOT
# ============================================================================

# Add column to identify which genes get black outlines (histones + top 10 DEGs)
data_filtered <- data_filtered %>%
  mutate(highlight = gene_symbol %in% genes_to_label$gene_symbol)

# Separate data for outlined vs regular points
data_outlined <- data_filtered %>% filter(highlight)
data_regular <- data_filtered %>% filter(!highlight)

# Add label formatting info to genes_to_label
genes_to_label <- genes_to_label %>%
  mutate(
    label_face = ifelse(gene_category == "Histone", "bold", "plain")
  )

p <- ggplot(data_filtered, aes(x = `Sum Signal`, y = log2FoldChange)) +
  
  # Plot regular points (no outline) - RASTERIZED
  geom_point_rast(data = data_regular, aes(color = gene_category), 
                  alpha = 0.6, size = 2.5, raster.dpi = 300) +
  
  # Plot outlined points (histones + top DEGs) - RASTERIZED
  geom_point_rast(data = data_outlined, aes(fill = gene_category), 
                  shape = 21, size = 2.5, color = "black", stroke = 0.8, 
                  alpha = 0.6, raster.dpi = 300) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray60", linewidth = 0.5) +
  
  # Color scheme - need both color and fill scales
  scale_color_manual(
    values = c("Histone" = histone_color, "Non-histone DEG" = deg_color),
    name = ""
  ) +
  scale_fill_manual(
    values = c("Histone" = histone_color, "Non-histone DEG" = deg_color),
    name = ""
  ) +
  
  # Add text labels - histones in bold, DEGs in plain
  geom_text_repel(
    data = filter(genes_to_label, gene_category == "Histone"),
    aes(label = display_name),
    size = 3,
    fontface = "bold",
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    min.segment.length = 0,
    seed = 42
  ) +
  geom_text_repel(
    data = filter(genes_to_label, gene_category != "Histone"),
    aes(label = display_name),
    size = 3,
    fontface = "plain",
    max.overlaps = 50,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    min.segment.length = 0,
    seed = 42
  ) +
  
  labs(x = "Sum Signal (CPM)", y = "Log2FC expression") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# ============================================================================
# SAVE
# ============================================================================

ggsave(output_png, p, width = 10, height = 8, dpi = 300)
ggsave(output_svg, p, width = 10, height = 8)

genes_to_label %>%
  select(gene_symbol, display_name, `Sum Signal`, log2FoldChange, padj, gene_category) %>%
  write.csv(output_csv, row.names = FALSE)

cat("\nDone! Files saved:\n")
cat(paste("  -", output_png, "\n"))
cat(paste("  -", output_svg, "\n"))
cat(paste("  -", output_csv, "\n"))

cat("\n=== Summary Statistics ===\n")
cat(paste("Total genes plotted:", nrow(data_filtered), "\n"))
cat(paste("  - Histone genes:", sum(data_filtered$gene_category == "Histone"), "\n"))
cat(paste("  - Non-histone DEGs:", sum(data_filtered$gene_category == "Non-histone DEG"), "\n"))
cat(paste("\nGenes labeled:", nrow(genes_to_label), "\n"))
cat("\nLog2FC range:", round(min(data_filtered$log2FoldChange), 2), "to", 
    round(max(data_filtered$log2FoldChange), 2), "\n")
cat("Sum Signal range:", round(min(data_filtered$`Sum Signal`), 2), "to", 
    round(max(data_filtered$`Sum Signal`), 2), "\n")
