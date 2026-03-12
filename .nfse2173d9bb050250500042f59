#!/usr/bin/env python3
"""
Merge deepTools signal stats Excel files with DESeq2 expression data.
Matches on Gene (FBgn ID) in signal stats to gene_id in DESeq2 results.
"""

import subprocess
import sys

try:
    import openpyxl
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", "openpyxl"])
    import openpyxl

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment

DESEQ2_FILE = "/path/to/input/deseq2.results.tsv"

SIGNAL_DIR = "/path/to/input/matrix_output_repMergeHis100"
OUTPUT_DIR = SIGNAL_DIR

SAMPLE_NAMES = [
    "sample1",
    "sample2",
    
]

DESEQ2_COLS = ["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "gene_symbol"]

def merge_and_write(signal_xlsx, deseq2_df, output_xlsx):
    # Read signal stats
    signal_df = pd.read_excel(signal_xlsx)
    print(f"  Signal stats: {len(signal_df)} genes")

    # Left join on Gene == gene_id so all signal stats genes are kept
    merged = signal_df.merge(
        deseq2_df,
        left_on="Gene",
        right_on="gene_id",
        how="left"
    )
    # Drop the redundant gene_id column after merge
    merged.drop(columns=["gene_id"], inplace=True)

    matched = merged[DESEQ2_COLS[0]].notna().sum()
    print(f"  Matched to DESeq2: {matched} genes")
    print(f"  Unmatched: {len(merged) - matched} genes")

    # Load existing workbook to preserve formatting style, then write fresh
    wb = load_workbook(signal_xlsx)
    ws = wb.active

    # Signal stat headers already exist in cols 1-11; append DESeq2 headers starting at col 12
    deseq2_header_fill = PatternFill(start_color="70AD47", end_color="70AD47", fill_type="solid")
    header_font = Font(bold=True, color="FFFFFF")
    header_align = Alignment(horizontal="center")

    for col_offset, col_name in enumerate(DESEQ2_COLS):
        col_idx = 12 + col_offset
        cell = ws.cell(row=1, column=col_idx, value=col_name)
        cell.font = header_font
        cell.fill = deseq2_header_fill
        cell.alignment = header_align
        # Set column widths
        ws.column_dimensions[cell.column_letter].width = 16

    # Write merged DESeq2 data rows
    for row_idx in range(len(merged)):
        excel_row = row_idx + 2  # +2 for 1-indexed + header
        for col_offset, col_name in enumerate(DESEQ2_COLS):
            col_idx = 12 + col_offset
            val = merged.iloc[row_idx][col_name]
            # Write None for NaN so Excel shows blank instead of nan
            cell = ws.cell(row=excel_row, column=col_idx, value=None if pd.isna(val) else val)
            # Format numeric columns
            if col_name in ("pvalue", "padj"):
                cell.number_format = "0.00E+00"
            elif col_name not in ("gene_symbol",) and pd.api.types.is_numeric_dtype(merged[col_name]):
                cell.number_format = "0.0000"

    wb.save(output_xlsx)
    print(f"  Output: {output_xlsx}\n")

def main():
    # Read DESeq2 file once
    deseq2_df = pd.read_csv(DESEQ2_FILE, sep="\t")
    print(f"DESeq2 file loaded: {len(deseq2_df)} genes")
    print(f"DESeq2 columns: {list(deseq2_df.columns)}\n")

    # Dynamically set DESEQ2_COLS to everything except gene_id (which is used for merging)
    global DESEQ2_COLS
    DESEQ2_COLS = [c for c in deseq2_df.columns if c != "gene_id"]

    for sample in SAMPLE_NAMES:
        signal_xlsx = f"{SIGNAL_DIR}/{sample}_signal_stats.xlsx"
        output_xlsx = f"{OUTPUT_DIR}/{sample}_signal_deseq2_merged.xlsx"
        print(f"Processing: {sample}")
        merge_and_write(signal_xlsx, deseq2_df, output_xlsx)

    print("All samples complete.")

if __name__ == "__main__":
    main()
