#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Normalize entire chrHis chromosome by dividing CPM by 100

This script reads a bigWig file and divides ALL signal on chrHis by 100
to account for the collapsed histone gene array (~100 copies).

All other chromosomes remain unchanged.

Requires: pyBigWig (install with: pip install pyBigWig)

Usage:
    python normalize_chrHis_by100.py input.bigWig output_normalized.bigWig
"""

import sys

# Check for pyBigWig
try:
    import pyBigWig
except ImportError:
    print("ERROR: pyBigWig not installed.", file=sys.stderr)
    print("Install with: pip install pyBigWig", file=sys.stderr)
    sys.exit(1)

NORMALIZATION_FACTOR = 100.0

def normalize_chromosome(bw_in, chrom, chrom_size, normalization_factor):
    """
    Extract and normalize values for entire chromosome
    
    Args:
        bw_in: input bigWig file handle
        chrom: chromosome name
        chrom_size: chromosome length
        normalization_factor: value to divide by (100)
    
    Returns:
        tuple: (starts, ends, values) for normalized data
    """
    # Get all intervals with values
    try:
        intervals_data = bw_in.intervals(chrom)
        if intervals_data is None:
            return None, None, None
    except RuntimeError:
        return None, None, None
    
    starts = []
    ends = []
    values = []
    
    for start, end, value in intervals_data:
        # Divide value by normalization factor
        normalized_value = value / normalization_factor
        
        starts.append(start)
        ends.append(end)
        values.append(normalized_value)
    
    return starts, ends, values

def create_normalized_bigwig(input_bw, output_bw):
    """
    Create normalized bigWig file with chrHis divided by 100
    """
    print(f"\n{'='*70}")
    print(f"Normalizing chrHis chromosome by dividing CPM by {NORMALIZATION_FACTOR}")
    print(f"{'='*70}")
    print(f"Input:  {input_bw}")
    print(f"Output: {output_bw}")
    print(f"\nNormalization:")
    print(f"  chrHis: ALL values / {NORMALIZATION_FACTOR}")
    print(f"  Other chromosomes: unchanged")
    print(f"{'='*70}\n")
    
    # Open input bigWig
    print("Opening input bigWig...")
    bw_in = pyBigWig.open(input_bw)
    
    # Get chromosome info
    chroms = bw_in.chroms()
    print(f"Found {len(chroms)} chromosomes")
    
    # Check if chrHis exists
    if 'chrHis' not in chroms:
        print("\nWARNING: chrHis not found in this bigWig file!")
        print(f"Available chromosomes: {', '.join(sorted(chroms.keys()))}")
        print("\nProceeding anyway - output will be identical to input")
    else:
        print(f"chrHis found: {chroms['chrHis']:,} bp")
    
    # Create output bigWig
    print("\nCreating output bigWig...")
    bw_out = pyBigWig.open(output_bw, "w")
    bw_out.addHeader(list(chroms.items()))
    
    # Process each chromosome
    print("\nProcessing chromosomes:")
    chrHis_processed = False
    
    for chrom, chrom_size in chroms.items():
        print(f"  {chrom:15s} ({chrom_size:12,} bp)...", end=' ')
        
        # Get chromosome data
        try:
            intervals_data = bw_in.intervals(chrom)
            if intervals_data is None:
                print("no data")
                continue
        except RuntimeError:
            print("no data")
            continue
        
        # For chrHis, normalize by dividing by 100
        # For all other chromosomes, keep original values
        if chrom == 'chrHis':
            starts, ends, values = normalize_chromosome(
                bw_in, chrom, chrom_size, NORMALIZATION_FACTOR
            )
            chrHis_processed = True
            print(f"{len(starts):,} intervals (NORMALIZED / {NORMALIZATION_FACTOR})")
        else:
            # Copy values unchanged
            starts = []
            ends = []
            values = []
            for start, end, value in intervals_data:
                starts.append(start)
                ends.append(end)
                values.append(value)
            print(f"{len(starts):,} intervals (unchanged)")
        
        # Write to output
        if starts is not None and len(starts) > 0:
            bw_out.addEntries(
                [chrom] * len(starts),
                starts,
                ends=ends,
                values=values
            )
    
    # Close files
    bw_in.close()
    bw_out.close()
    
    print(f"\n{'='*70}")
    if chrHis_processed:
        print(f"SUCCESS: Normalization complete!")
        print(f"\nWhat was done:")
        print(f"  - chrHis: ALL signal divided by {NORMALIZATION_FACTOR}")
        print(f"  - All other chromosomes: unchanged")
    else:
        print(f"WARNING: chrHis was not found or had no data")
        print(f"Output file created but is identical to input")
    print(f"{'='*70}")
    print(f"\nOutput file: {output_bw}")
    print(f"\nThis normalized bigWig can now be used with deepTools!")
    print(f"{'='*70}\n")

def main():
    if len(sys.argv) != 3:
        print(__doc__)
        print("\nExample:")
        print("  python normalize_chrHis_by100.py input.bigWig output_normalized.bigWig")
        sys.exit(1)
    
    input_bigwig = sys.argv[1]
    output_bigwig = sys.argv[2]
    
    # Validate inputs
    try:
        bw_test = pyBigWig.open(input_bigwig)
        bw_test.close()
    except Exception as e:
        print(f"ERROR: Cannot open input bigWig file: {input_bigwig}", file=sys.stderr)
        print(f"Error details: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Run normalization
    create_normalized_bigwig(input_bigwig, output_bigwig)

if __name__ == '__main__':
    main()