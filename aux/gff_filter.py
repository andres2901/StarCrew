#!/usr/bin/env python3
"""
GFF Intron Density Filter.

This module filters genes in GFF files by analyzing the relationship 
between gene length and intron counts. It extracts gene IDs that fall
below specified thresholds for intron density (introns/kb) and total 
intron count.
"""

import argparse
import sys
from collections import defaultdict


def parse_attributes(attributes_str):
    """
    Parse the key-value pairs from the GFF attributes column (column 9).

    Args:
        attributes_str (str): The semicolon-separated attributes string.

    Returns:
        dict: A dictionary of attribute keys and values.
    """
    attributes = {}
    if attributes_str:
        for part in attributes_str.split(';'):
            if '=' in part:
                key, value = part.split('=', 1)
                attributes[key.strip()] = value.strip()
    return attributes


def filter_and_extract_gene_ids(input_file, output_file, max_density, max_count):
    """
    Filter a GFF file based on intron density and total count.

    Args:
        input_file (str): Path to input GFF.
        output_file (str): Path for the resulting list of Gene IDs.
        max_density (float): Max introns allowed per 1000 bp.
        max_count (int): Max total introns allowed per gene.
    """
    gene_lengths = {}
    mrna_to_gene_map = {}
    intron_counts = defaultdict(int)

    # PASS 1: Collect gene lengths and map mRNA to gene IDs
    print(f"Parsing GFF: '{input_file}' (Pass 1/2)...")
    try:
        with open(input_file, 'r') as infile:
            for line in infile:
                if not line.strip() or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                attrs = parse_attributes(parts[8])

                if feature_type == 'gene':
                    gene_id = attrs.get('ID')
                    if gene_id:
                        # end - start + 1
                        gene_lengths[gene_id] = int(parts[4]) - int(parts[3]) + 1
                elif feature_type in ['mRNA', 'transcript']:
                    mrna_id = attrs.get('ID')
                    parent_id = attrs.get('Parent')
                    if mrna_id and parent_id:
                        mrna_to_gene_map[mrna_id] = parent_id
    except FileNotFoundError:
        sys.exit(f"Error: Input file '{input_file}' not found.")

    # PASS 2: Collect intron counts linked to genes
    print("Analyzing introns (Pass 2/2)...")
    with open(input_file, 'r') as infile:
        for line in infile:
            if not line.strip() or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            if parts[2] == 'intron':
                attrs = parse_attributes(parts[8])
                parent_id = attrs.get('Parent')
                if parent_id in mrna_to_gene_map:
                    gene_id = mrna_to_gene_map[parent_id]
                    intron_counts[gene_id] += 1

    # Apply Filters
    print("Applying filters...")
    passing_gene_ids = []
    for gene_id, length in gene_lengths.items():
        if length == 0:
            continue
        
        # Calculation: (Count / Length) * 1000
        density = (intron_counts[gene_id] / length) * 1000
        
        if density <= max_density and intron_counts[gene_id] <= max_count:
            passing_gene_ids.append(gene_id)

    # Output Results
    print(f"Total genes analyzed: {len(gene_lengths)}")
    print(f"Genes passed filter: {len(passing_gene_ids)}")
    
    try:
        with open(output_file, 'w') as outfile:
            for gene_id in sorted(passing_gene_ids):
                outfile.write(f"{gene_id}\n")
        print(f"Unique IDs written to '{output_file}'. Done.")
    except Exception as e:
        sys.exit(f"Error writing to output: {e}")


def main():
    """CLI entry point for the GFF filter script."""
    parser = argparse.ArgumentParser(
        description="Filter GFF files by gene intron density and count."
    )
    parser.add_argument('-i', '--input', required=True, 
                        help="Path to the input GFF file.")
    parser.add_argument('-o', '--output', required=True, 
                        help="Path to the output text file.")
    parser.add_argument('-d', '--density', type=float, default=6.0, 
                        help="Max introns per 1000 bp (default: 6.0).")
    parser.add_argument('-m', '--max-introns', type=int, default=30, 
                        help="Max total introns per gene (default: 30).")

    args = parser.parse_args()
    filter_and_extract_gene_ids(args.input, args.output, args.density, args.max_introns)


if __name__ == "__main__":
    main()