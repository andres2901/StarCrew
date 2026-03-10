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


def parse_attributes(attributes_str: str) -> dict:
    """Parse key-value pairs from the GFF attributes column (column 9).

    Args:
        attributes_str: Semicolon-separated attributes string from
                        column 9 of a GFF file.

    Returns:
        Dictionary of attribute keys mapped to their values.
    """

    attributes = {}
    if attributes_str:
        for part in attributes_str.split(';'):
            if '=' in part:
                key, value = part.split('=', 1)
                attributes[key.strip()] = value.strip()
    return attributes


def filter_and_extract_gene_ids(
    input_file: str,
    output_file: str,
    max_density: float,
    max_count: int
) -> None:
    """Filter a GFF file based on intron density and total intron count.

    Performs two passes over the input file: the first collects gene
    lengths and maps mRNA IDs to their parent gene IDs; the second
    counts introns per gene. Genes passing both thresholds are written
    to the output file as a sorted list of IDs.

    Args:
        input_file: Path to the input GFF file.
        output_file: Path for the output file containing filtered gene IDs.
        max_density: Maximum introns allowed per 1000 bp.
        max_count: Maximum total introns allowed per gene.

    Raises:
        FileNotFoundError: If the input file does not exist.
        PermissionError: If the output file cannot be written.
        OSError: If any other I/O error occurs during writing.
    """

    gene_lengths = {}
    mrna_to_gene_map = {}
    intron_counts = defaultdict(int)

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

    try:
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
    except FileNotFoundError:
        sys.exit(f"Error: Input file '{input_file}' not found.")

    passing_gene_ids = []
    for gene_id, length in gene_lengths.items():
        if length == 0:
            continue
        
        density = (intron_counts[gene_id] / length) * 1000
        
        if density <= max_density and intron_counts[gene_id] <= max_count:
            passing_gene_ids.append(gene_id)

    print(f"Total genes analyzed: {len(gene_lengths)}")
    print(f"Genes passed filter: {len(passing_gene_ids)}")
    
    try:
        with open(output_file, 'w') as outfile:
            for gene_id in sorted(passing_gene_ids):
                outfile.write(f"{gene_id}\n")
        print(f"Unique IDs written to '{output_file}'. Done.")
    except PermissionError:
        sys.exit(f"Error: No write permission for: {output_file}")
    except OSError as e:
        sys.exit(f"Error writing output file: {e}")


def main() -> None:
    """Parse command-line arguments and launch the GFF intron density filter."""

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