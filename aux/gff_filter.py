import argparse
import sys
from collections import defaultdict

def parse_attributes(attributes_str):
    """Parses the key-value pairs from the GFF attributes string."""
    attributes = {}
    if attributes_str:
        for part in attributes_str.split(';'):
            if '=' in part:
                key, value = part.split('=', 1)
                attributes[key.strip()] = value.strip()
    return attributes

def filter_and_extract_gene_ids(input_file, output_file, max_intron_density, max_introns_per_gene):
    """Filters a GFF file based on a gene's intron density and total intron count."""
    gene_lengths = {}
    mrna_to_gene_map = {}
    intron_counts = defaultdict(int)

    print(f"Parsing GFF file: '{input_file}' in pass 1...")
    # PASS 1: Collect gene lengths and map mRNA to gene IDs
    try:
        with open(input_file, 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 8:
                    continue
                
                feature_type = parts[2]
                attributes = parse_attributes(parts[8])
                
                if feature_type == 'gene':
                    gene_id = attributes.get('ID')
                    start = int(parts[3])
                    end = int(parts[4])
                    if gene_id:
                        gene_lengths[gene_id] = end - start + 1
                elif feature_type in ['mRNA', 'transcript']:
                    mrna_id = attributes.get('ID')
                    parent_id = attributes.get('Parent')
                    if mrna_id and parent_id:
                        mrna_to_gene_map[mrna_id] = parent_id
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        return

    print("Pass 1 complete. Analyzing intron counts in pass 2...")

    # PASS 2: Collect intron counts and link to genes using the map
    try:
        with open(input_file, 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 8:
                    continue
                
                feature_type = parts[2]
                attributes = parse_attributes(parts[8])
                
                if feature_type == 'intron':
                    parent_id = attributes.get('Parent')
                    if parent_id and parent_id in mrna_to_gene_map:
                        gene_id = mrna_to_gene_map[parent_id]
                        intron_counts[gene_id] += 1
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        return

    print("Pass 2 complete. Applying filter based on intron density and count...")
    
    passing_gene_ids = set()
    total_genes = 0
    passed_filter_count = 0
    
    for gene_id, length in gene_lengths.items():
        total_genes += 1
        
        if length == 0:
            print(f"Warning: Gene '{gene_id}' has a length of 0. Skipping density calculation.")
            continue
        
        # Calculate intron density
        intron_density = (intron_counts[gene_id] / length) * 1000
        
        # Check against both filters
        passes_density_filter = intron_density <= max_intron_density
        passes_count_filter = intron_counts[gene_id] <= max_introns_per_gene
        
        if passes_density_filter and passes_count_filter:
            passing_gene_ids.add(gene_id)
            passed_filter_count += 1
    
    print(f"Filtering complete. Total genes analyzed: {total_genes}")
    print(f"Genes passed filter: {passed_filter_count}")
    print(f"Writing unique gene IDs to output file: '{output_file}'...")
    
    try:
        with open(output_file, 'w') as outfile:
            for gene_id in sorted(list(passing_gene_ids)):
                outfile.write(gene_id + '\n')
        print("Done.")
    except Exception as e:
        print(f"Error: Could not write to output file '{output_file}'. Reason: {e}", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter GFF files based on a gene's intron density and total intron count.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Path to the input GFF file.")
    parser.add_argument('-o', '--output', type=str, required=True, help="Path to the output text file.")
    parser.add_argument('-d', '--density', type=float, default=6.0, help="The maximum number of introns per 1000 bp allowed. (default: 6.0).")
    parser.add_argument('-m', '--max-introns', type=int, default=30, help="The maximum total number of introns per gene. (default: 30).")
 
    args = parser.parse_args()
    
    filter_and_extract_gene_ids(args.input, args.output, args.density, args.max_introns)
