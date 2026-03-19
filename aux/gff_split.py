#!/usr/bin/env python3
"""
Split GFF file base on contigs

This script indexes a master GFF into memory and split it based on contigID
from the fasta file that it comes from.
"""

import argparse
import sys
import os


def load_gff_features(gff_path: str, fasta_path: str) -> tuple[dict, list]:
    """Read and index a GFF file into memory by contig ID.

    Identify the seqID of the contigs from the fasta file that will be loaded.
    Performs a single pass over the GFF file, separating header lines from
    feature lines and grouping features by their SeqID (column 1).

    Args:
        gff_path: Path to the input GFF file.
        fasta_path: Path to the input fasta file.

    Returns:
        Tuple (indexed_features, header_lines) where:
            - indexed_features: Dictionary mapping contig ID to a list
                                of feature rows (each row is a list of
                                9 string columns).
            - header_lines: List of comment/header lines preserved
                            for output.

    Raises:
        FileNotFoundError: If the file does not exist.
        PermissionError: If the file cannot be read.
        OSError: If any other I/O error occurs.
    """
    
    indexed_features = {}
    header_lines = []
    
    try:
        with open(fasta_path, 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    indexed_features[line.replace(">","")] = []

    except FileNotFoundError:
        sys.exit(f"Error: Fasta file '{fasta_path}' not found.")
    except PermissionError:
        sys.exit(f"Error: No read permission for '{fasta_path}'.")
    except OSError as e:
        sys.exit(f"Error reading Fasta file: {e}")
    
    try:
        with open(gff_path, 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('#'):
                    header_lines.append(line + '\n')
                    continue
                
                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                f_seqid = parts[0]
                if f_seqid not in indexed_features:
                    indexed_features[f_seqid] = []
                    
                indexed_features[f_seqid].append(parts)

    except FileNotFoundError:
        sys.exit(f"Error: GFF file '{gff_path}' not found.")
    except PermissionError:
        sys.exit(f"Error: No read permission for '{gff_path}'.")
    except OSError as e:
        sys.exit(f"Error reading GFF file: {e}")

    print(f"Finished loading. Found {len(indexed_features)} contigs.")

    return indexed_features, header_lines


def process_contigs(
    all_features: dict,
    header_lines: list,
    output_dir: str
) -> None:
    """write GFF features for each SeqID.

    Write GFF files per seqID identify in the fasta file and update the ID
    and Parent ID to include the seqID. It can returns GFF with no features
    if no gene model is identify for a specific SeqID.

    Args:
        all_features: Indexed GFF features as returned by load_gff_features().
        header_lines: Header lines to include in the output file.
        output_dir: Path to the output directory to write the GFF files.

    Raises:
        PermissionError: If the output file cannot be written.
        OSError: If any other I/O error occurs during writing.
    """

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        
    try:
        contig_list = list(all_features.keys())
        
        for contig in contig_list:
            contig_features = all_features[contig]
            output_path = os.path.join(output_dir, f"{contig}.gff")
            with open(output_path, 'w') as outfile:
                outfile.writelines(header_lines)
            
                for f in contig_features:
                    f[8] = str(f[8]).replace("ID=","ID="+ contig + ".")
                    f[8] = str(f[8]).replace("Parent=","Parent="+ contig + ".")
                    outfile.write('\t'.join(f) + '\n')

            print(f"-> '{contig}': Wrote {len(contig_features)} features.")
        
    except PermissionError:
        print(f"Error: No write permission for '{output_path}'.", file=sys.stderr)
    except OSError as e:
        print(f"Error writing '{contig}': {e}", file=sys.stderr)


def main() -> None:
    """Parse command-line arguments and launch the GFF divider pipeline."""

    parser = argparse.ArgumentParser(
        description="High-performance GFF spliter per contig."
    )
    parser.add_argument('-g', '--gff', required=True, help='Master GFF input file')
    parser.add_argument('-f', '--fasta', required=True, help='Master fasta input file')
    parser.add_argument('-o', '--output-dir', required=True, help='Directory for output GFFs')
    
    args = parser.parse_args()
    all_features, headers = load_gff_features(args.gff,args.fasta)
    process_contigs(all_features, headers, args.output_dir)
    

if __name__ == '__main__':
    main()
