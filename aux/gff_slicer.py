#!/usr/bin/env python3
"""
GFF Slicer with Strict Containment.

This script indexes a master GFF into memory and slices it based on a 
coordinate file. It ensures features are fully contained within the slice 
and correctly re-sorts features from 1 to N after coordinate transformation.
"""

import argparse
import sys
import os


def load_gff_features(input_path: str) -> tuple[dict, list]:
    """Read and index a GFF file into memory by contig ID.

    Performs a single pass over the file, separating header lines from
    feature lines and grouping features by their SeqID (column 1).

    Args:
        input_path: Path to the input GFF file.

    Returns:
        Tuple (indexed_features, header_lines) where:
            - indexed_features: Dictionary mapping contig ID to a list
                                of feature rows (each row is a list of
                                9 string columns).
            - header_lines: List of comment/header lines preserved
                            for output.

    Raises:
        FileNotFoundError: If the GFF file does not exist.
        PermissionError: If the file cannot be read.
        OSError: If any other I/O error occurs.
    """
    
    indexed_features = {}
    header_lines = []
    
    try:
        with open(input_path, 'r') as infile:
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
        sys.exit(f"Error: GFF file '{input_path}' not found.")
    except PermissionError:
        sys.exit(f"Error: No read permission for '{input_path}'.")
    except OSError as e:
        sys.exit(f"Error reading GFF file: {e}")

    print(f"Finished loading. Found {len(indexed_features)} contigs.")

    return indexed_features, header_lines


def process_section(
    all_features: dict,
    header_lines: list,
    output_path: str,
    section_name: str,
    target_seqid: str,
    slice_start: int,
    slice_end: int,
    target_strand: str
) -> None:
    """Slice, transform, and write GFF features for a single genomic region.

    Filters features using strict containment within the slice boundaries,
    applies coordinate transformation relative to the slice origin, flips
    strand orientation for reverse-complement regions, re-sorts features
    by start position, and writes the result to a new GFF file.

    Args:
        all_features: Indexed GFF features as returned by load_gff_features().
        header_lines: Header lines to include in the output file.
        output_path: Path to the output GFF file to write.
        section_name: Name assigned to this slice, used as the SeqID
                      in the output GFF.
        target_seqid: Contig ID to slice from the indexed features.
        slice_start: Start coordinate of the slice (1-based, inclusive).
        slice_end: End coordinate of the slice (1-based, inclusive).
        target_strand: Strand orientation of the slice ('+' or '-').

    Raises:
        PermissionError: If the output file cannot be written.
        OSError: If any other I/O error occurs during writing.
    """

    if slice_start >= slice_end:
        print(f"Error: Start >= End for '{section_name}'. Skipping.", file=sys.stderr)
        return

    if target_seqid not in all_features:
        print(f"Warning: Contig '{target_seqid}' not found. Skipping.", file=sys.stderr)
        return

    try:
        processed_features = []
        contig_features = all_features[target_seqid]
        
        for parts in contig_features:
            current_parts = parts[:] 
            
            try:
                f_start = int(current_parts[3])
                f_end = int(current_parts[4])
                f_strand = current_parts[6]
            except ValueError:
                continue 

            if f_start < slice_start or f_end > slice_end:
                continue

            if target_strand == '+':
                new_start = f_start - slice_start + 1
                new_end = f_end - slice_start + 1
                new_strand = f_strand
            else:
                if f_strand == '+':
                    new_strand = '-'
                elif f_strand == '-':
                    new_strand = '+'
                else:
                    new_strand = f_strand
                
                n_coord_start = slice_end - f_end + 1
                n_coord_end = slice_end - f_start + 1
                new_start = min(n_coord_start, n_coord_end)
                new_end = max(n_coord_start, n_coord_end)

            current_parts[0] = section_name
            current_parts[3] = new_start
            current_parts[4] = new_end
            current_parts[6] = new_strand
            
            processed_features.append(current_parts)

        processed_features.sort(key=lambda x: x[3])

        with open(output_path, 'w') as outfile:
            outfile.write("## GFF Slicer Output\n")
            outfile.write(f"## Section: {section_name}\n")
            outfile.write(f"## Original Region: {target_seqid}:{slice_start}-{slice_end} ({target_strand})\n")
            outfile.writelines(header_lines)
            
            for f in processed_features:
                f[3], f[4] = str(f[3]), str(f[4])
                outfile.write('\t'.join(f) + '\n')

        print(f"-> '{section_name}': Wrote {len(processed_features)} features.")
        
    except PermissionError:
        print(f"Error: No write permission for '{output_path}'.", file=sys.stderr)
    except OSError as e:
        print(f"Error writing '{section_name}': {e}", file=sys.stderr)


def process_coordinate_file(
    coord_file_path: str,
    input_gff_path: str,
    output_dir: str
) -> None:
    """Orchestrate batch GFF slicing from a coordinate file.

    Loads the master GFF once into memory and iterates over each entry
    in the coordinate file, calling process_section() for each valid row.

    Args:
        coord_file_path: Path to the TSV coordinate file. Each row must
                         contain: Name, Contig, Start, End, Strand.
        input_gff_path: Path to the master GFF file to slice from.
        output_dir: Path to the directory where output GFF files will
                    be written. Created if it does not exist.

    Raises:
        FileNotFoundError: If the coordinate file does not exist.
        PermissionError: If the coordinate file cannot be read.
        OSError: If any other I/O error occurs.
    """

    all_features, header_lines = load_gff_features(input_gff_path)
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    try:
        with open(coord_file_path, 'r') as cf:
            for line in cf:
                line = line.strip()
                if not line or line.startswith('#'): 
                    continue

                fields = line.split()
                if len(fields) != 5: 
                    continue

                name, contig, start_s, end_s, strand = fields
                out_file = os.path.join(output_dir, f"{name}.gff")
                
                process_section(all_features, header_lines, out_file, 
                                name, contig, int(start_s), int(end_s), strand)
        
        print("\nBatch processing complete.")

    except FileNotFoundError:
        sys.exit(f"Error: Coordinate file '{coord_file_path}' not found.")
    except PermissionError:
        sys.exit(f"Error: No read permission for '{coord_file_path}'.")
    except OSError as e:
        sys.exit(f"Error reading coordinate file: {e}")


def main() -> None:
    """Parse command-line arguments and launch the GFF slicing pipeline."""

    parser = argparse.ArgumentParser(
        description="High-performance GFF slicer with strict containment and sorting."
    )
    parser.add_argument('-c', '--coord-file', required=True, help='TSV: Name, Contig, Start, End, Strand')
    parser.add_argument('-i', '--input-gff', required=True, help='Master GFF input file')
    parser.add_argument('-o', '--output-dir', required=True, help='Directory for output GFFs')
    
    args = parser.parse_args()
    process_coordinate_file(args.coord_file, args.input_gff, args.output_dir)


if __name__ == '__main__':
    main()