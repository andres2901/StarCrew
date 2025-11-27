import argparse
import sys
import os

def load_gff_features(input_path):
    """
    Reads the entire GFF file into memory once.
    Returns: A tuple containing (dict of {contig_id: [list_of_feature_parts]}, list of header lines).
    """
    print(f"STEP 1: Loading and indexing GFF features from '{input_path}'...")
    
    # Dictionary to hold features, indexed by Contig ID (SeqID)
    indexed_features = {}
    header_lines = []
    
    try:
        with open(input_path, 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('#'):
                    # Capture header lines to be written to all output files
                    header_lines.append(line + '\n')
                    continue
                
                # Parse the GFF line (TSV format)
                parts = line.split('\t')
                if len(parts) < 9:
                    print(f"Warning: Skipping malformed GFF line: {line[:50]}...", file=sys.stderr)
                    continue

                f_seqid = parts[0]
                
                # Store the entire list of parts for later processing
                if f_seqid not in indexed_features:
                    indexed_features[f_seqid] = []
                indexed_features[f_seqid].append(parts)

    except Exception as e:
        print(f"Error during GFF loading: {e}", file=sys.stderr)
        return {}, []

    print(f"Finished loading. Found {len(indexed_features)} contigs.")
    return indexed_features, header_lines


def process_section(all_features, header_lines, output_path, section_name, target_seqid, slice_start, slice_end, target_strand):
    """
    Applies filtering and coordinate transformation for a single section using 
    the in-memory feature data.
    """
    
    # 1. Input validation and setup
    if slice_start >= slice_end:
        print(f"Error for section '{section_name}': Start position ({slice_start}) must be strictly less than the end position ({slice_end}). Skipping.", file=sys.stderr)
        return

    # Check if the target contig exists in the loaded data
    if target_seqid not in all_features:
        print(f"Warning: Contig '{target_seqid}' not found in the input GFF file. Skipping section '{section_name}'.", file=sys.stderr)
        return

    # Calculate the length of the new sliced sequence
    SLICE_LENGTH = slice_end - slice_start + 1
    
    print(f"-> Processing section: {section_name} ({target_seqid}:{slice_start}-{slice_end} [{target_strand}])")

    try:
        with open(output_path, 'w') as outfile:
            features_processed = 0
            
            # Write custom header
            outfile.write(f"## GFF Slicer Output\n")
            outfile.write(f"## Section: {section_name}\n")
            outfile.write(f"## Original Region: {target_seqid}:{slice_start}-{slice_end} ({target_strand})\n")
            
            # Write original file headers
            outfile.writelines(header_lines)
            
            # Get only the features for the current contig
            contig_features = all_features[target_seqid]
            
            # --- Feature Processing Loop (In-Memory) ---
            for parts in contig_features:
                
                # Note: f_seqid is already known to be target_seqid since we indexed it.
                f_start_str, f_end_str, f_strand = parts[3], parts[4], parts[6]
                
                try:
                    f_start = int(f_start_str)
                    f_end = int(f_end_str)
                except ValueError:
                    # Should be handled during initial load, but included for safety
                    continue 

                # 2. Filtering: Check if the feature overlaps the slice region.
                # The feature must end after the slice starts AND start before the slice ends.
                if f_start < slice_start or f_end > slice_end:
                    continue

                # --- Overlap Found. Begin coordinate transformation ---
                
                features_processed += 1
                new_start, new_end, new_strand = None, None, None
                
                # Case A: Forward Slice (target_strand is '+')
                if target_strand == '+':
                    new_start = f_start - slice_start + 1
                    new_end = f_end - slice_start + 1
                    new_strand = f_strand
                    
                # Case B: Reverse Complement Slice (target_strand is '-')
                elif target_strand == '-':
                    # 1. Flip the feature's strand
                    if f_strand == '+':
                        new_strand = '-'
                    elif f_strand == '-':
                        new_strand = '+'
                    else:
                        new_strand = f_strand 

                    # 2. Calculate new coordinates (which are inherently flipped)
                    n_coord_start = slice_end - f_end + 1
                    n_coord_end = slice_end - f_start + 1

                    # GFF format requires start <= end
                    new_start = min(n_coord_start, n_coord_end)
                    new_end = max(n_coord_start, n_coord_end)

                # 3. Clipping: Adjust coordinates for partial overlaps
                new_start = max(1, new_start)
                new_end = min(SLICE_LENGTH, new_end)
                
                # --- Output new line ---
                
                parts[0] = section_name      # SeqID (renamed)
                parts[3] = str(new_start)    # Start
                parts[4] = str(new_end)      # End
                parts[6] = new_strand        # Strand
                
                outfile.write('\t'.join(parts) + '\n')

            print(f"-> Wrote {features_processed} features for section '{section_name}'.")

    except Exception as e:
        print(f"An unexpected error occurred while writing section '{section_name}': {e}", file=sys.stderr)


def process_coordinate_file(coord_file_path, input_gff_path, output_dir):
    """
    Reads the coordinate file, loads the master GFF once, and initiates 
    slicing for each defined section.
    """
    
    # 0. Load the master GFF file into memory (the biggest performance gain)
    all_features, header_lines = load_gff_features(input_gff_path)
    if not all_features:
        print("Aborting due to empty or unreadable GFF file data.", file=sys.stderr)
        return
    
    # 1. Ensure the output directory exists
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    print("STEP 2: Processing sections defined in coordinate file...")

    try:
        sections_processed = 0
        with open(coord_file_path, 'r') as cf:
            for i, line in enumerate(cf):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                # Expected format: SectionName ContigID Start End Strand
                fields = line.split()
                if len(fields) != 5:
                    print(f"Warning: Skipping line {i+1} in coordinate file. Expected 5 fields, got {len(fields)}: {line}", file=sys.stderr)
                    continue

                section_name, contig_id, start_str, end_str, strand = fields

                try:
                    start = int(start_str)
                    end = int(end_str)
                except ValueError:
                    print(f"Warning: Skipping section '{section_name}'. Start/End coordinates must be integers: {line}", file=sys.stderr)
                    continue

                if strand not in ['+', '-']:
                    print(f"Warning: Skipping section '{section_name}'. Strand must be '+' or '-': {line}", file=sys.stderr)
                    continue

                # Define the unique output path for this section
                output_filename = f"{section_name}.gff"
                output_path = os.path.join(output_dir, output_filename)
                
                # Process this section using the in-memory data
                process_section(all_features, header_lines, output_path, section_name, contig_id, start, end, strand)
                sections_processed += 1
        
        print(f"\nCompleted processing of {sections_processed} sections.")

    except FileNotFoundError:
        print(f"Error: Coordinate file not found at '{coord_file_path}'.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during coordinate file processing: {e}", file=sys.stderr)


def main():
    """Sets up command-line arguments and calls the main processing function."""
    parser = argparse.ArgumentParser(
        description="Batch slices GFF features based on a coordinate file, handling reverse complements."
    )
    
    parser.add_argument('-c', '--coord-file', required=True, dest='coordinates_file',
                        help='Path to the TSV/TXT file defining the sections to slice (SectionName, ContigID, Start, End, Strand).')
    parser.add_argument('-i', '--input-gff', required=True, dest='input_gff',
                        help='Path to the master input GFF/GTF annotation file.')
    parser.add_argument('-o', '--output-dir', required=True, dest='output_folder',
                        help='Path to the directory where all sliced GFF files will be saved (one file per section).')

    # This argument is only kept for backward compatibility/help clarity but is ignored
    parser.add_argument('-t', '--target-strand', choices=['+', '-'], default='+',
                        help=argparse.SUPPRESS) # Hide in help output since it's ignored

    args = parser.parse_args()
    
    process_coordinate_file(args.coordinates_file, args.input_gff, args.output_folder)

if __name__ == '__main__':
    main()