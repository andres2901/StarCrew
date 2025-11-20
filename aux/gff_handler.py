import argparse
import sys
import os

def slice_gff(input_path, output_path, section_name, target_seqid, slice_start, slice_end, target_strand):
    """
    Reads a GFF file, filters features to a specific slice region, 
    and recalculates coordinates, handling both forward and reverse complement strands.
    """
    
    # 1. Input validation and setup
    if slice_start >= slice_end:
        print(f"Error for section '{section_name}': Start position ({slice_start}) must be strictly less than the end position ({slice_end}). Skipping.", file=sys.stderr)
        return

    # Calculate the length of the new sliced sequence
    SLICE_LENGTH = slice_end - slice_start + 1
    
    print(f"\nProcessing section: {section_name} ({target_seqid}:{slice_start}-{slice_end} [{target_strand}])")
    print(f"Writing to: {output_path}")

    try:
        # Open input file and output file (in the specified output directory)
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            features_processed = 0
            original_header_lines = []
            header_processed = False
            
            # Write a new header indicating the slicing operation
            outfile.write(f"## GFF Slicer Output\n")
            outfile.write(f"## Section: {section_name}\n")
            outfile.write(f"## Original Region: {target_seqid}:{slice_start}-{slice_end} ({target_strand})\n")
            
            # --- Main File Processing Loop ---
            for line in infile:
                # 2. Handle comments and empty lines
                if line.startswith('#'):
                    # Save header/comment lines from the original file
                    if not header_processed:
                        original_header_lines.append(line)
                    continue
                
                # Write original headers if we find the first feature line
                if not header_processed:
                    # Write preserved original header lines only once
                    outfile.writelines(original_header_lines)
                    header_processed = True
                
                if not line.strip():
                    continue

                # Parse the GFF line (TSV format)
                try:
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue # Skip malformed lines

                    f_seqid, f_source, f_type, f_start_str, f_end_str, f_score, f_strand, f_phase, f_attributes = parts
                    
                    # Ensure coordinates are integers
                    f_start = int(f_start_str)
                    f_end = int(f_end_str)

                except ValueError:
                    # Skip lines where start/end aren't valid numbers
                    print(f"Warning: Skipping line with non-numeric coordinates: {line.strip()}", file=sys.stderr)
                    continue
                
                # 3. Filtering: Check if the feature belongs to the target contig 
                # and overlaps the slice region.
                if f_seqid != target_seqid:
                    continue
                    
                # The feature must end after the slice starts AND start before the slice ends.
                if f_end < slice_start or f_start > slice_end:
                    continue

                # --- Overlap Found. Begin coordinate transformation ---
                
                features_processed += 1
                new_start, new_end, new_strand = None, None, None
                
                # Case A: Forward Slice (target_strand is '+')
                if target_strand == '+':
                    # New sequence starts at slice_start. Simple subtraction.
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
                        new_strand = f_strand # Preserve unknown strands

                    # 2. Calculate new coordinates (which are inherently flipped)
                    # N_start is derived from F_end (new coordinate starts at slice_end)
                    n_coord_start = slice_end - f_end + 1
                    # N_end is derived from F_start
                    n_coord_end = slice_end - f_start + 1

                    # GFF format requires start <= end, so we normalize the flipped coordinates.
                    new_start = min(n_coord_start, n_coord_end)
                    new_end = max(n_coord_start, n_coord_end)

                # 4. Clipping: Adjust coordinates for partial overlaps
                # Features that start before the slice (new_start < 1) are reset to 1.
                new_start = max(1, new_start)
                # Features that end after the slice (new_end > SLICE_LENGTH) are clipped.
                new_end = min(SLICE_LENGTH, new_end)
                
                # --- Output new line ---
                
                parts[3] = str(new_start) # Start
                parts[4] = str(new_end)   # End
                parts[6] = new_strand     # Strand
                
                # IMPORTANT: Replace the SeqID with the custom section_name for the new coordinate system
                parts[0] = section_name
                
                outfile.write('\t'.join(parts) + '\n')

            print(f"-> Wrote {features_processed} features for section '{section_name}'.")

    except Exception as e:
        print(f"An unexpected error occurred while processing section '{section_name}': {e}", file=sys.stderr)


def process_coordinate_file(coord_file_path, input_gff_path, output_dir):
    """
    Reads the coordinate file and initiates slicing for each defined section.
    """
    # 1. Ensure the output directory exists
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
        
    # Check if the input GFF file exists
    if not os.path.exists(input_gff_path):
        print(f"Error: Input GFF file not found at '{input_gff_path}'. Aborting.", file=sys.stderr)
        return

    try:
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
                
                # Perform the slicing operation for this specific section
                slice_gff(input_gff_path, output_path, section_name, contig_id, start, end, strand)

    except FileNotFoundError:
        print(f"Error: Coordinate file not found at '{coord_file_path}'.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during file processing: {e}", file=sys.stderr)


def main():
    """Sets up command-line arguments and calls the main processing function."""
    parser = argparse.ArgumentParser(
        description="Batch slices GFF features based on a coordinate file, handling reverse complements."
    )
    
    # Updated command line arguments to use flags
    parser.add_argument('-c', '--coord-file', required=True, dest='coordinates_file',
                        help='Path to the TSV/TXT file defining the sections to slice (SectionName, ContigID, Start, End, Strand).')
    parser.add_argument('-i', '--input-gff', required=True, dest='input_gff',
                        help='Path to the master input GFF/GTF annotation file.')
    parser.add_argument('-o', '--output-dir', required=True, dest='output_folder',
                        help='Path to the directory where all sliced GFF files will be saved (one file per section).')

    # Slice parameters for the coordinate file
    parser.add_argument('-t', '--target-strand', choices=['+', '-'], default='+',
                        help="Target strand for single slice operations. NOTE: This is IGNORED when using a coordinate file, as the strand is read from the file.")

    args = parser.parse_args()
    
    # The dest= argument ensures these variables are named correctly in the args object
    process_coordinate_file(args.coordinates_file, args.input_gff, args.output_folder)

if __name__ == '__main__':
    main()
