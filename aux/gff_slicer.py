import argparse
import sys
import os

def load_gff_features(input_path):
    """
    Reads the entire GFF file into memory once to optimize performance for batch processing.
    Returns: A tuple (dict of {contig_id: [list_of_feature_parts]}, list of header lines).
    """
    print(f"STEP 1: Loading and indexing GFF features from '{input_path}'...")
    
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
                    continue

                f_seqid = parts[0]
                if f_seqid not in indexed_features:
                    indexed_features[f_seqid] = []
                
                # Store the original parts list for later processing
                indexed_features[f_seqid].append(parts)

    except Exception as e:
        print(f"Error during GFF loading: {e}", file=sys.stderr)
        return {}, []

    print(f"Finished loading. Found {len(indexed_features)} contigs.")
    return indexed_features, header_lines


def process_section(all_features, header_lines, output_path, section_name, target_seqid, slice_start, slice_end, target_strand):
    """
    Slices features using STRICT CONTAINMENT. 
    Correctly flips strands and ensures output is sorted by start position (1 to N).
    """
    
    if slice_start >= slice_end:
        print(f"Error: Start ({slice_start}) must be < End ({slice_end}) for section '{section_name}'.", file=sys.stderr)
        return

    if target_seqid not in all_features:
        print(f"Warning: Contig '{target_seqid}' not found. Skipping '{section_name}'.", file=sys.stderr)
        return

    try:
        processed_features = []
        contig_features = all_features[target_seqid]
        
        for parts in contig_features:
            # IMPORTANT: Copy the list so modification doesn't affect other slices in memory
            current_parts = parts[:] 
            
            try:
                f_start = int(current_parts[3])
                f_end = int(current_parts[4])
                f_strand = current_parts[6]
            except ValueError:
                continue 

            # --- 1. STRICT CONTAINMENT FILTER ---
            # Feature must be 100% inside the slice boundaries [slice_start, slice_end]
            if f_start < slice_start or f_end > slice_end:
                continue

            # --- 2. COORDINATE TRANSFORMATION & STRAND FLIPPING ---
            if target_strand == '+':
                # Forward orientation: simple relative offset
                new_start = f_start - slice_start + 1
                new_end = f_end - slice_start + 1
                new_strand = f_strand
            else:
                # Reverse orientation ('-'):
                # A: Flip the feature strand correctly (e.g. + becomes -)
                if f_strand == '+':
                    new_strand = '-'
                elif f_strand == '-':
                    new_strand = '+'
                else:
                    new_strand = f_strand
                
                # B: Recalculate coordinates relative to the slice_end
                n_coord_start = slice_end - f_end + 1
                n_coord_end = slice_end - f_start + 1
                
                new_start = min(n_coord_start, n_coord_end)
                new_end = max(n_coord_start, n_coord_end)

            # Update the parts list
            current_parts[0] = section_name
            current_parts[3] = new_start  # Store as int for sorting
            current_parts[4] = new_end    # Store as int for sorting
            current_parts[6] = new_strand
            
            processed_features.append(current_parts)

        # --- 3. RE-SORT BY NEW START POSITION ---
        # Ensures GFF is organized 1 -> N even for negative strand slices
        processed_features.sort(key=lambda x: x[3])

        # --- 4. WRITE TO FILE ---
        with open(output_path, 'w') as outfile:
            outfile.write(f"## GFF Slicer Output\n")
            outfile.write(f"## Section: {section_name}\n")
            outfile.write(f"## Original Region: {target_seqid}:{slice_start}-{slice_end} ({target_strand})\n")
            outfile.writelines(header_lines)
            
            for f in processed_features:
                # Convert coordinates back to strings for joining
                f[3] = str(f[3])
                f[4] = str(f[4])
                outfile.write('\t'.join(f) + '\n')

        print(f"-> Section '{section_name}': Wrote {len(processed_features)} features.")

    except Exception as e:
        print(f"Error processing '{section_name}': {e}", file=sys.stderr)


def process_coordinate_file(coord_file_path, input_gff_path, output_dir):
    """Reads coordinate file and initiates batch processing."""
    all_features, header_lines = load_gff_features(input_gff_path)
    if not all_features: 
        return
    
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
                process_section(all_features, header_lines, os.path.join(output_dir, f"{name}.gff"), 
                                name, contig, int(start_s), int(end_s), strand)
        
        print("\nBatch processing complete.")

    except Exception as e:
        print(f"Error in coordinate file: {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Batch GFF slicer with strict containment, strand flipping, and sorting.")
    parser.add_argument('-c', '--coord-file', required=True, help='TSV: Name, Contig, Start, End, Strand')
    parser.add_argument('-i', '--input-gff', required=True, help='Master GFF input file')
    parser.add_argument('-o', '--output-dir', required=True, help='Directory for output GFFs')
    args = parser.parse_args()
    
    process_coordinate_file(args.coord_file, args.input_gff, args.output_dir)

if __name__ == '__main__':
    main()