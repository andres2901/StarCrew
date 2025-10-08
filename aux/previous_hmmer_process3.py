import os
import glob
import pandas as pd
import argparse

# --- Variables that can be easily changed ---
EXON_RANGE = (3, 10)
POS_RANGE_KB = 20 

def parse_hmm_file(file_path):
    """Parses a single HMMER domtblout file and returns unique IDs and a dataframe of all scores."""
    parsed_data = []
    if not os.path.exists(file_path):
        return set(), pd.DataFrame(), pd.DataFrame()

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 19:
                hmm_id = parts[0].strip()
                target_length = parts[2].strip()
                score2 = parts[17]
                score3 = parts[18]
                parsed_data.append([hmm_id, target_length, score2, score3])
    
    if not parsed_data:
        return set(), pd.DataFrame(), pd.DataFrame()

    df_raw = pd.DataFrame(parsed_data, columns=['ID', 'Length', 'Score_2', 'Score_3'])
    df_raw['Length'] = pd.to_numeric(df_raw['Length'])
    df_raw['Score_2'] = pd.to_numeric(df_raw['Score_2'])
    df_raw['Score_3'] = pd.to_numeric(df_raw['Score_3'])

    df_agg = df_raw.groupby('ID').agg(
        min_score2=('Score_2', 'min'),
        max_score3=('Score_3', 'max')
    ).reset_index()
    
    return set(df_raw['ID'].values), df_raw, df_agg

def get_fasta_length(fasta_path):
    """Parses a single FASTA file and returns the total sequence length."""
    length = 0
    if not os.path.exists(fasta_path):
        return None
    with open(fasta_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                length += len(line.strip())
    return length

def find_fasta_path(folder, base_name, extensions=['.fa', '.fasta', '.fna']):
    """Finds a FASTA file in a folder with a given base name and a list of possible extensions."""
    for ext in extensions:
        file_path = os.path.join(folder, f'{base_name}{ext}')
        if os.path.exists(file_path):
            return file_path
    return None

def get_gff_data(gff_path, fasta_path):
    """Parses a GFF file and returns a DataFrame of gene features and their exon counts, using FASTA for length."""
    if not os.path.exists(gff_path):
        return pd.DataFrame(), "GFF file not found."
    
    fasta_length = get_fasta_length(fasta_path)
    if not fasta_length:
        return pd.DataFrame(), f"FASTA file not found or is empty: {fasta_path}"

    element_start = 1
    element_end = fasta_length
    print(f"Using FASTA file for element length: {fasta_length}bp")

    gff_data = []
    mrna_exons = {}
    mrna_to_gene_map = {}
    
    with open(gff_path, 'r') as gff_f:
        for gff_line in gff_f:
            if gff_line.startswith('#'):
                continue
            gff_parts = gff_line.split('\t')
            if len(gff_parts) >= 9:
                try:
                    attributes = gff_parts[8]
                    feature_type = gff_parts[2]
                    
                    id_part = None
                    parent_part = None
                    
                    for part in attributes.split(';'):
                        if part.startswith('ID='):
                            id_part = part.split('ID=')[1].strip()
                        elif part.startswith('Parent='):
                            parent_part = part.split('Parent=')[1].strip()
                    
                    start_pos = int(gff_parts[3])
                    end_pos = int(gff_parts[4])

                    if feature_type == 'gene' and id_part:
                        gff_data.append({
                            'ID': id_part,
                            'Start_Pos': start_pos,
                            'End_Pos': end_pos,
                            'Strand': gff_parts[6],
                            'Exon_Count': 0
                        })
                    elif feature_type == 'mRNA' and id_part and parent_part:
                        mrna_to_gene_map[id_part] = parent_part
                    elif feature_type == 'exon' and parent_part:
                        if parent_part in mrna_exons:
                            mrna_exons[parent_part] += 1
                        else:
                            mrna_exons[parent_part] = 1

                except (ValueError, IndexError):
                    continue
    
    if not gff_data:
        return pd.DataFrame(), "No 'gene' features found in GFF file or parsing failed."
    
    gff_df = pd.DataFrame(gff_data)
    gff_df['Element_Start'] = element_start
    gff_df['Element_End'] = element_end
    
    # Map exon counts from mRNA to the corresponding gene
    final_gene_exons = {}
    for mrna_id, exon_count in mrna_exons.items():
        gene_id = mrna_to_gene_map.get(mrna_id)
        if gene_id:
            if gene_id in final_gene_exons:
                final_gene_exons[gene_id] += exon_count
            else:
                final_gene_exons[gene_id] = exon_count

    gff_df['Exon_Count'] = gff_df['ID'].map(final_gene_exons).fillna(0).astype(int)
    gff_df.sort_values(by='Start_Pos', inplace=True)
    gff_df['Relative_Pos'] = range(1, len(gff_df) + 1)
    
    return gff_df, "Success"

def check_filters(gene_id, gff_df, exon_range, pos_range_kb):
    """Checks if a gene meets the exon count and location criteria."""
    gene_info = gff_df[gff_df['ID'] == gene_id].iloc[0]
    element_start = gene_info['Element_Start']
    element_end = gene_info['Element_End']
    pos_range = pos_range_kb * 1000 # Convert to base pairs
    
    # Check Exon Count
    if not (exon_range[0] <= gene_info['Exon_Count'] <= exon_range[1]):
        return False, f"ID '{gene_id}' has {gene_info['Exon_Count']} exons, not in the required range {exon_range}."

    # Check Location (absolute position)
    is_start = (gene_info['Strand'] == '+' and (gene_info['Start_Pos'] - element_start) <= pos_range)
    is_end = (gene_info['Strand'] == '-' and (element_end - gene_info['End_Pos']) <= pos_range)

    if not (is_start or is_end):
        return False, f"ID '{gene_id}' is not at the beginning (for + strand) or end (for - strand)."
    
    return True, "All filters passed."

def run_gff_analysis(hmm_ids, gff_df, pos_range_kb):
    """Performs the GFF analysis on a given set of IDs."""
    if gff_df.empty:
        return None, "GFF data is empty."
    
    common_ids = set(hmm_ids).intersection(set(gff_df['ID'].values))
    if not common_ids:
        return None, "No common IDs found between HMM and GFF files."

    candidates = gff_df[gff_df['ID'].isin(common_ids)].copy()
    
    element_start = candidates['Element_Start'].iloc[0]
    element_end = candidates['Element_End'].iloc[0]
    pos_range = pos_range_kb * 1000 # Convert to base pairs

    # Apply location and strand filters
    candidates['is_beginning'] = (candidates['Strand'] == '+')
    candidates['is_end'] = (candidates['Strand'] == '-')
    
    beginning_candidates = candidates[candidates['is_beginning']].copy()
    end_candidates = candidates[candidates['is_end']].copy()
    
    # Apply absolute position filter for the first 10kb of genes on the positive strand
    beginning_candidates['dist_from_start'] = beginning_candidates['Start_Pos'] - element_start
    beginning_candidates = beginning_candidates[beginning_candidates['dist_from_start'] <= pos_range]
    
    # Apply absolute position filter for the last 10kb of genes on the negative strand
    end_candidates['dist_from_end'] = element_end - end_candidates['End_Pos']
    end_candidates = end_candidates[end_candidates['dist_from_end'] <= pos_range]

    # Combine the candidates for tie-breaking
    all_candidates = pd.concat([beginning_candidates, end_candidates])
    
    if all_candidates.empty:
        return None, "No genes passed location and strand filters."

    # Sort and select the best one
    all_candidates.sort_values(by=['is_beginning', 'Relative_Pos', 'Exon_Count'], ascending=[False, True, True], inplace=True)
    
    best_gene = all_candidates.iloc[0]['ID']

    return best_gene, "Success"

def process_hmm_files(hmm_folder1, hmm_folder2, hmm_folder3, gff_folder, fasta_folder, output_file, empty_output_file, min_common, min_length):
    hmm_files1 = glob.glob(os.path.join(hmm_folder1, '*.txt'))
    
    if not hmm_files1:
        print(f"No .txt files found in {hmm_folder1}")
        return

    all_results = []
    
    with open(empty_output_file, 'w') as empty_f:
        for hmm_path1 in hmm_files1:
            base_name = os.path.basename(hmm_path1).split('.')[0]
            print(f"\nProcessing {base_name}...")

            final_selected_id = None
            final_reason = ""
            
            fasta_path = find_fasta_path(fasta_folder, base_name)
            if not fasta_path:
                final_reason = f"No FASTA file found for {base_name} with .fa, .fasta, or .fna extensions."
                print(f"No ID could be selected for {base_name}. Reason: {final_reason}. Writing to empty file.")
                empty_f.write(f"{base_name}\n")
                continue

            # Get GFF data and element length from FASTA
            gff_df, gff_reason = get_gff_data(os.path.join(gff_folder, f'{base_name}.gff'), fasta_path)
            if gff_df.empty:
                final_reason = gff_reason
                print(f"No ID could be selected for {base_name}. Reason: {final_reason}. Writing to empty file.")
                empty_f.write(f"{base_name}\n")
                continue
            
            ids1_full, df1_raw, df1_agg = parse_hmm_file(hmm_path1)
            
            if df1_raw.empty:
                final_reason = "First HMM file is empty. Skipping."
                print(f"No ID could be selected for {base_name}. Reason: {final_reason}. Writing to empty file.")
                empty_f.write(f"{base_name}\n")
                continue
            
            df1_filtered = df1_raw[df1_raw['Length'] >= min_length]
            ids1 = set(df1_filtered['ID'].values)
            
            if not ids1:
                final_reason = f"No IDs found in HMM folder 1 after filtering with min_length >= {min_length}."
            else:
                hmm_path2 = os.path.join(hmm_folder2, f'{base_name}.txt')
                hmm_path3 = os.path.join(hmm_folder3, f'{base_name}.txt')
                ids2, _, _ = parse_hmm_file(hmm_path2)
                ids3, _, _ = parse_hmm_file(hmm_path3)
                
                candidates_for_analysis = set()
                if min_common == 3:
                    candidates_for_analysis = ids1.intersection(ids2).intersection(ids3)
                elif min_common == 2:
                    candidates_for_analysis = ids1.intersection(ids2)
                else: # min_common == 1
                    candidates_for_analysis = ids1
                
                if len(candidates_for_analysis) == 1:
                    final_selected_id = list(candidates_for_analysis)[0]
                    final_reason = f"Unique ID found at level {min_common}."
                elif len(candidates_for_analysis) > 1:
                    result, gff_reason_analysis = run_gff_analysis(candidates_for_analysis, gff_df, POS_RANGE_KB)
                    if result:
                        final_selected_id = result
                        final_reason = f"Selected via GFF tie-breaker. Reason: {gff_reason_analysis}"
                    else:
                        final_reason = gff_reason_analysis
                else:
                    final_reason = f"No common IDs found at level {min_common}."
            
            # --- FINAL VALIDATION STEP ---
            if final_selected_id:
                is_valid, validation_reason = check_filters(final_selected_id, gff_df, EXON_RANGE, POS_RANGE_KB)
                if is_valid:
                    hmm_scores = df1_agg[df1_agg['ID'] == final_selected_id].iloc[0]
                    # NOTE: Keeping the score data here, but it will be filtered out on save.
                    all_results.append({
                        'ID': final_selected_id,
                        'HMM_Min_Score2': hmm_scores['min_score2'],
                        'HMM_Max_Score3': hmm_scores['max_score3']
                    })
                    print(f"Success: Selected ID '{final_selected_id}'. Reason: {final_reason}")
                else:
                    print(f"No ID could be selected for {base_name}. Reason: {validation_reason}. Writing to empty file.")
                    empty_f.write(f"{base_name}\n")
            else:
                print(f"No ID could be selected for {base_name}. Reason: {final_reason}. Writing to empty file.")
                empty_f.write(f"{base_name}\n")
    
    if all_results:
        final_df = pd.DataFrame(all_results)
        # MODIFICATION: Only write the 'ID' column to the output file.
        final_df['ID'].to_csv(output_file, index=False, header=False)
        print(f"\nSuccessfully wrote a list of {len(final_df)} IDs to {output_file}")
    else:
        print("\nNo results were collected to write to the output file.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process HMMSEARCH domtblout and GFF files and save to a single structured table.')
    parser.add_argument('--hmm1', dest='hmm_folder1', required=True, help='Path to the first HMMER folder.')
    parser.add_argument('--hmm2', dest='hmm_folder2', required=True, help='Path to the second HMMER folder.')
    parser.add_argument('--hmm3', dest='hmm_folder3', required=True, help='Path to the third HMMER folder.')
    parser.add_argument('--gff', dest='gff_folder', required=True, help='Path to the folder containing .gff files.')
    parser.add_argument('--fasta', dest='fasta_folder', required=True, help='Path to the folder containing FASTA files for element lengths.')
    parser.add_argument('--output', dest='output_file', required=True, help='Path and filename for the single output TSV file.')
    parser.add_argument('--empty', dest='empty_output_file', required=True, help='Path and filename for the file to list empty results.')
    parser.add_argument('--min_common', type=int, default=1, choices=[1, 2, 3], help='Minimum number of folders with common results to finalize a selection (1, 2, or 3).')
    parser.add_argument('--min_length', type=int, default=600, help='Minimum value for the third column (Length) of the first HMMER file (default: 600).')

    args = parser.parse_args()
    
    process_hmm_files(args.hmm_folder1, args.hmm_folder2, args.hmm_folder3, args.gff_folder, args.fasta_folder, args.output_file, args.empty_output_file, args.min_common, args.min_length)