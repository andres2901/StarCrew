import os
import glob
import pandas as pd
import argparse

# --- Variable that can be easily changed ---
EXON_RANGE = (4, 10)

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
                # Assuming the third column is the length of the target
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

def get_gff_data(gff_path):
    """Parses a GFF file and returns a DataFrame of gene features and their exon counts."""
    if not os.path.exists(gff_path):
        return pd.DataFrame(), "GFF file not found."
    
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

                    if feature_type == 'gene' and id_part:
                        gff_data.append({
                            'ID': id_part,
                            'Start_Pos': int(gff_parts[3]),
                            'Strand': gff_parts[6],
                            'Exon_Count': 0 # Initialize to 0 for now
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

def check_filters(gene_id, gff_df, exon_range):
    """Checks if a gene meets the exon count and location criteria."""
    gene_info = gff_df[gff_df['ID'] == gene_id].iloc[0]
    total_genes = len(gff_df)
    
    # Check Exon Count
    if not (exon_range[0] <= gene_info['Exon_Count'] <= exon_range[1]):
        return False, f"ID '{gene_id}' has {gene_info['Exon_Count']} exons, not in the required range {exon_range}."

    # Check Location
    position_percent = (gene_info['Relative_Pos'] / total_genes) * 100
    is_start = (gene_info['Strand'] == '+' and position_percent <= 20)
    is_end = (gene_info['Strand'] == '-' and position_percent >= 80)

    if not (is_start or is_end):
        return False, f"ID '{gene_id}' is not at the beginning (for + strand) or end (for - strand)."
    
    return True, "All filters passed."

def run_gff_analysis(hmm_ids, gff_df):
    """Performs the GFF analysis on a given set of IDs."""
    if gff_df.empty:
        return None, "GFF data is empty."
    
    common_ids = set(hmm_ids).intersection(set(gff_df['ID'].values))
    if not common_ids:
        return None, "No common IDs found between HMM and GFF files."

    candidates = gff_df[gff_df['ID'].isin(common_ids)].copy()

    # Apply location and strand filters
    candidates['is_beginning'] = (candidates['Strand'] == '+')
    candidates['is_end'] = (candidates['Strand'] == '-')
    
    beginning_candidates = candidates[candidates['is_beginning']].copy()
    end_candidates = candidates[candidates['is_end']].copy()
    
    # Apply relative position filter for the first 10% of genes on the positive strand
    total_genes = len(gff_df)
    beginning_candidates['percent_pos'] = (beginning_candidates['Relative_Pos'] / total_genes) * 100
    beginning_candidates = beginning_candidates[beginning_candidates['percent_pos'] <= 10]
    
    # Apply relative position filter for the last 10% of genes on the negative strand
    end_candidates['percent_pos'] = (end_candidates['Relative_Pos'] / total_genes) * 100
    end_candidates = end_candidates[end_candidates['percent_pos'] >= 90]

    # Combine the candidates for tie-breaking
    all_candidates = pd.concat([beginning_candidates, end_candidates])
    
    if all_candidates.empty:
        return None, "No genes passed location and strand filters."

    # Sort and select the best one
    all_candidates.sort_values(by=['is_beginning', 'Relative_Pos', 'Exon_Count'], ascending=[False, True, True], inplace=True)
    
    best_gene = all_candidates.iloc[0]['ID']

    return best_gene, "Success"

def process_hmm_files(hmm_folder1, hmm_folder2, hmm_folder3, gff_folder, output_file, empty_output_file, min_common, min_length):
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
                
                unique_candidate = None
                unique_candidate_level = 0
                gff_candidates = set()
                
                common_ids_123 = ids1.intersection(ids2).intersection(ids3)
                common_ids_12 = ids1.intersection(ids2)
                
                if len(common_ids_123) == 1:
                    unique_candidate = list(common_ids_123)[0]
                    unique_candidate_level = 3
                elif len(common_ids_12) == 1:
                    unique_candidate = list(common_ids_12)[0]
                    unique_candidate_level = 2
                elif len(ids1) == 1:
                    unique_candidate = list(ids1)[0]
                    unique_candidate_level = 1

                # Set the list to be used for GFF tie-breaker, if needed
                if len(common_ids_123) > 1:
                    gff_candidates = common_ids_123
                    gff_candidate_level = 3
                elif len(common_ids_12) > 1:
                    gff_candidates = common_ids_12
                    gff_candidate_level = 2
                elif len(ids1) > 1:
                    gff_candidates = ids1
                    gff_candidate_level = 1

                # Step 2: Make the final decision based on the findings
                if unique_candidate is not None and unique_candidate_level >= min_common:
                    final_selected_id = unique_candidate
                    final_reason = f"Finalized with unique ID at level {unique_candidate_level}."
                elif len(gff_candidates) > 1 and gff_candidate_level >= min_common:
                    gff_df, gff_reason = get_gff_data(os.path.join(gff_folder, f'{base_name}.gff'))
                    if gff_df.empty:
                        final_reason = gff_reason
                    else:
                        result, gff_reason = run_gff_analysis(candidates_for_analysis, gff_df)
                        if result:
                            final_selected_id = result
                            final_reason = f"Selected via GFF tie-breaker. Reason: {gff_reason}"
                        else:
                            final_reason = gff_reason
                else:
                    final_reason = "No unique ID found that meets the criteria"

            
            # --- FINAL VALIDATION STEP ---
            if final_selected_id:
                gff_df, _ = get_gff_data(os.path.join(gff_folder, f'{base_name}.gff'))
                if not gff_df.empty:
                    is_valid, validation_reason = check_filters(final_selected_id, gff_df, EXON_RANGE)
                    if is_valid:
                        hmm_scores = df1_agg[df1_agg['ID'] == final_selected_id].iloc[0]
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
                    print(f"No ID could be selected for {base_name}. Reason: GFF file for final validation is empty. Writing to empty file.")
                    empty_f.write(f"{base_name}\n")
            else:
                print(f"No ID could be selected for {base_name}. Reason: {final_reason}. Writing to empty file.")
                empty_f.write(f"{base_name}\n")
    
    if all_results:
        final_df = pd.DataFrame(all_results)
        final_df.to_csv(output_file, sep='\t', index=False, header=False)
        print(f"\nSuccessfully wrote a combined table with {len(final_df)} rows to {output_file}")
    else:
        print("\nNo results were collected to write to the output file.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process HMMSEARCH domtblout and GFF files and save to a single structured table.')
    parser.add_argument('--hmm1', dest='hmm_folder1', required=True, help='Path to the first HMMER folder.')
    parser.add_argument('--hmm2', dest='hmm_folder2', required=True, help='Path to the second HMMER folder.')
    parser.add_argument('--hmm3', dest='hmm_folder3', required=True, help='Path to the third HMMER folder.')
    parser.add_argument('--gff', dest='gff_folder', required=True, help='Path to the folder containing .gff files.')
    parser.add_argument('--output', dest='output_file', required=True, help='Path and filename for the single output TSV file.')
    parser.add_argument('--empty', dest='empty_output_file', required=True, help='Path and filename for the file to list empty results.')
    parser.add_argument('--min_common', type=int, default=1, choices=[1, 2, 3], help='Minimum number of folders with common results to finalize a selection (1, 2, or 3).')
    parser.add_argument('--min_length', type=int, default=500, help='Minimum value for the third column (Length) of the first HMMER file (default: 600).')

    args = parser.parse_args()
    
    process_hmm_files(args.hmm_folder1, args.hmm_folder2, args.hmm_folder3, args.gff_folder, args.output_file, args.empty_output_file, args.min_common, args.min_length)