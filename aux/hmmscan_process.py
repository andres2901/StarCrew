#!/usr/bin/env python3
"""
HMM and GFF Integration Pipeline.

Identifies candidate genes by finding consensus between three HMMER runs 
and validating them against constraints  form a GFF file (exon count and 
genomic position relative to element boundaries).
"""

import os
import glob
import pandas as pd
import argparse

# --- Variables that can be easily changed ---
EXON_RANGE = (2, 11)

def parse_hmm_file(file_path):
    """Parses a single HMMER domtblout file and returns unique IDs"""
    parsed_data = []
    if not os.path.exists(file_path):
        return set(), pd.DataFrame()

    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 19:
                hmm_id = parts[3].strip()
                target_length = parts[5].strip()
                parsed_data.append([hmm_id, target_length])
    
    if not parsed_data:
        return set(), pd.DataFrame()

    df_raw = pd.DataFrame(parsed_data, columns=['ID', 'Length'])
    df_raw['Length'] = pd.to_numeric(df_raw['Length'], errors='coerce')
    
    return set(df_raw['ID'].values), df_raw

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
    """Finds a FASTA file in a folder with matching base name."""
    for ext in extensions:
        file_path = os.path.join(folder, f'{base_name}{ext}')
        if os.path.exists(file_path):
            return file_path
    return None

def get_gff_data(gff_path, fasta_path):
    """Parses a GFF file and returns a DataFrame of gene features and exon counts."""
    if not os.path.exists(gff_path):
        return pd.DataFrame(), "GFF file not found."
    
    fasta_length = get_fasta_length(fasta_path)
    if not fasta_length:
        return pd.DataFrame(), f"FASTA file not found or is empty: {fasta_path}"

    gff_data = []
    mrna_exons = {}
    mrna_to_gene_map = {}
    
    with open(gff_path, 'r') as gff_f:
        for gff_line in gff_f:
            if gff_line.startswith('#') or not gff_line.strip():
                continue
            gff_parts = gff_line.split('\t')
            if len(gff_parts) >= 9:
                try:
                    attributes = gff_parts[8]
                    feature_type = gff_parts[2]
                    id_part, parent_part = None, None
                    
                    for item in attributes.split(';'):
                        if item.startswith('ID='):
                            id_part = item.split('ID=')[1].strip()
                        elif item.startswith('Parent='):
                            parent_part = item.split('Parent=')[1].strip()
                    
                    if feature_type == 'gene' and id_part:
                        gff_data.append({
                            'ID': id_part,
                            'Start_Pos': int(gff_parts[3]),
                            'End_Pos': int(gff_parts[4]),
                            'Strand': gff_parts[6],
                            'Exon_Count': 0,
                            'Element_Start': 1,
                            'Element_End': fasta_length
                        })
                    elif feature_type == 'mRNA' and id_part and parent_part:
                        mrna_to_gene_map[id_part] = parent_part
                    elif feature_type == 'exon' and parent_part:
                        mrna_exons[parent_part] = mrna_exons.get(parent_part, 0) + 1
                except (ValueError, IndexError):
                    continue
    
    if not gff_data:
        return pd.DataFrame(), "No 'gene' features found in GFF."
    
    gff_df = pd.DataFrame(gff_data)
    final_gene_exons = {}
    for mrna_id, exon_count in mrna_exons.items():
        gene_id = mrna_to_gene_map.get(mrna_id)
        if gene_id:
            final_gene_exons[gene_id] = final_gene_exons.get(gene_id, 0) + exon_count

    gff_df['Exon_Count'] = gff_df['ID'].map(final_gene_exons).fillna(0).astype(int)
    gff_df.sort_values(by='Start_Pos', inplace=True)
    gff_df['Relative_Pos'] = range(1, len(gff_df) + 1)
    
    return gff_df, "Success"

def check_filters(gene_id, gff_df, exon_range, pos_range_kb):
    """Checks if a gene meets the exon count and location criteria."""
    matches = gff_df[gff_df['ID'] == gene_id]
    if matches.empty:
        return False, f"ID '{gene_id}' not found in GFF."
        
    gene_info = matches.iloc[0]
    pos_range = pos_range_kb * 1000 
    
    #  Dynamic range adjustment for consistency
    seq_length = gene_info['Element_End'] - gene_info['Element_Start']
    if seq_length < (2 * pos_range):
        pos_range = seq_length / 2

    if not (exon_range[0] <= gene_info['Exon_Count'] <= exon_range[1]):
        return False, f"ID '{gene_id}' has {gene_info['Exon_Count']} exons, not in the required range {exon_range}."

    is_start = (gene_info['Strand'] == '+' and (gene_info['Start_Pos'] - gene_info['Element_Start']) <= pos_range)
    is_end = (gene_info['Strand'] == '-' and (gene_info['Element_End'] - gene_info['End_Pos']) <= pos_range)

    if not (is_start or is_end):
        return False, f"ID '{gene_id}' is not at the beginning (for + strand) or end (for - strand)."
    
    return True, "All filters passed."

def run_gff_analysis(hmm_ids, gff_df, pos_range_kb):
    """Performs tie-breaking based on location and strand."""
    common_ids = set(hmm_ids).intersection(set(gff_df['ID'].values))
    if not common_ids:
        return None, "No common IDs found between HMM and GFF files."

    candidates = gff_df[gff_df['ID'].isin(common_ids)].copy()
    pos_range = pos_range_kb * 1000

    # Dynamic range adjustment based on sequence length
    # We assume all candidates share the same element bounds
    element_start = candidates['Element_Start'].iloc[0]
    element_end = candidates['Element_End'].iloc[0]
    seq_length = element_end - element_start
    
    if seq_length < (2 * pos_range):
        pos_range = seq_length / 2

    cond_plus = (candidates['Strand'] == '+') & (candidates['Start_Pos'] - candidates['Element_Start'] <= pos_range)
    cond_minus = (candidates['Strand'] == '-') & (candidates['Element_End'] - candidates['End_Pos'] <= pos_range)
    
    all_candidates = candidates[cond_plus | cond_minus].copy()
    if all_candidates.empty:
        return None, "No genes passed location and strand filters."

    all_candidates['is_plus'] = all_candidates['Strand'] == '+'
    all_candidates.sort_values(by=['is_plus', 'Relative_Pos', 'Exon_Count'], ascending=[False, True, True], inplace=True)
    
    return all_candidates.iloc[0]['ID'], "Success"

def process_hmm_files(hmm_folder1, hmm_folder2, hmm_folder3, gff_folder, fasta_folder, output_file, empty_output_file, min_common, min_length, range_kb):
    hmm_files1 = glob.glob(os.path.join(hmm_folder1, '*.txt'))
    all_results = []
    
    with open(empty_output_file, 'w') as empty_f:
        for hmm_path1 in hmm_files1:
            base_name = os.path.basename(hmm_path1).split('.')[0]
            print(f"Processing {base_name}...")

            final_selected_id = None
            final_reason = ""

            fasta_path = find_fasta_path(fasta_folder, base_name)
            if not fasta_path:
                final_reason = f"No FASTA file found for {base_name}."
                print(f"No ID selected for {base_name}. Reason: {final_reason}")
                empty_f.write(f"{base_name}\n")
                continue

            gff_df, gff_reason = get_gff_data(os.path.join(gff_folder, f'{base_name}.gff'), fasta_path)
            if gff_df.empty:
                print(f"No ID selected for {base_name}. Reason: {gff_reason}")
                empty_f.write(f"{base_name}\n")
                continue
            
            ids1_full, df1_raw = parse_hmm_file(hmm_path1)
            ids1 = set(df1_raw[df1_raw['Length'] >= min_length]['ID'].values) if not df1_raw.empty else set()
            ids2, _ = parse_hmm_file(os.path.join(hmm_folder2, f'{base_name}.txt'))
            ids3, _ = parse_hmm_file(os.path.join(hmm_folder3, f'{base_name}.txt'))

            common_ids_123 = ids1.intersection(ids2).intersection(ids3)
            common_ids_12 = ids1.intersection(ids2)

            if len(common_ids_123) == 1:
                final_selected_id = list(common_ids_123)[0]
                final_reason = "Finalized with unique ID at level 3."
            elif len(common_ids_123) > 1:
                final_selected_id, g_res = run_gff_analysis(common_ids_123, gff_df, range_kb)
                if final_selected_id:
                    final_reason = f"Selected via GFF tie-breaker at level 3. ({g_res})"
                else:
                    final_reason = f"{g_res}"
            elif len(common_ids_12) == 1 and min_common <= 2:
                final_selected_id = list(common_ids_12)[0]
                final_reason = "Finalized with unique ID at level 2."
            elif len(common_ids_12) > 1 and min_common <= 2:
                final_selected_id, g_res = run_gff_analysis(common_ids_12, gff_df, range_kb)
                if final_selected_id:
                    final_reason = f"Selected via GFF tie-breaker at level 2. ({g_res})"
                else:
                    final_reason = f"{g_res}"
            elif len(ids1) == 1 and min_common == 1:
                final_selected_id = list(ids1)[0]
                final_reason = "Finalized with unique ID at level 1."
            elif len(ids1) > 1 and min_common == 1:
                final_selected_id, g_res = run_gff_analysis(ids1, gff_df, range_kb)
                if final_selected_id:
                    final_reason = f"Selected via GFF tie-breaker at level 1. ({g_res})"
                else:
                    final_reason = f"{g_res}"
            else:
                final_reason = f"No common IDs found at level {min_common}."

            if final_selected_id:
                is_valid, validation_reason = check_filters(final_selected_id, gff_df, EXON_RANGE, range_kb)
                if is_valid:
                    all_results.append({'ID': final_selected_id})
                    print(f"  Success: Selected ID '{final_selected_id}'. Reason: {final_reason}")
                else:
                    print(f"  No ID selected for {base_name}. Reason: {validation_reason}")
                    empty_f.write(f"{base_name}\n")
            else:
                print(f"  No ID selected for {base_name}. Reason: {final_reason}")
                empty_f.write(f"{base_name}\n")
    
    if all_results:
        pd.DataFrame(all_results)['ID'].to_csv(output_file, index=False, header=False)
        print(f"\nSuccessfully wrote {len(all_results)} IDs to {output_file}")

def main():
    """Main CLI execution."""
    parser = argparse.ArgumentParser(description='Process HMMSEARCH domtblout and GFF files.')
    parser.add_argument('--hmm1', dest='hmm_folder1', required=True, help='Path to first HMMER folder.')
    parser.add_argument('--hmm2', dest='hmm_folder2', required=True, help='Path to second HMMER folder.')
    parser.add_argument('--hmm3', dest='hmm_folder3', required=True, help='Path to third HMMER folder.')
    parser.add_argument('--gff', dest='gff_folder', required=True, help='Path to GFF folder.')
    parser.add_argument('--fasta', dest='fasta_folder', required=True, help='Path to FASTA folder.')
    parser.add_argument('--output', dest='output_file', required=True, help='Output filename.')
    parser.add_argument('--empty', dest='empty_output_file', required=True, help='Log for empty results.')
    parser.add_argument('--min_common', type=int, default=2, choices=[1, 2, 3])
    parser.add_argument('--min_length', type=int, default=250)
    parser.add_argument('--range_kb', type=int, default=20)

    args = parser.parse_args()
    process_hmm_files(args.hmm_folder1, args.hmm_folder2, args.hmm_folder3, 
                      args.gff_folder, args.fasta_folder, args.output_file, 
                      args.empty_output_file, args.min_common, args.min_length, args.range_kb)

if __name__ == "__main__":
    main()