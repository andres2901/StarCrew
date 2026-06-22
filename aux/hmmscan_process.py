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
import sys


# CONSTANT
EXON_RANGE = (2, 11)


def parse_hmm_file(file_path: str) -> tuple[set, pd.DataFrame]:
    """Parse a single HMMER domtblout file and extract unique hit IDs.

    Args:
        file_path: Path to the HMMER domtblout output file.

    Returns:
        Tuple (id_set, df_raw) where:
            - id_set: Set of unique hit IDs found in the file.
            - df_raw: DataFrame with columns ['ID', 'Length'] for
                      all parsed hits.
    """

    parsed_data = []
    if not os.path.exists(file_path):
        return set(), pd.DataFrame()

    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 19:
                    hmm_id = parts[3].strip()
                    target_length = parts[5].strip()
                    parsed_data.append([hmm_id, target_length])
    except PermissionError:
        return set(), pd.DataFrame()
    except OSError as e:
        print(f"Error reading HMM file '{file_path}': {e}", file=sys.stderr)
        return set(), pd.DataFrame()
    
    if not parsed_data:
        return set(), pd.DataFrame()

    df_raw = pd.DataFrame(parsed_data, columns=['ID', 'Length'])
    df_raw['Length'] = pd.to_numeric(df_raw['Length'], errors='coerce')
    
    return set(df_raw['ID'].values), df_raw


def get_fasta_length(fasta_path: str) -> int | None:
    """Calculate the total sequence length of a FASTA file.

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        Total number of nucleotide characters across all sequences,
        or None if the file does not exist or cannot be read.
    """

    length = 0
    if not os.path.exists(fasta_path):
        return None
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    length += len(line.strip())
    except PermissionError:
        return None
    except OSError as e:
        print(f"Error reading FASTA file '{fasta_path}': {e}", file=sys.stderr)
        return None
    
    return length


def find_fasta_path(
    folder: str,
    base_name: str,
    extensions: list[str] | None = None
) -> str | None:
    """Search for a FASTA file matching a base name in a folder.

    Args:
        folder: Path to the folder to search in.
        base_name: Base filename without extension.
        extensions: List of extensions to try, in order.
                    Defaults to ['.fa', '.fasta', '.fna'].

    Returns:
        Full path to the first matching file found, or None if
        no match exists.
    """

    if extensions is None:
        extensions = ['.fa', '.fasta', '.fna']

    for ext in extensions:
        file_path = os.path.join(folder, f'{base_name}{ext}')
        if os.path.exists(file_path):
            return file_path
    return None


def get_gff_data(gff_path: str, fasta_path: str) -> tuple[pd.DataFrame, str]:
    """Parse a GFF file and return a DataFrame of gene features with exon counts.

    Performs a single pass over the GFF, collecting gene coordinates,
    mapping mRNA IDs to parent genes, and counting exons per gene.
    Uses the FASTA file to determine the full element length.

    Args:
        gff_path: Path to the input GFF file.
        fasta_path: Path to the corresponding FASTA file, used to
                    determine total sequence length.

    Returns:
        Tuple (gff_df, message) where:
            - gff_df: DataFrame with columns ID, Start_Pos, End_Pos,
                      Strand, Exon_Count, Element_Start, Element_End,
                      and Relative_Pos. Empty if parsing fails.
            - message: Status string, 'Success' or an error description.
    """

    if not os.path.exists(gff_path):
        return pd.DataFrame(), "GFF file not found."
    
    fasta_length = get_fasta_length(fasta_path)
    if not fasta_length:
        return pd.DataFrame(), f"FASTA file not found or is empty: {fasta_path}"

    gff_data = []
    mrna_exons = {}
    mrna_to_gene_map = {}

    try:
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
    except PermissionError:
        return pd.DataFrame(), f"No read permission for '{gff_path}'."
    except OSError as e:
        return pd.DataFrame(), f"Error reading GFF file: {e}"

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


def check_filters(
    gene_id: str,
    gff_df: pd.DataFrame,
    exon_range: tuple[int, int],
    pos_range_kb: int
) -> tuple[bool, str]:
    """Validate a single gene against exon count and location filters.

    Args:
        gene_id: Gene ID to look up in the GFF DataFrame.
        gff_df: DataFrame as returned by get_gff_data().
        exon_range: Tuple (min, max) defining the accepted exon count range.
        pos_range_kb: Maximum allowed distance in kb from the element
                      boundary, depending on strand orientation.

    Returns:
        Tuple (passed, message) where:
            - passed: True if the gene passes all filters.
            - message: Description of the result or the reason for failure.
    """

    matches = gff_df[gff_df['ID'] == gene_id]
    if matches.empty:
        return False, f"ID '{gene_id}' not found in GFF."
        
    gene_info = matches.iloc[0]
    pos_range = pos_range_kb * 1000 
    
    seq_length = gene_info['Element_End'] - gene_info['Element_Start']
    if seq_length < (2 * pos_range):
        pos_range = seq_length / 2

    if not (exon_range[0] <= gene_info['Exon_Count'] <= exon_range[1]):
        return False, (f"ID '{gene_id}' has {gene_info['Exon_Count']} exons, not in the required range {exon_range}.", f"Exon_count")

    is_start = (gene_info['Strand'] == '+' and (gene_info['Start_Pos'] - gene_info['Element_Start']) <= pos_range)
    is_end = (gene_info['Strand'] == '-' and (gene_info['Element_End'] - gene_info['End_Pos']) <= pos_range)

    if not (is_start or is_end):
        return False, (f"ID '{gene_id}' is not at the beginning (for + strand) or end (for - strand).", f"Boundary")
    
    return True, "All filters passed."


def run_gff_analysis(
    hmm_ids: set,
    gff_df: pd.DataFrame,
    pos_range_kb: int
) -> tuple[str | None, str]:
    """Select the best candidate gene from a set of HMM hits using GFF data.

    Filters candidates by location and strand, then ranks them by strand
    preference, relative genomic position, and exon count.

    Args:
        hmm_ids: Set of gene IDs identified by HMMER.
        gff_df: DataFrame as returned by get_gff_data().
        pos_range_kb: Maximum allowed distance in kb from the element
                      boundary for a candidate to be considered.

    Returns:
        Tuple (selected_id, message) where:
            - selected_id: ID of the top-ranked candidate, or None if
                           no candidate passes the filters.
            - message: Description of the result or reason for failure.
    """

    common_ids = set(hmm_ids).intersection(set(gff_df['ID'].values))
    if not common_ids:
        return None, (f"No common IDs found between HMM and GFF files.", f"File")

    candidates = gff_df[gff_df['ID'].isin(common_ids)].copy()
    pos_range = pos_range_kb * 1000

    element_start = candidates['Element_Start'].iloc[0]
    element_end = candidates['Element_End'].iloc[0]
    seq_length = element_end - element_start
    
    if seq_length < (2 * pos_range):
        pos_range = seq_length / 2

    cond_plus = (candidates['Strand'] == '+') & (candidates['Start_Pos'] - candidates['Element_Start'] <= pos_range)
    cond_minus = (candidates['Strand'] == '-') & (candidates['Element_End'] - candidates['End_Pos'] <= pos_range)
    
    all_candidates = candidates[cond_plus | cond_minus].copy()
    if all_candidates.empty:
        return None, (f"No genes passed location and strand filters.", f"Boundary")

    all_candidates['is_plus'] = all_candidates['Strand'] == '+'
    all_candidates.sort_values(by=['is_plus', 'Relative_Pos', 'Exon_Count'], ascending=[False, True, True], inplace=True)
    
    return all_candidates.iloc[0]['ID'], "Success"


def _effective_pos_range(gene_info: pd.Series, pos_range_kb: int) -> float:
    """Compute the position-range threshold (bp), shrunk for short elements.

    Mirrors the logic in check_filters()/run_gff_analysis(): if the element
    is shorter than twice the requested kb range, the range is capped at
    half the element length so the two boundary windows never overlap.
    """
    pos_range = pos_range_kb * 1000
    seq_length = gene_info['Element_End'] - gene_info['Element_Start']
    if seq_length < (2 * pos_range):
        pos_range = seq_length / 2
    return pos_range


def _gene_failure_distance(
    gene_info: pd.Series,
    exon_range: tuple[int, int],
    pos_range_kb: int,
    hit_length: float | None = None,
    min_length: int | None = None
) -> tuple[float, str]:
    """Score how far a single gene is from passing the acceptance filters.
 
    A distance of 0 means the gene would pass the exon-count, position, and
    (when applicable) length filters. Larger values mean further from
    passing. The exon distance is in exon-count units; the position
    distance is converted to kb; the length distance is normalized by /100
    so all three terms sit on roughly comparable scales and can be summed.
 
    Args:
        gene_info: A single row (Series) from gff_df.
        exon_range: Tuple (min, max) accepted exon count.
        pos_range_kb: Distance threshold in kb from element boundary.
        hit_length: The HMM target length recorded for this gene's ID in
                    the primary (length-filtering) HMM run. None if the ID
                    never appeared there (length unknown / not applicable).
        min_length: The minimum required hit length (--min_length). If
                    None, the length filter is skipped entirely.
 
    Returns:
        Tuple (distance, reason) where distance is a non-negative float
        (0 = would pass) and reason is a short code (or combination of
        codes joined with " + ") describing which filter(s) the gene
        fails: "Exon_count", "Boundary", "Length", or "Confident_level"
        if none of those apply.
    """
 
    exon_count = gene_info['Exon_Count']
    if exon_count < exon_range[0]:
        exon_dist = exon_range[0] - exon_count
    elif exon_count > exon_range[1]:
        exon_dist = exon_count - exon_range[1]
    else:
        exon_dist = 0
 
    pos_range = _effective_pos_range(gene_info, pos_range_kb)
    if gene_info['Strand'] == '+':
        actual_dist = gene_info['Start_Pos'] - gene_info['Element_Start']
    else:
        actual_dist = gene_info['Element_End'] - gene_info['End_Pos']
    pos_dist_bp = max(0.0, actual_dist - pos_range)
    pos_dist_kb = pos_dist_bp / 1000.0
 
    length_dist = 0.0
    if min_length is not None and hit_length is not None and hit_length < min_length:
        length_dist = (min_length - hit_length) / 100.0
 
    total = float(exon_dist) + pos_dist_kb + length_dist
 
    reasons = []
    if exon_dist > 0:
        reasons.append("Exon_count")
    if pos_dist_kb > 0:
        reasons.append("Boundary")
    if length_dist > 0:
        reasons.append("Length")
 
    reason = "+".join(reasons) if reasons else "Confident_level"
 
    return total, reason


def find_closest_candidate(
    candidate_ids: set,
    gff_df: pd.DataFrame,
    exon_range: tuple[int, int],
    pos_range_kb: int,
    length_map: dict | None = None,
    min_length: int | None = None
) -> tuple[str | None, str]:
    """Find the candidate gene closest to passing the acceptance filters.
 
    Among all candidate IDs that have a matching entry in gff_df, returns
    the one with the smallest combined exon/position/length "distance to
    passing" score, along with a short code for why it still fails.
 
    Args:
        candidate_ids: Set of gene IDs to evaluate (e.g. union of every
                        HMM hit at any consensus level for this element,
                        INCLUDING ones dropped by the length filter).
        gff_df: DataFrame as returned by get_gff_data().
        exon_range: Tuple (min, max) accepted exon count.
        pos_range_kb: Distance threshold in kb from element boundary.
        length_map: Optional dict mapping gene ID -> HMM hit length, built
                    from the primary HMM run's raw output. Used to score
                    and report the "Length" failure reason.
        min_length: The minimum required hit length (--min_length).
 
    Returns:
        Tuple (closest_id, reason). closest_id is None if no candidate ID
        has a matching entry in gff_df, in which case reason explains that.
    """
 
    if gff_df.empty or not candidate_ids:
        return None, "No candidate genes available for comparison."
 
    matches = gff_df[gff_df['ID'].isin(candidate_ids)].copy()
    if matches.empty:
        return None, "None of the identified IDs were found in the GFF."
 
    length_map = length_map or {}
 
    best_id = None
    best_reason = ""
    best_dist = None
 
    for _, row in matches.iterrows():
        hit_length = length_map.get(row['ID'])
        dist, reason = _gene_failure_distance(
            row, exon_range, pos_range_kb, hit_length, min_length
        )
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_id = row['ID']
            best_reason = reason
 
    return best_id, best_reason


def process_hmm_files(
    hmm_folder1: str,
    hmm_folder2: str,
    hmm_folder3: str,
    gff_folder: str,
    fasta_folder: str,
    output_file: str,
    empty_output_file: str,
    not_passed_file: str,
    min_common: int,
    min_length: int,
    range_kb: int
) -> None:
    """Orchestrate the full HMM and GFF integration pipeline.

    Iterates over all HMM files in the first folder, resolves consensus
    IDs across three HMMER runs at decreasing stringency levels, validates
    each candidate against GFF filters, and writes accepted IDs to the
    output file. Elements that did not yield an accepted gene are logged
    to empty_output_file (basenames only) and to not_passed_file with the
    rejection reason and the closest-passing candidate gene.

    Args:
        hmm_folder1: Path to the first HMMER results folder.
        hmm_folder2: Path to the second HMMER results folder.
        hmm_folder3: Path to the third HMMER results folder.
        gff_folder: Path to the folder containing GFF files.
        fasta_folder: Path to the folder containing FASTA files.
        output_file: Path to the output file for accepted gene IDs.
        empty_output_file: Path to the log file for elements with no result.
        not_passed_file: Path to the tab-separated log of rejected elements,
                          one line per element: "<Name>\\t<Reason>\\t<Closest gene>".
        min_common: Minimum consensus level required (1, 2, or 3 HMM runs).
        min_length: Minimum hit length to consider from the first HMM run.
        range_kb: Distance threshold in kb from element boundary.

    Raises:
        PermissionError: If the output file cannot be written.
        OSError: If any other I/O error occurs during writing.
    """

    hmm_files1 = glob.glob(os.path.join(hmm_folder1, '*.txt'))
    all_results = []
    
    with open(empty_output_file, 'w') as empty_f, \
         open(not_passed_file, 'w') as not_passed_f:

        for hmm_path1 in hmm_files1:
            base_name = os.path.basename(hmm_path1).split('.')[0]
            print(f"Processing {base_name}...")

            final_selected_id = None
            final_reason = ""

            # ── helper to log a rejection consistently ───────────────────────
            def log_rejection(reason: str, closest_gene: str = "NA") -> None:
                empty_f.write(f"{base_name}\n")
                not_passed_f.write(f"{base_name}\t{reason}\t{closest_gene}\n")

            fasta_path = find_fasta_path(fasta_folder, base_name)
            if not fasta_path:
                final_reason = f"No FASTA file found for {base_name}."
                print(f"No ID selected for {base_name}. Reason: {final_reason}")
                log_rejection("File")
                continue

            gff_df, gff_reason = get_gff_data(os.path.join(gff_folder, f'{base_name}.gff'), fasta_path)
            if gff_df.empty:
                print(f"No ID selected for {base_name}. Reason: {gff_reason}")
                log_rejection("File")
                continue
            
            ids1_full, df1_raw = parse_hmm_file(hmm_path1)
            ids1 = set(df1_raw[df1_raw['Length'] >= min_length]['ID'].values) if not df1_raw.empty else set()
            ids2, _ = parse_hmm_file(os.path.join(hmm_folder2, f'{base_name}.txt'))
            ids3, _ = parse_hmm_file(os.path.join(hmm_folder3, f'{base_name}.txt'))

            id_length_map = {}
            if not df1_raw.empty:
                id_length_map = df1_raw.groupby('ID')['Length'].max().to_dict()
 
            # Pool of every ID identified at any consensus level
            all_seen_ids = ids1_full | ids2 | ids3

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
                    final_reason = g_res
            elif len(common_ids_12) == 1 and min_common <= 2:
                final_selected_id = list(common_ids_12)[0]
                final_reason = "Finalized with unique ID at level 2."
            elif len(common_ids_12) > 1 and min_common <= 2:
                final_selected_id, g_res = run_gff_analysis(common_ids_12, gff_df, range_kb)
                if final_selected_id:
                    final_reason = f"Selected via GFF tie-breaker at level 2. ({g_res})"
                else:
                    final_reason = g_res
            elif len(ids1) == 1 and min_common == 1:
                final_selected_id = list(ids1)[0]
                final_reason = "Finalized with unique ID at level 1."
            elif len(ids1) > 1 and min_common == 1:
                final_selected_id, g_res = run_gff_analysis(ids1, gff_df, range_kb)
                if final_selected_id:
                    final_reason = f"Selected via GFF tie-breaker at level 1. ({g_res})"
                else:
                    final_reason = g_res
            else:
                final_reason = f"No common IDs found at level {min_common}."

            if final_selected_id:
                is_valid, validation_reason = check_filters(final_selected_id, gff_df, EXON_RANGE, range_kb)
                if is_valid:
                    all_results.append({'ID': final_selected_id})
                    print(f"  Success: Selected ID '{final_selected_id}'. Reason: {final_reason}")
                else:
                    print(f"  No ID selected for {base_name}. Reason: {validation_reason[0]}")
                    log_rejection(validation_reason[1], final_selected_id)
            else:
                print(f"  No ID selected for {base_name}. Reason: {final_reason}")
                closest_id, closest_reason = find_closest_candidate(
                    all_seen_ids, gff_df, EXON_RANGE, range_kb,
                    length_map=id_length_map, min_length=min_length
                )
                if closest_id:
                    log_rejection(closest_reason, closest_id)
                else:
                    log_rejection("No_match", "NA")
    
    if all_results:
        try:
            pd.DataFrame(all_results)['ID'].to_csv(output_file, index=False, header=False)
            print(f"\nSuccessfully wrote {len(all_results)} IDs to {output_file}")
        except PermissionError:
            sys.exit(f"Error: No write permission for '{output_file}'.")
        except OSError as e:
            sys.exit(f"Error writing output file: {e}")


def main() -> None:
    """Parse command-line arguments and launch the HMM-GFF pipeline."""

    parser = argparse.ArgumentParser(description='Process HMMSEARCH domtblout and GFF files.')
    parser.add_argument('--hmm1', dest='hmm_folder1', required=True, help='Path to first HMMER folder.')
    parser.add_argument('--hmm2', dest='hmm_folder2', required=True, help='Path to second HMMER folder.')
    parser.add_argument('--hmm3', dest='hmm_folder3', required=True, help='Path to third HMMER folder.')
    parser.add_argument('--gff', dest='gff_folder', required=True, help='Path to GFF folder.')
    parser.add_argument('--fasta', dest='fasta_folder', required=True, help='Path to FASTA folder.')
    parser.add_argument('--output', dest='output_file', required=True, help='Output filename.')
    parser.add_argument('--empty', dest='empty_output_file', required=True, help='Log for empty results.')
    parser.add_argument('--not_passed', dest='not_passed_file', required=True,
                         help='Tab-separated log of rejected elements: "<Name>\\t<Reason>\\t<Closest gene>".')
    parser.add_argument('--min_common', type=int, default=2, choices=[1, 2, 3])
    parser.add_argument('--min_length', type=int, default=250)
    parser.add_argument('--range_kb', type=int, default=20)

    args = parser.parse_args()
    process_hmm_files(args.hmm_folder1, args.hmm_folder2, args.hmm_folder3,
                      args.gff_folder, args.fasta_folder, args.output_file, 
                      args.empty_output_file, args.not_passed_file,
                      args.min_common, args.min_length, args.range_kb)


if __name__ == "__main__":
    main()
    