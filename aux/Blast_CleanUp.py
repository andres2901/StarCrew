#!/usr/bin/env python3
"""
BLAST result cleanup and reciprocal hit removal.

This script processes BLAST output to remove reciprocal hits, filter by
identity and fragment size, merge overlapping intervals, and validate
the final hit coverage.
"""

import argparse
import sys


def process_blast_hits(args):
    """
    Filter and merge BLAST hits based on user-defined thresholds.

    Args:
        args (argparse.Namespace): Validated command-line arguments.

    Returns:
        dict: A dictionary of cleaned hits ready for output.
        dict: A dictionary containing statistics of the filtering process.
    """
    dirty = {}
    clean = {}
    
    stats = {
        "removed": 0,
        "accepted": 0,
        "self_hit": 0,
        "filtered_size": 0,
        "filtered_pident": 0,
        "filtered_coverage": 0
    }

    # Loop through the input file
    try:
        with open(args.file_in, "r") as f_in:
            for line_num, line_content in enumerate(f_in, 1):
                line = line_content.strip().split()

                if len(line) < 11:
                    print(f"Warning: Line {line_num} skipped (insufficient columns).")
                    stats["removed"] += 1
                    continue

                seq1, seq2 = line[0], line[1]

                try:
                    pident_val = float(line[3])
                    qstart_val, qend_val, qlen_orig = int(line[5]), int(line[6]), int(line[7])
                    sstart_val, send_val, slen_orig = int(line[8]), int(line[9]), int(line[10])
                except ValueError as e:
                    print(f"Error parsing line {line_num}: {e}")
                    stats["removed"] += 1
                    continue

                # Normalize coordinates
                q_span = (min(qstart_val, qend_val), max(qstart_val, qend_val))
                s_span = (min(sstart_val, send_val), max(sstart_val, send_val))

                if seq1 == seq2:
                    stats["self_hit"] += 1
                    continue

                curr_q_len = q_span[1] - q_span[0] + 1
                curr_s_len = s_span[1] - s_span[0] + 1

                if curr_q_len < args.MINIMUM_FRAGMENT_SIZE or curr_s_len < args.MINIMUM_FRAGMENT_SIZE:
                    stats["filtered_size"] += 1
                    continue

                if pident_val < args.MINIMUM_PIDENT:
                    stats["filtered_pident"] += 1
                    continue

                hit_key = f"{seq1}\t{seq2}"
                reciprocal_key = f"{seq2}\t{seq1}"

                if reciprocal_key in dirty:
                    stats["removed"] += 1
                else:
                    fragment_data = {
                        'pident': pident_val, 'q_start': q_span[0], 'q_end': q_span[1],
                        's_start': s_span[0], 's_end': s_span[1], 'q_len': qlen_orig,
                        's_len': slen_orig, 'q_hit_len': curr_q_len, 's_hit_len': curr_s_len
                    }
                    if hit_key not in dirty:
                        dirty[hit_key] = []
                    dirty[hit_key].append(fragment_data)

    except FileNotFoundError:
        sys.exit(f"Error: File {args.file_in} not found.")

    # Process intervals and merge
    for hit_key, fragments in dirty.items():
        qlen_orig = fragments[0]['q_len']
        slen_orig = fragments[0]['s_len']

        # Merge Query Intervals
        fragments.sort(key=lambda x: x['q_start'])
        merged_q = merge_intervals(fragments, 'q_start', 'q_end')

        # Merge Subject Intervals
        fragments.sort(key=lambda x: x['s_start'])
        merged_s = merge_intervals(fragments, 's_start', 's_end')

        # Filter intervals by size
        filtered_q = [i for i in merged_q if (i[1] - i[0] + 1) >= args.MINIMUM_MERGED_SIZE]
        filtered_s = [i for i in merged_s if (i[1] - i[0] + 1) >= args.MINIMUM_MERGED_SIZE]

        if not filtered_q or not filtered_s:
            stats["filtered_size"] += 1
            continue

        total_q_sum = sum(end - start + 1 for start, end in filtered_q)
        total_s_sum = sum(end - start + 1 for start, end in filtered_s)

        total_weighted_pident = sum(f['pident'] * f['s_hit_len'] for f in fragments)
        total_s_hit_len = sum(f['s_hit_len'] for f in fragments)
        final_pident = round(total_weighted_pident / total_s_hit_len, 3) if total_s_hit_len > 0 else 0.0

        q_cover = (100.0 * total_q_sum / qlen_orig) if qlen_orig > 0 else 0.0
        s_cover = (100.0 * total_s_sum / slen_orig) if slen_orig > 0 else 0.0

        if q_cover < args.MINIMUM_HIT_COVERAGE and s_cover < args.MINIMUM_HIT_COVERAGE:
            stats["filtered_coverage"] += 1
            continue

        clean[hit_key] = {
            'pident': final_pident, 'qlen': qlen_orig, 'slen': slen_orig,
            'total_q_hit_len_sum': total_q_sum, 'total_s_hit_len_sum': total_s_sum
        }
        stats["accepted"] += 1

    return clean, stats


def merge_intervals(fragments, start_attr, end_attr):
    """Function to merge overlapping numeric intervals."""
    if not fragments:
        return []
    merged = []
    curr_start = fragments[0][start_attr]
    curr_end = fragments[0][end_attr]

    for i in range(1, len(fragments)):
        if fragments[i][start_attr] <= curr_end:
            curr_end = max(curr_end, fragments[i][end_attr])
        else:
            merged.append((curr_start, curr_end))
            curr_start = fragments[i][start_attr]
            curr_end = fragments[i][end_attr]
    merged.append((curr_start, curr_end))
    return merged


def main():
    """Main CLI execution entry point."""
    parser = argparse.ArgumentParser(
        description="Clean up BLAST results by removing reciprocal hits and filtering."
    )
    parser.add_argument("-f", "--file", dest="file_in", required=True, help="Input BLAST file")
    parser.add_argument("-o", "--output", dest="file_out", required=True, help="Output file")
    parser.add_argument("-fs", "--fragmentSize", type=int, default=2000, dest="MINIMUM_FRAGMENT_SIZE")
    parser.add_argument("-i", "--identity", type=float, default=70.0, dest="MINIMUM_PIDENT")
    parser.add_argument("-ms", "--mergeSize", type=int, default=2000, dest="MINIMUM_MERGED_SIZE")
    parser.add_argument("-c", "--coverage", type=float, default=20.0, dest="MINIMUM_HIT_COVERAGE")

    args = parser.parse_args()

    print(f"\n  Processing: {args.file_in}")
    print(f"  Min Fragment: {args.MINIMUM_FRAGMENT_SIZE} bp | Min Identity: {args.MINIMUM_PIDENT}%")

    clean_hits, stats = process_blast_hits(args)

    with open(args.file_out, "w") as outfile:
        for hit_key, data in clean_hits.items():
            seq1, seq2 = hit_key.split('\t')
            row = [seq1, seq2, str(data['pident']), str(data['qlen']), 
                   str(data['slen']), str(data['total_q_hit_len_sum']), 
                   str(data['total_s_hit_len_sum'])]
            print('\t'.join(row), file=outfile)

    print(f"  Writing cleaned file with {stats['accepted']} hits")
    print(f"    Removed reciprocal: {stats['removed']}")
    print(f"    Self matching:      {stats['self_hit']}")
    print(f"    Filtered by size:   {stats['filtered_size']}")
    print(f"    Filtered by ident:  {stats['filtered_pident']}")
    print(f"    Filtered by cover:  {stats['filtered_coverage']}")
    print("Done\n")


if __name__ == "__main__":
    main()
    