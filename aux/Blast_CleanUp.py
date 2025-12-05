#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Clean up of blast result by removing reciprocal hits (A-B = B-A), filtering hits by size and percentage identity, merging fragments by interval, and filtering the final element pairs by coverage.")

# Add the arguments to the parser
parser.add_argument("-f", "--file", dest="file_in", required=True,
                    help="Input file. This assumes a file with the following columns: 'qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen'.")
parser.add_argument("-o", "--output", dest="file_out", required=True,
                    help="Output file. Returns the filtered file to the specified location.")
parser.add_argument('-fs', '--fragmentSize', type=int, default=2000, dest='MINIMUM_FRAGMENT_SIZE',
                    help="The minimum fragment size of a blast alignment to be used (default: 2000) [range: 1000, 5000].")
parser.add_argument('-i', '--identity', type=float, default=70.0, dest='MINIMUM_PIDENT',
                    help="The minimum percentage of identity of a blast alignment to be used (default: 70.0) [range: 60.0, 90.0].")
parser.add_argument('-ms', '--mergeSize', type=int, default=2000, dest='MINIMUM_MERGED_SIZE',
                    help="The minimum merge fragment size to be used (default: 5000) [range: 2000, 10000].")
parser.add_argument('-c', '--coverage', type=float, default=20.0, dest='MINIMUM_HIT_COVERAGE',
                    help="The minimum coverage of the filter merge fragments to pass (default: 20.0) [range: 10.0, 50.0].")
args = parser.parse_args()

# Setting minimum filtering values
#MINIMUM_FRAGMENT_SIZE = 2000
#MINIMUM_MERGED_SIZE = 5000
#MINIMUM_HIT_COVERAGE = 20.0
#MINIMUM_PIDENT = 70.0

print("")
print("  Setting hit names")
print(f"  Minimum fragment size: {args.MINIMUM_FRAGMENT_SIZE} bp")
print(f"  Minimum identity percentage: {args.MINIMUM_PIDENT}%")
print(f"  Minimum merged size: {args.MINIMUM_MERGED_SIZE} bp")
print(f"  Minimum coverage: {args.MINIMUM_HIT_COVERAGE}%")

dirty = {}
clean = {}
removed = 0
accepted = 0
selfHit = 0
filtered_by_size = 0
filtered_by_pident = 0
filtered_by_coverage = 0

# This loops through the input file
for line_num, line_content in enumerate(open(args.file_in), 1):
    line = line_content.strip().split()
    
    if len(line) < 11:
        print(f"Warning: Line {line_num} in {args.file_in} has fewer than 11 columns and will be skipped: {line_content.strip()}")
        removed += 1
        continue

    seq1 = line[0]
    seq2 = line[1]
    
    try:
        evalue_val = line[2]
        pident_val = float(line[3])
        bitscore_val = line[4]
        
        qstart_val = int(line[5])
        qend_val = int(line[6])
        qlen_orig = int(line[7])
        
        sstart_val = int(line[8])
        send_val = int(line[9])
        slen_orig = int(line[10])

    except ValueError as e:
        print(f"Error parsing numerical values on line {line_num}: {line_content.strip()}. Error: {e}")
        removed += 1
        continue

    # Normalize coordinates to always be start < end
    qstart_span = min(qstart_val, qend_val)
    qend_span = max(qstart_val, qend_val)
    sstart_span = min(sstart_val, send_val)
    send_span = max(sstart_val, send_val)

    if seq1 == seq2:
        selfHit += 1
        continue
    
    current_qhit_len = qend_span - qstart_span + 1
    current_shit_len = send_span - sstart_span + 1
    
    # Filter by minimum fragment size
    if current_qhit_len < args.MINIMUM_FRAGMENT_SIZE or current_shit_len < args.MINIMUM_FRAGMENT_SIZE:
        filtered_by_size += 1
        continue

     # Filter by minimum pident
    if pident_val < args.MINIMUM_PIDENT:
        filtered_by_pident += 1
        continue

    hit_key = f"{seq1}\t{seq2}"
    reciprocal_hit_key = f"{seq2}\t{seq1}"

    if reciprocal_hit_key in dirty:
        removed += 1 
    
    elif hit_key in dirty:
        fragment_data = {
            'pident': pident_val,
            'q_start': qstart_span,
            'q_end': qend_span,
            's_start': sstart_span,
            's_end': send_span,
            'q_len': qlen_orig,
            's_len': slen_orig,
            'q_hit_len': current_qhit_len,
            's_hit_len': current_shit_len
        }
        dirty[hit_key].append(fragment_data)
    
    else:
        fragment_data = {
            'pident': pident_val,
            'q_start': qstart_span,
            'q_end': qend_span,
            's_start': sstart_span,
            's_end': send_span,
            'q_len': qlen_orig,
            's_len': slen_orig,
            'q_hit_len': current_qhit_len,
            's_hit_len': current_shit_len
        }
        dirty[hit_key] = [fragment_data]

# Process each hit_key to merge fragments and apply filtering criteria
for hit_key in dirty:
    fragments = dirty[hit_key]
    
    qlen_orig = fragments[0]['q_len']
    slen_orig = fragments[0]['s_len']

    # Sort fragments by q_start to handle overlaps
    fragments.sort(key=lambda x: x['q_start'])
    
    merged_q_intervals = []
    
    if fragments:
        # Merge Query intervals
        current_q_start = fragments[0]['q_start']
        current_q_end = fragments[0]['q_end']
        
        for i in range(1, len(fragments)):
            next_q_start = fragments[i]['q_start']
            next_q_end = fragments[i]['q_end']
            
            if next_q_start <= current_q_end:
                current_q_end = max(current_q_end, next_q_end)
            else:
                merged_q_intervals.append((current_q_start, current_q_end))
                current_q_start = next_q_start
                current_q_end = next_q_end
        merged_q_intervals.append((current_q_start, current_q_end))
        
        # Merge Subject intervals
        # Sort fragments by s_start to handle overlaps on the subject
        fragments.sort(key=lambda x: x['s_start'])
        merged_s_intervals = []
        
        current_s_start = fragments[0]['s_start']
        current_s_end = fragments[0]['s_end']
        
        for i in range(1, len(fragments)):
            next_s_start = fragments[i]['s_start']
            next_s_end = fragments[i]['s_end']
            
            if next_s_start <= current_s_end:
                current_s_end = max(current_s_end, next_s_end)
            else:
                merged_s_intervals.append((current_s_start, current_s_end))
                current_s_start = next_s_start
                current_s_end = next_s_end
        merged_s_intervals.append((current_s_start, current_s_end))

    # Filter out individual merged intervals that are too small
    
    filtered_q_intervals = []
    for start, end in merged_q_intervals:
        if (end - start + 1) >= args.MINIMUM_MERGED_SIZE:
            filtered_q_intervals.append((start, end))
        else:
            filtered_by_size += 1

    filtered_s_intervals = []
    for start, end in merged_s_intervals:
        if (end - start + 1) >= args.MINIMUM_MERGED_SIZE:
            filtered_s_intervals.append((start, end))
        else:
            filtered_by_size += 1
    
    # If there are no valid merged intervals left after filtering, skip this hit pair
    if not filtered_q_intervals or not filtered_s_intervals:
        continue

    # Calculate total lengths from the 'filtered' merged intervals
    total_q_hit_len_sum = sum(end - start + 1 for start, end in filtered_q_intervals)
    total_s_hit_len_sum = sum(end - start + 1 for start, end in filtered_s_intervals)
    
    # Calculate weighted pident from the remaining fragments
    total_weighted_pident_sum = sum(frag['pident'] * frag['s_hit_len'] for frag in fragments)
    final_pident = round(total_weighted_pident_sum / sum(frag['s_hit_len'] for frag in fragments), 3) if sum(frag['s_hit_len'] for frag in fragments) > 0 else 0.0

    final_combined_hit = {
        'qlen': qlen_orig,
        'slen': slen_orig,
        'pident': final_pident,
        'total_q_hit_len_sum': total_q_hit_len_sum,
        'total_s_hit_len_sum': total_s_hit_len_sum
    }

    if final_combined_hit:
        q_cover = (100.0 * final_combined_hit['total_q_hit_len_sum'] / final_combined_hit['qlen']) if final_combined_hit['qlen'] > 0 else 0.0
        s_cover = (100.0 * final_combined_hit['total_s_hit_len_sum'] / final_combined_hit['slen']) if final_combined_hit['slen'] > 0 else 0.0

        if q_cover < args.MINIMUM_HIT_COVERAGE and s_cover < args.MINIMUM_HIT_COVERAGE:
            filtered_by_coverage += 1
            continue

        clean[hit_key] = final_combined_hit
        accepted += 1
        
print("  Writing cleaned file with '", accepted, "' hits", sep="")

with open(args.file_out, "w") as outfile:
    for hit_key in clean.keys():
        seq1, seq2 = hit_key.split('\t')
        
        output_line_elements = [
            seq1,
            seq2,
            str(clean[hit_key]['pident']),
            str(clean[hit_key]['qlen']),
            str(clean[hit_key]['slen']),
            str(clean[hit_key]['total_q_hit_len_sum']),
            str(clean[hit_key]['total_s_hit_len_sum'])
        ]
        print('\t'.join(output_line_elements), file=outfile)
        
print("    Removed hits:\t", removed)
print("    Self matching:\t", selfHit)
print("    Filtered by min fragment size:\t", filtered_by_size)
print("    Filtered by min identity:\t", filtered_by_pident)
print("    Filtered by min coverage:\t", filtered_by_coverage)
print("Done")
print("")
