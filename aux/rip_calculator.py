#!/usr/bin/env python3
"""
RIP Triple-Threshold Window Tool.

Calculates the percentage of a sequence affected by Repeat-Induced Point (RIP) 
mutations using a sliding window approach and a triple-threshold overlap logic.
"""

import argparse
from Bio import SeqIO
import sys
import os

def calculate_rip_affected_bp(sequence, window_size, step_size, comp_thresh, prod_thresh, subst_thresh):
    """
    Calculates affected BP by checking for 3 consecutive positive windows.
    Returns: (total_length, affected_bp, percentage)
    """
    sequence = sequence.upper()
    seq_len = len(sequence)
    affected_bases = [False] * seq_len
    window_results = []
    
    # Step 1: Analyze Windows
    for start in range(0, seq_len - window_size + 1, step_size):
        end = start + window_size
        subseq = sequence[start:end]
        
        # Dinucleotide counts
        c_TA = subseq.count('TA')
        c_AT = subseq.count('AT')
        c_CA = subseq.count('CA')
        c_TG = subseq.count('TG')
        c_AC = subseq.count('AC')
        c_GT = subseq.count('GT')
        
        # Calculate Indices
        prod_idx = c_TA / c_AT if c_AT > 0 else 0.0
        subs_denom = c_AC + c_GT
        subs_idx = (c_CA + c_TG) / subs_denom if subs_denom > 0 else 2.0 # Default high to fail thresh
        comp_idx = prod_idx - subs_idx
        
        # Check Thresh
        is_pos = (comp_idx > comp_thresh and 
                  prod_idx > prod_thresh and 
                  subs_idx < subst_thresh)
        
        window_results.append({'start': start, 'end': end, 'is_pos': is_pos})
        
    # Step 2: Triple Overlap Logic
    for i in range(len(window_results) - 2):
        w1, w2, w3 = window_results[i], window_results[i+1], window_results[i+2]
        
        if w1['is_pos'] and w2['is_pos'] and w3['is_pos']:
            # The intersection starts at the start of the 3rd window
            # and ends at the end of the 1st window.
            intersect_start = w3['start']
            intersect_end = min(w1['end'], w2['end'], w3['end'])
            
            if intersect_start < intersect_end:
                for j in range(intersect_start, intersect_end):
                    affected_bases[j] = True
                    
    total_affected = sum(affected_bases)
    pct = (total_affected / seq_len) * 100 if seq_len > 0 else 0.0
    return seq_len, total_affected, pct

def main():
    """CLI Entry Point."""
    parser = argparse.ArgumentParser(description="RIP detection with triple-consecutive window overlap.")
    parser.add_argument('fasta', help="Input FASTA file")
    parser.add_argument('-w', '--window', type=int, default=1000, help="Window size (1000)")
    parser.add_argument('-s', '--step', type=int, default=250, help="Step size (250)")
    parser.add_argument('-tc', '--comp', type=float, default=0.0, help="Min Composite Index")
    parser.add_argument('-tp', '--prod', type=float, default=0.8, help="Min Product Index")
    parser.add_argument('-ts', '--subst', type=float, default=1.05, help="Max Substrate Index")

    args = parser.parse_args()

    if (args.step * 3) > (args.window - 100):
        sys.exit(f"Error: Step size ({args.step}) is too large for window size ({args.window}). "
                 f"(Step * 3) must be <= (Window - 100).")

    if not os.path.exists(args.fasta):
        sys.exit(f"Error: {args.fasta} not found.")

    print("ID\tLength\tRIP_BP\tRIP_Pct")
    
    try:
        for record in SeqIO.parse(args.fasta, "fasta"):
            seq_str = str(record.seq)
            # Threshold to ensure 3 windows are even possible
            if len(seq_str) < (args.window + 2 * args.step):
                print(f"{record.id}\t{len(seq_str)}\t0\t0.00")
                continue
                
            length, bp, pct = calculate_rip_affected_bp(
                seq_str, args.window, args.step, args.comp, args.prod, args.subst
            )
            print(f"{record.id}\t{length}\t{bp}\t{pct:.2f}")
            
    except Exception as e:
        sys.exit(f"Processing Error: {e}")

if __name__ == "__main__":
    main()