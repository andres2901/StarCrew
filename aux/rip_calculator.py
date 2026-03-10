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


def calculate_rip_affected_bp(
    sequence: str,
    window_size: int,
    step_size: int,
    comp_thresh: float,
    prod_thresh: float,
    subst_thresh: float
) -> tuple[int, int, float]:
    """Calculate the number of base pairs affected by RIP mutations.

    Slides a window over the sequence and evaluates three RIP indices
    per window. Marks the intersection region of three consecutive
    positive windows as RIP-affected.

    Args:
        sequence: Nucleotide sequence to analyze.
        window_size: Length of the sliding window in base pairs.
        step_size: Number of base pairs to advance between windows.
        comp_thresh: Minimum composite index value for a positive window.
        prod_thresh: Minimum product index value for a positive window.
        subst_thresh: Maximum substrate index value for a positive window.

    Returns:
        Tuple (total_length, affected_bp, percentage) where:
            - total_length: Total length of the input sequence.
            - affected_bp: Number of base pairs flagged as RIP-affected.
            - percentage: Percentage of the sequence that is RIP-affected.
    """

    sequence = sequence.upper()
    seq_len = len(sequence)
    affected_bases = [False] * seq_len
    window_results = []
    
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
        
        prod_idx = c_TA / c_AT if c_AT > 0 else 0.0
        subs_denom = c_AC + c_GT
        subs_idx = (c_CA + c_TG) / subs_denom if subs_denom > 0 else 2.0 # Default high to fail thresh
        comp_idx = prod_idx - subs_idx
        
        is_pos = (comp_idx > comp_thresh and 
                  prod_idx > prod_thresh and 
                  subs_idx < subst_thresh)
        
        window_results.append({'start': start, 'end': end, 'is_pos': is_pos})
        
    for i in range(len(window_results) - 2):
        w1, w2, w3 = window_results[i], window_results[i+1], window_results[i+2]
        
        if w1['is_pos'] and w2['is_pos'] and w3['is_pos']:
            intersect_start = w3['start']
            intersect_end = min(w1['end'], w2['end'], w3['end'])
            
            if intersect_start < intersect_end:
                for j in range(intersect_start, intersect_end):
                    affected_bases[j] = True
                    
    total_affected = sum(affected_bases)
    pct = (total_affected / seq_len) * 100 if seq_len > 0 else 0.0
    return seq_len, total_affected, pct


def main() -> None:
    """Parse command-line arguments and launch the RIP detection pipeline."""
    
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
            
    except FileNotFoundError:
        sys.exit(f"Error: FASTA file '{args.fasta}' not found.")
    except PermissionError:
        sys.exit(f"Error: No read permission for '{args.fasta}'.")
    except ValueError as e:
        sys.exit(f"Error parsing FASTA file '{args.fasta}': {e}")
    except OSError as e:
        sys.exit(f"Error reading file '{args.fasta}': {e}")


if __name__ == "__main__":
    main()