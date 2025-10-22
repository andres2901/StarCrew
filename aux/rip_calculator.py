import argparse
from Bio import SeqIO
import sys

def calculate_consecutive_rip_percentage(sequence, window_size, step_size, comp_thresh, prod_thresh, subst_thresh):
    """
    Calculates the percentage of RIP-affected sequence, counting ONLY the overlap 
    region of THREE consecutive, overlapping windows that satisfy the triple thresholds.
    """
    sequence = sequence.upper()
    sequence_length = len(sequence)
    
    # Array to track which base pairs are covered by consecutive positive windows
    affected_bases = [False] * sequence_length 
    
    # Store the result of the threshold check for all windows
    window_results = []
    
    # --- Step 1: Analyze All Windows and Store Results ---
    
    # Iterate through the sequence with the specified step size
    for start in range(0, sequence_length - window_size + 1, step_size):
        end = start + window_size
        window = sequence[start:end]
        
        # Observed Dinucleotide Counts
        count_TA = window.count('TA')
        count_AT = window.count('AT')
        count_CA = window.count('CA')
        count_TG = window.count('TG')
        count_AC = window.count('AC')
        count_GT = window.count('GT')
        
        # --- Index Calculation ---
        rip_product_index = count_TA / count_AT if count_AT > 0 else float('nan')
        substrate_denominator = count_AC + count_GT
        rip_substrate_index = (count_CA + count_TG) / substrate_denominator if substrate_denominator > 0 else float('nan')
            
        rip_composite_index = float('nan')
        if not (rip_product_index == float('nan') or rip_substrate_index == float('nan')):
            rip_composite_index = rip_product_index - rip_substrate_index
        
        # --- Triple Threshold Check ---
        is_positive = (
            not rip_composite_index == float('nan') and 
            rip_composite_index > comp_thresh and           # Composite > Comp Threshold
            rip_product_index > prod_thresh and             # Product > Prod Threshold (Enrichment)
            rip_substrate_index < subst_thresh              # Substrate < Substrate Threshold (Depletion)
        )
        
        # Store the window result along with its coordinates
        window_results.append({
            'start': start, 
            'end': end,         
            'is_positive': is_positive
        })
        
    # --- Step 2: Identify CONSECUTIVE TRIPLE Positive Windows and Mark ONLY THE OVERLAP ---
    
    # Iterate up to the third-to-last window to allow for a triple check (i, i+1, i+2)
    for i in range(len(window_results) - 2):
        w1 = window_results[i]
        w2 = window_results[i+1]
        w3 = window_results[i+2]
        
        # Condition: ALL THREE consecutive windows must be positive
        if w1['is_positive'] and w2['is_positive'] and w3['is_positive']:
            
            # The affected region is the INTERSECTION (overlap) of the three windows.
            
            # Intersection Start: The latest starting coordinate (i.e., the start of w3)
            affected_start = w3['start']   
            
            # Intersection End: The earliest ending coordinate (i.e., the end of w1)
            affected_end = w1['end']       
            
            # Example: W1 (0-1000), W2 (500-1500), W3 (1000-2000)
            # Intersection Start (w3['start']): 1000
            # Intersection End (w1['end']): 1000
            # This logic needs adjustment for the correct overlap region of three windows.
            # 
            # Correct intersection of W1(S1, E1), W2(S2, E2), W3(S3, E3) is:
            # Start: max(S1, S2, S3) -> S3
            # End: min(E1, E2, E3) -> E1
            
            # Recalculating affected_end to find the EARLIEST end coordinate among the three
            affected_end = min(w1['end'], w2['end'], w3['end'])
            
            # Mark all base pairs within this intersection region as True
            if affected_start < affected_end:
                for j in range(affected_start, affected_end):
                    if j < sequence_length:
                        affected_bases[j] = True
                    
    # --- Step 3: Count Uniquely Affected Base Pairs ---
    total_affected_bp = sum(affected_bases)
    
    # Calculate percentage
    rip_percentage = (total_affected_bp / sequence_length) * 100 if sequence_length > 0 else 0.0
    
    return sequence_length, total_affected_bp, rip_percentage

# --- Main Execution Block ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate the percentage of RIP-affected sequence using a triple threshold (Composite, Product, Substrate) and counting ONLY the overlap of THREE consecutive positive windows.")
    
    parser.add_argument('fasta_file', help="Input FASTA file path.")
    parser.add_argument('-w', '--window', type=int, default=1000, help="Sliding window size (bp). Default: 1000")
    parser.add_argument('-s', '--step', type=int, default=500, help="Sliding window step size (bp). Default: 500")
    
    # Triple Threshold Arguments
    parser.add_argument('-tc', '--comp_thresh', type=float, default=0.0, help="Minimum Composite Index threshold (Product - Substrate). Default: 0.0")
    parser.add_argument('-tp', '--prod_thresh', type=float, default=0.86, help="Minimum Product Index threshold (TpA/ApT). Default: 0.86")
    parser.add_argument('-ts', '--subst_thresh', type=float, default=1.03, help="Maximum Substrate Index threshold ((CpA+TpG)/(ApC+GpT)). Default: 1.03")
    
    args = parser.parse_args()
    
    # Header for the output table
    print(f"Scaffold_ID\tLength (bp)\tRIP_Triple_Overlap_BP\tRIP_Triple_Overlap_Percentage (%)")
    print("-" * 85)

    try:
        # Read the FASTA file using Biopython
        for record in SeqIO.parse(args.fasta_file, "fasta"):
            sequence = str(record.seq)
            
            # Skip sequences shorter than the windows required
            min_length = args.window + 2 * args.step
            if len(sequence) < min_length:
                print(f"Warning: Sequence '{record.id}' length ({len(sequence)} bp) is too short for 3 consecutive windows. Skipping.", file=sys.stderr)
                print(f"{record.id}\t{len(sequence)}\t0\t0.00")
                continue
            
            # Calculate affected base pairs and percentage
            sequence_length, total_affected_bp, rip_percentage = \
                calculate_consecutive_rip_percentage(sequence, args.window, args.step, args.comp_thresh, args.prod_thresh, args.subst_thresh)
            
            # Print final result
            print(f"{record.id}\t{sequence_length}\t{total_affected_bp}\t{rip_percentage:.2f}")
                
    except FileNotFoundError:
        print(f"Error: The file '{args.fasta_file}' was not found.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
