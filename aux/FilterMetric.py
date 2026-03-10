#!/usr/bin/env python3
"""
Metric filtering approach of synteny results.

Process MSCANX results against original DIAMOND data and filter pairs
with metric values below a threshold.
"""

import os
import csv
import glob
import pandas as pd
import argparse
import sys


# CONSTANTS
PIDENT_THRESHOLD = 60.0
PIDENT_THRESHOLD_BONUS = 95.0
MIN_LENGTH = 100


def get_pairs(file_path: str) -> list[list]:
    """Process the pair file.

    Read and store the pair information

    Args:
        file_path: Path to the file that contain the syntenet pair information.

    Returns:
        list (pair_info) where:
            - pair_info: List of the pair information.
                         Each element is [element01, element02, GCP, ECP01, ECP02].

    Raises:
        FileNotFoundError: If the file do not exist.
        PermissionError: If the file does not have reading permission.
    """
    
    full_info = []
    try:
        with open(file_path, 'r') as f:
            reader = csv.reader(f, delimiter =";")
            for row in reader:
                if reader.line_num == 1:
                    continue
                full_info.append([row[0].strip(),row[1].strip(),float(row[2].strip()),float(row[3].strip()),float(row[4].strip())])
    except FileNotFoundError:
        sys.exit(f"Error: File not found: {file_path}") 
    except PermissionError:
        sys.exit(f"Error: No read permission for: {file_path}") 

    return full_info


def get_anchors(pairs: list[list],folder_path: str) -> dict[tuple, list]:
    """Process the pair file.

    Read and store the pair information

    Args:
        pairs: List of the pair information.
               Each element is [element01, element02, GCP, ECP01, ECP02].
        folder_path: Path to the folder where .collinearity files are stored.

    Returns:
        dict (pair, anchors) where:
            - pair: Pair information as (element01, element02).
            - anchors: List of anchor pairs from the .collinearity file.
                       Each elements is [Gene01, Gene02]

    Raises:
        FileNotFoundError: If the file do not exist.
    """
    
    anchor_pairs = dict()
    for pair in pairs:
        if pair[0].lower() < pair[1].lower():
            filename = pair[0] + "_" + pair[1] + ".collinearity"
        else:
            filename = pair[1] + "_" + pair[0] + ".collinearity"
        
        anchors = []
        try:
            with open(os.path.join(folder_path, f'{filename}'), 'r') as file:
                for line in file:
                    if line.startswith('#'):
                        continue
                    parts = line.split("\t")
                    if len(parts) == 4:
                        anchors.append(sorted([parts[1],parts[2]]))
        
            anchor_pairs.update({(pair[0],pair[1]) : anchors})  
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {filename}.tsv.")
    
    return anchor_pairs


def get_diamond_results(anchors: dict[tuple, list],folder_path: str) -> dict[tuple, list]:
    """Process the pair file.

    Read and store the pair information

    Args:
        anchors: Dictionary with anchor list per pair.
        folder_path: Path to the folder where DIAMOND results are stored.

    Returns:
        dict (pair, matches) where:
            - pair: Pair information as (element01, element02).
            - matches: List of DIAMOND matches filtered by threshold.
                       Each element is [qseqid, sseqid, pident, length].

    Raises:
        FileNotFoundError: If the file do not exist.
    """
    
    diamond_pairs = dict()
    for pair in anchors.keys():
        matches = []
        
        filename = pair[0] + "_" + pair[1]
        try:
            with open(os.path.join(folder_path, f'{filename}.tsv'), 'r') as file:
                for line in file:
                    parts = line.split()
                    if len(parts) == 12:
                        if (sorted([parts[0],parts[1]])) in anchors.get(pair) and float(parts[2]) >= PIDENT_THRESHOLD and int(parts[3]) >= MIN_LENGTH:
                            matches.append([parts[0],parts[1], float(parts[2]), int(parts[3])])
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {filename}.tsv.")   

        filename = pair[1] + "_" + pair[0]
        try:
            with open(os.path.join(folder_path, f'{filename}.tsv'), 'r') as file:
                for line in file:
                    parts = line.split()
                    if len(parts) == 12:
                        if (sorted([parts[0],parts[1]])) in anchors.get(pair) and float(parts[2]) >= PIDENT_THRESHOLD and int(parts[3]) >= MIN_LENGTH:
                            matches.append([parts[0],parts[1], float(parts[2]), int(parts[3])])
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {filename}.tsv.")
        
        diamond_pairs.update({pair : matches})  
    
    return diamond_pairs
    

def process_metric(anchors: list, diamond: list) -> tuple[float, float]:
    """Calculate the quality metric for a syntenic pair.

    Combines the number of DIAMOND matches with a bonus for high
    sequence identity, normalized over the total number of anchors.

    Args:
        anchors: List of anchor pairs from the .collinearity file.
        diamond: List of DIAMOND matches filtered by threshold.
                 Each element is [qseqid, sseqid, pident, length].

    Returns:
        Tuple (final_metric, ratio) where:
            - final_metric: raw metric value (base + bonus).
            - ratio: normalized value between 0 and 1.

    Raises:
        ZeroDivisionError: If anchors is empty.
    """
    
    anchor_number = len(anchors)
    
    base_value = len(diamond)
    
    bonus = 0.0
    
    for match in diamond:
        if float(match[2]) >= PIDENT_THRESHOLD_BONUS:
            bonus += int(match[3])*0.2
    
    final_metric = base_value + bonus
    try:
        ratio = min(1, final_metric/anchor_number)
    except ZeroDivisionError:
        raise ZeroDivisionError(f"Cannot compute ratio: anchor list is empty for this pair.")
    return final_metric, ratio


def process_filter(collinearity_file: str,
                   syntenet_folder: str,
                   diamond_folder: str,
                   threshold: int,
                   output_file: str
) -> None:
    """Process the filter process

    Manage all the information and process each step of the process

    Args:
        collinearity_file: Path to the file that contain the syntenet pair information.
        syntenet_folder: Path to the folder where .collinearity files are stored.
        diamond_folder: Path to the folder where DIAMOND results are stored.
        threshold: Minimum quality metric value to accept a pair.
        output_file: Path to the file where the accepted pair information is going to be stored.

    Returns:
        None
    """
    
    input_info = get_pairs(collinearity_file)
    
    anchors = get_anchors(input_info,syntenet_folder)
    
    diamond = get_diamond_results(anchors,diamond_folder)
    
    all_results = []
    
    for pair in input_info:
        metric, ratio = process_metric(anchors.get((pair[0],pair[1])) ,diamond.get((pair[0],pair[1])))
        if metric >= threshold:
            all_results.append([pair[0],pair[1],pair[2]*ratio,pair[3]*ratio,pair[4]*ratio])
    
    if all_results:
        df_results = pd.DataFrame(all_results, columns=['element01','element02','General_percentage','element01_percentage','element02_percentage'])
        try:
            df_results.to_csv(output_file, index=False, header=True, sep=';',decimal='.', float_format="%.2f")
            print(f"\nSuccessfully wrote {len(all_results)} pairs to {output_file}")
        except PermissionError:
            sys.exit(f"Error: No write permission for '{output_file}'.")
        except OSError as e:
            sys.exit(f"Error writing output file: {e}")
      

def main() -> None:
    """Parse command-line arguments and launch the filtering pipeline."""

    parser = argparse.ArgumentParser(description='Process collinearity results, filter and update results.')
    parser.add_argument('--collinearity', dest='collinearity', required=True, help='Path to collinearity results.')
    parser.add_argument('--syntenet', dest='syntenet', required=True, help='Path to results of syntenet')
    parser.add_argument('--diamond', dest='diamond', required=True, help='Path to DIAMOND results.')
    parser.add_argument('--threshold', dest='threshold', required=True, type=int, help='Minimum metric value to pass')
    parser.add_argument('--output', dest='output_file', required=False, default = './output.csv', help='Path to output file.')
    
    args = parser.parse_args()
    process_filter(args.collinearity,args.syntenet,args.diamond,args.threshold,args.output_file)
                   
 
if __name__ == "__main__":
    main()