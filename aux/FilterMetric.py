#!/usr/bin/env python3
"""
Metric filteric approach of synteny

Process MSCANX results against original DIAMOND data and filter pairs
with metric values below a threshold.
"""

import os
import csv
import glob
import pandas as pd
import argparse

# DIAMOND filter threshold

pident = 60.0
pident_bonus = 95.0
length = 100

def get_pairs(file_path):
    """Identify the pairs to analyze"""
    
    Collinearity_pairs = []
    full_info = []
    with open(file_path, 'r') as f:
        reader = csv.reader(f, delimiter =";")
        for row in reader:
            if reader.line_num == 1:
                continue
            full_info.append([row[0].strip(),row[1].strip(),float(row[2].strip()),float(row[3].strip()),float(row[4].strip())])
    
    return full_info

def get_anchors(Pairs,folder_path):
    """Extract anchor points per pair"""
    
    anchor_pairs = dict()
    for pair in Pairs:
        if pair[0].lower() < pair[1].lower():
            file = pair[0] + "_" + pair[1] + ".collinearity"
        else:
            file = pair[1] + "_" + pair[0] + ".collinearity"
        
        anchors = []
        with open(os.path.join(folder_path, f'{file}'), 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                parts = line.split("\t")
                if len(parts) == 4:
                    anchors.append(sorted([parts[1],parts[2]]))
        
            anchor_pairs.update({(pair[0],pair[1]) : anchors})  
    
    return anchor_pairs

def get_diamond_results(anchors,folder_path):
    """Extract DIAMOND results per pair"""
    
    Diamond_pairs = dict()
    for pair in anchors.keys():
        matches = []
        
        file = pair[0] + "_" + pair[1]
        with open(os.path.join(folder_path, f'{file}.tsv'), 'r') as file:
            for line in file:
                parts = line.split()
                if len(parts) == 12:
                    if (sorted([parts[0],parts[1]])) in anchors.get(pair) and float(parts[2]) >= pident and int(parts[3]) >= length:
                        matches.append([parts[0],parts[1], parts[2], parts[3]])
                
        file = pair[1] + "_" + pair[0]
        with open(os.path.join(folder_path, f'{file}.tsv'), 'r') as file:
            for line in file:
                parts = line.split()
                if len(parts) == 12:
                    if (sorted([parts[0],parts[1]])) in anchors.get(pair) and float(parts[2]) >= pident and int(parts[3]) >= length:
                        matches.append([parts[0],parts[1], float(parts[2]), int(parts[3])])
        
        Diamond_pairs.update({pair : matches})  
    
    return Diamond_pairs
    
def process_metric(anchors, diamond):
    """Process the metric"""
    
    anchor_number = len(anchors)
    
    base_value = len(diamond)
    
    bonus = float()
    
    for match in diamond:
        if float(match[2]) >= pident_bonus:
            bonus += int(match[3])*0.2
    
    final_metric = base_value + bonus
    ratio = min(1, final_metric/anchor_number)
    
    return final_metric, ratio


def process_filter(collinearity_file,syntenet_folder,diamond_folder,threshold,output_file):
    """Process the collinerity pairs and return the filter output"""
    
    input_info = get_pairs(collinearity_file)
    
    anchors = get_anchors(input_info,syntenet_folder)
    
    Diamond = get_diamond_results(anchors,diamond_folder)
    
    all_results = []
    
    for pair in input_info:
        Metric, ratio = process_metric(anchors.get((pair[0],pair[1])) ,Diamond.get((pair[0],pair[1])))
        if Metric >= threshold:
            all_results.append([pair[0],pair[1],pair[2]*ratio,pair[3]*ratio,pair[4]*ratio])
    
    if all_results:
        df_results = pd.DataFrame(all_results, columns=['element01','element02','General_percentage','element01_percentage','element02_percentage'])
        df_results.to_csv(output_file, index=False, header=True, sep=';',decimal='.', float_format="%.2f")
        print(f"\nSuccessfully wrote {len(all_results)} pairs to {output_file}")
      
def main():
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