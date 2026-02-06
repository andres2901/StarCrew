#!/usr/bin/env python3
"""
Metadata Enrichment Tool.

Merges a primary data file with metadata for two distinct columns, 
interleaves the resulting metadata columns, and handles missing values.
"""

import pandas as pd
import argparse
import sys
import os

def merge_and_process_files(data_file_path, metadata_file_path, output_file_path):
    """
    Performs dual-column enrichment and interleaved reordering.
    """
    try:
        data_df = pd.read_csv(
            data_file_path, 
            sep=';', 
            dtype={'element01': str, 'element02': str}
        )
        
        # Metadata file: Auto-detect separator, treat all as strings
        try:
            metadata_df = pd.read_csv(metadata_file_path, sep=None, engine='python', dtype=str)
        except Exception:
            metadata_df = pd.read_csv(metadata_file_path, sep='\t', dtype=str)

        # 2. Setup Join Keys
        if len(metadata_df.columns) < 2:
            sys.exit("Error: Metadata file must have at least two columns.")
            
        # Standardize the join key name (second column)
        join_key = "ElementID_updated"
        metadata_df.columns.values[1] = join_key

        # Identify columns to be added (everything from index 2 onwards)
        meta_payload_cols = metadata_df.columns[2:].tolist()
        
        if not meta_payload_cols:
            print("Warning: No metadata columns found to merge. Saving raw data as TSV.")
            data_df.to_csv(output_file_path, sep='\t', index=False)
            return

        # 3. Sequential Merging
        # We perform two left-joins to enrich both element01 and element02
        
        # Merge for element01
        meta01 = metadata_df[[join_key] + meta_payload_cols].copy()
        meta01.columns = [join_key] + [f"element01_{c}" for c in meta_payload_cols]
        
        merged_df = pd.merge(
            data_df, meta01, left_on='element01', right_on=join_key, how='left'
        ).drop(columns=[join_key])
        
        # Merge for element02
        meta02 = metadata_df[[join_key] + meta_payload_cols].copy()
        meta02.columns = [join_key] + [f"element02_{c}" for c in meta_payload_cols]
        
        merged_df = pd.merge(
            merged_df, meta02, left_on='element02', right_on=join_key, how='left'
        ).drop(columns=[join_key])

        # 4. Cleanup and Reorder
        # Fill missing metadata with N/A
        new_merged_cols = [f"element01_{c}" for c in meta_payload_cols] + \
                          [f"element02_{c}" for c in meta_payload_cols]
        merged_df[new_merged_cols] = merged_df[new_merged_cols].fillna('N/A')

        # Interleave columns: [Originals..., element01_ColA, element02_ColA, element01_ColB...]
        interleaved_cols = []
        for c in meta_payload_cols:
            interleaved_cols.extend([f"element01_{c}", f"element02_{c}"])
        
        final_order = data_df.columns.tolist() + interleaved_cols
        merged_df = merged_df[final_order]

        # 5. Save
        merged_df.to_csv(output_file_path, sep='\t', index=False)
        print(f"Success! Processed data saved to: {output_file_path}")

    except Exception as e:
        sys.exit(f"An unexpected error occurred: {e}")

def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Enrich pair-data with metadata for both elements."
    )
    parser.add_argument('-d', '--data-file', required=True, help="Data file (semicolon-sep)")
    parser.add_argument('-m', '--metadata-file', required=True, help="Metadata file")
    parser.add_argument('-o', '--output-file', required=True, help="Output TSV path")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.data_file) or not os.path.exists(args.metadata_file):
        sys.exit("Error: One or more input files do not exist.")

    merge_and_process_files(args.data_file, args.metadata_file, args.output_file)

if __name__ == "__main__":
    main()