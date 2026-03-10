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


def merge_and_process_files(
    data_file_path: str,
    metadata_file_path: str,
    output_file_path: str
) -> None:
    """Enrich a paired data file with metadata for both element columns.

    Performs two sequential left joins to attach metadata to element01
    and element02 independently, interleaves the resulting columns, fills
    missing values with 'N/A', and writes the result as a TSV file.

    Args:
        data_file_path: Path to the semicolon-separated input data file.
                        Must contain 'element01' and 'element02' columns.
        metadata_file_path: Path to the metadata file. Separator is
                            auto-detected. The second column is used as
                            the join key.
        output_file_path: Path to the output TSV file to write.

    Raises:
        FileNotFoundError: If the data file does not exist.
        PermissionError: If a file cannot be read or written.
        pd.errors.ParserError: If either input file cannot be parsed.
        OSError: If any other I/O error occurs during writing.
    """

    try:
        data_df = pd.read_csv(
            data_file_path, 
            sep=';', 
            dtype={'element01': str, 'element02': str}
        )
        
        try:
            metadata_df = pd.read_csv(metadata_file_path, sep=None, engine='python', dtype=str)
        except pd.errors.ParserError:
            try:
                metadata_df = pd.read_csv(metadata_file_path, sep='\t', dtype=str)
            except pd.errors.ParserError as e:
                sys.exit(f"Error parsing metadata file '{metadata_file_path}': {e}")
        except PermissionError:
            sys.exit(f"Error: No read permission for '{metadata_file_path}'.")

        if len(metadata_df.columns) < 2:
            sys.exit("Error: Metadata file must have at least two columns.")
            
        join_key = "ElementID_updated"
        metadata_df.columns.values[1] = join_key

        meta_payload_cols = metadata_df.columns[2:].tolist()
        
        if not meta_payload_cols:
            print("Warning: No metadata columns found to merge. Saving raw data as TSV.")
            data_df.to_csv(output_file_path, sep='\t', index=False)
            return

        meta01 = metadata_df[[join_key] + meta_payload_cols].copy()
        meta01.columns = [join_key] + [f"element01_{c}" for c in meta_payload_cols]
        
        merged_df = pd.merge(
            data_df, meta01, left_on='element01', right_on=join_key, how='left'
        ).drop(columns=[join_key])
        
        meta02 = metadata_df[[join_key] + meta_payload_cols].copy()
        meta02.columns = [join_key] + [f"element02_{c}" for c in meta_payload_cols]
        
        merged_df = pd.merge(
            merged_df, meta02, left_on='element02', right_on=join_key, how='left'
        ).drop(columns=[join_key])

        new_merged_cols = [f"element01_{c}" for c in meta_payload_cols] + \
                          [f"element02_{c}" for c in meta_payload_cols]
        merged_df[new_merged_cols] = merged_df[new_merged_cols].fillna('N/A')

        interleaved_cols = []
        for c in meta_payload_cols:
            interleaved_cols.extend([f"element01_{c}", f"element02_{c}"])
        
        final_order = data_df.columns.tolist() + interleaved_cols
        merged_df = merged_df[final_order]

        try:
            merged_df.to_csv(output_file_path, sep='\t', index=False)
        except PermissionError:
            sys.exit(f"Error: No write permission for '{output_file_path}'.")
        except OSError as e:
            sys.exit(f"Error writing output file: {e}")

        print(f"Success! Processed data saved to: {output_file_path}")

    except FileNotFoundError:
        sys.exit(f"Error: Data file '{data_file_path}' not found.")
    except PermissionError:
        sys.exit(f"Error: No read permission for '{data_file_path}'.")
    except pd.errors.ParserError as e:
        sys.exit(f"Error parsing data file '{data_file_path}': {e}")


def main() -> None:
    """Parse command-line arguments and launch the metadata enrichment tool."""
    
    parser = argparse.ArgumentParser(
        description="Enrich pair-data with metadata for both elements."
    )
    parser.add_argument('-d', '--data-file', required=True, help="Data file (semicolon-sep)")
    parser.add_argument('-m', '--metadata-file', required=True, help="Metadata file")
    parser.add_argument('-o', '--output-file', required=True, help="Output TSV path")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.data_file):
        sys.exit(f"Error: Data file '{args.data_file}' not found.")
    if not os.path.exists(args.metadata_file):
        sys.exit(f"Error: Metadata file '{args.metadata_file}' not found.")

    merge_and_process_files(args.data_file, args.metadata_file, args.output_file)


if __name__ == "__main__":
    main()