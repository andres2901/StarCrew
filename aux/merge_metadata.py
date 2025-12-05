import pandas as pd
import argparse
import sys

def merge_and_process_files(data_file_path, metadata_file_path, output_file_path):
    """
    Loads two files, performs merges, renames, reorders columns, fills missing 
    values with 'N/A', and saves the result as tab-separated, treating IDs as strings.
    """
    try:
        # --- 1. Load Data and Metadata Files (Ensuring IDs are read as strings) ---
        
        # Data file (Semicolon separated)
        # Ensure element01 and element02 are read as strings
        data_df = pd.read_csv(
            data_file_path, 
            sep=';', 
            dtype={'element01': str, 'element02': str}
        )
        
        # Metadata file 
        # Read all columns as string to treat ElementID_updated correctly
        try:
            metadata_df = pd.read_csv(metadata_file_path, sep=None, engine='python', on_bad_lines='skip', dtype=str)
        except:
            metadata_df = pd.read_csv(metadata_file_path, sep='\t', dtype=str)

        # --- 2. Prepare Metadata Columns for Merging ---
        
        if len(metadata_df.columns) < 2:
            print("Error: Metadata file must have at least two columns.")
            sys.exit(1)
            
        # Rename the second column to the specific join key, treating it as string
        metadata_df.columns.values[1] = "ElementID_updated"

        # The data that needs to be added starts in the third column (index 2)
        metadata_cols_to_merge = metadata_df.columns[2:].tolist()
        
        if not metadata_cols_to_merge:
            print("Warning: Metadata file contains only key columns (first two). Nothing to merge.")
            # Still save the original data file as output
            data_df.to_csv(output_file_path, sep='\t', index=False)
            print(f"Saved original data file (tab-separated) to: {output_file_path}")
            return
            
        # Create two sets of renamed columns
        metadata_cols_element01 = {col: f"element01_{col}" for col in metadata_cols_to_merge}
        metadata_cols_element02 = {col: f"element02_{col}" for col in metadata_cols_to_merge}

        # List of the final new columns we expect to create (for reordering later)
        new_cols_01 = list(metadata_cols_element01.values())
        new_cols_02 = list(metadata_cols_element02.values())
        
        # --- 3. Perform the Merges ---
        
        metadata_for_join = metadata_df[['ElementID_updated'] + metadata_cols_to_merge].copy()

        # 3a. Merge 1: Based on 'element01' and 'ElementID_updated'
        metadata_01 = metadata_for_join.rename(columns=metadata_cols_element01)
        
        merged_df = pd.merge(
            data_df, 
            metadata_01, 
            left_on='element01', 
            right_on='ElementID_updated', 
            how='left'
        ).drop(columns=['ElementID_updated'])
        
        
        # 3b. Merge 2: Based on 'element02' and 'ElementID_updated'
        metadata_02 = metadata_for_join.rename(columns=metadata_cols_element02)

        merged_df = pd.merge(
            merged_df, 
            metadata_02, 
            left_on='element02', 
            right_on='ElementID_updated', 
            how='left'
        ).drop(columns=['ElementID_updated'])


        # --- 4. Fill Missing Values with 'N/A' ---
        
        # Identify all the newly merged columns
        merged_metadata_cols = new_cols_01 + new_cols_02
        
        # Replace NaN in only the newly merged columns with 'N/A'
        # Note: If any original data file columns also have NaN and need to be filled, 
        # you would need to adjust the logic here. This targets only the merged columns.
        merged_df[merged_metadata_cols] = merged_df[merged_metadata_cols].fillna('N/A')

        
        # --- 5. Reorder Columns for Final Output ---
        
        original_cols = data_df.columns.tolist()
        
        # Create the desired final column order (element01_X, element02_X, element01_Y, element02_Y, ...)
        ordered_new_cols = []
        for col_name in metadata_cols_to_merge:
            ordered_new_cols.append(f"element01_{col_name}")
            ordered_new_cols.append(f"element02_{col_name}")

        final_column_order = original_cols + ordered_new_cols
        
        # Apply the new column order
        merged_df = merged_df[final_column_order]

        # --- 6. Save the Output File ---
        
        merged_df.to_csv(output_file_path, sep='\t', index=False)
        
        print(f"Processing complete. Merged, cleaned, and saved to: {output_file_path}")

    except FileNotFoundError:
        print("Error: One or both input files were not found. Please check the paths.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print("Error: One or both input files are empty.")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge a data file with a metadata file, fill missing values with 'N/A', and output tab-separated data."
    )
    
    parser.add_argument(
        '-d', '--data-file', 
        required=True, 
        help="Path to the primary data file (semicolon-separated)."
    )
    parser.add_argument(
        '-m', '--metadata-file', 
        required=True, 
        help="Path to the metadata file (with ElementID_updated as the second column)."
    )
    parser.add_argument(
        '-o', '--output-file', 
        required=True, 
        help="Path for the output file (tab-separated)."
    )
    
    args = parser.parse_args()
    
    merge_and_process_files(args.data_file, args.metadata_file, args.output_file)