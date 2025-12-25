import json
import sys

def load_json(file_path):
    """Loads JSON data safely."""
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
        print(f"JSON file '{file_path}' loaded successfully.")
        return data
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: The file '{file_path}' is not a valid JSON format.")
        sys.exit(1)

def save_json(file_path, data):
    """Saves modified data back to the JSON file."""
    try:
        with open(file_path, 'w') as file:
            json.dump(data, file, indent=4)
        print(f"Update complete. The new JSON file is saved at {file_path}.")
    except Exception as e:
        print(f"An error occurred while writing to the file: {e}")
        sys.exit(1)

# --- MODIFICATION HANDLERS ---

def handle_orthofinder_cmd_update(data):
    """Handles the diamond search_cmd modification (Example 1)."""
    DIAMOND_KEY = 'diamond'
    SEARCH_CMD_KEY = 'search_cmd'
    
    NEW_SEARCH_CMD = (
        "diamond blastp --threads METHODTHREAD --ignore-warnings -d DATABASE -q INPUT -o OUTPUT --query-cover 80 --matrix SCOREMATRIX --gapopen GAPOPEN --gapextend GAPEXTEND --fast -p 1 --hit-membuf --quiet -e 1e-5 --compress 1"
    )

    if DIAMOND_KEY in data and SEARCH_CMD_KEY in data[DIAMOND_KEY]:
        data[DIAMOND_KEY][SEARCH_CMD_KEY] = NEW_SEARCH_CMD
        print("\nModified 'diamond' -> 'search_cmd' value.")
        return True
    else:
        print("\nWarning: Could not find 'diamond' or 'search_cmd' for update. Skipping.")
        return False

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    
    # Simple argument check for required arguments
    if len(sys.argv) < 3:
        print("Usage:")
        print("OrthoFinder Command Update: python json_updater.py <file_path> --cmd-update")
        sys.exit(1)

    file_path = sys.argv[1]
    mode = sys.argv[2]
    
    # 1. Load the JSON data
    config_data = load_json(file_path)
    
    modified = False
    
    # 2. Select and run the correct modification handler based on the mode
    if mode == '--cmd-update':
        print("\n--- Running OrthoFinder Command Update Mode ---")
        modified = handle_orthofinder_cmd_update(config_data)
        
    else:
        print(f"Error: Unknown mode '{mode}'. Use '--cmd-update'.")
        sys.exit(1)

    # 3. Save the JSON data if modifications were made
    if modified:
        save_json(file_path, config_data)
    else:
        print("\nOperation finished, but no relevant changes were applied.")
