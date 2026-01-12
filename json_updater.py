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
    """Updates search_cmd for diamond and diamond_ultra_sensitive with distinct values."""
    
    # Define a dictionary where each key has its own unique command string
    # I have updated 'diamond_ultra_sensitive' to include the --ultra-sensitive flag as an example
    modifications = {
        "diamond": (
            "diamond blastp --threads METHODTHREAD --ignore-warnings -d DATABASE -q INPUT -o OUTPUT "
            "--query-cover 80 --matrix SCOREMATRIX --gapopen GAPOPEN --gapextend GAPEXTEND "
            "--fast --hit-membuf --quiet -e 1e-5 --compress 1"
        ),
        "diamond_ultra_sens": (
            "diamond blastp --threads METHODTHREAD --ignore-warnings -d DATABASE -q INPUT -o OUTPUT "
            "--query-cover 50 --matrix SCOREMATRIX --gapopen GAPOPEN --gapextend GAPEXTEND "
            "--ultra-sensitive --hit-membuf --quiet -e 0.001 --compress 1"
        )
    }

    any_modified = False

    for key, new_command in modifications.items():
        if key in data and 'search_cmd' in data[key]:
            data[key]['search_cmd'] = new_command
            print(f"Modified '{key}' -> 'search_cmd' with unique parameters.")
            any_modified = True
        else:
            print(f"Warning: Key '{key}' not found. Skipping.")

    return any_modified

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    
    if len(sys.argv) < 3:
        print("Usage: python json_updater.py <file_path> --cmd-update")
        sys.exit(1)

    file_path = sys.argv[1]
    mode = sys.argv[2]
    
    config_data = load_json(file_path)
    modified = False
    
    if mode == '--cmd-update':
        print("\n--- Running OrthoFinder Command Update Mode ---")
        modified = handle_orthofinder_cmd_update(config_data)
    else:
        print(f"❌ Error: Unknown mode '{mode}'. Use '--cmd-update'.")
        sys.exit(1)

    if modified:
        save_json(file_path, config_data)
    else:
        print("\nOperation finished, but no changes were applied.")
