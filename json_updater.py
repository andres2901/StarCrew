"""
JSON Updater Module.

This script updates OrthoFinder command configurations within a JSON file.
"""

import json
import sys


def load_json(file_path: str) -> dict:
    """Load JSON data safely from a file.

    Args:
        file_path: Path to the JSON file.

    Returns:
        Dictionary containing the loaded JSON data.

    Raises:
        FileNotFoundError: If the file does not exist.
        json.JSONDecodeError: If the file is not valid JSON.
    """

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


def save_json(file_path: str, data: dict) -> None:
    """Save a dictionary back to a JSON file.

    Args:
        file_path: Path to the output JSON file.
        data: Dictionary to serialize and save.

    Raises:
        PermissionError: If the file cannot be written.
        OSError: If any other I/O error occurs.
    """

    try:
        with open(file_path, 'w') as file:
            json.dump(data, file, indent=2)
            print(f"Update complete. The new JSON file is saved at {file_path}.")
    except PermissionError:
        sys.exit(f"Error: No write permission for: '{file_path}'.")
    except OSError as e:
        sys.exit(f"Error writing file '{file_path}': {e}")


def handle_orthofinder_cmd_update(data: dict) -> bool:
    """Update search_cmd for diamond and diamond_ultra_sensitive entries.

    Replaces the search_cmd value for both diamond methods with updated
    command strings. Warns if an expected key is not found in the data.

    Args:
        data: Configuration dictionary loaded from the JSON file.

    Returns:
        True if at least one modification was applied, False otherwise.
    """

    modifications = {
        "diamond": (
            "diamond blastp --threads METHODTHREAD --ignore-warnings -d DATABASE -q INPUT -o OUTPUT "
            "--query-cover 70 --subject-cover 70 --matrix SCOREMATRIX --gapopen GAPOPEN --gapextend GAPEXTEND "
            "--fast --hit-membuf --quiet -e 1e-5 --compress 1"
        ),
        "diamond_ultra_sens": (
            "diamond blastp --threads METHODTHREAD --ignore-warnings -d DATABASE -q INPUT -o OUTPUT "
            "--query-cover 50 --subject-cover 50 --matrix SCOREMATRIX --gapopen GAPOPEN --gapextend GAPEXTEND "
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
            print(f"Warning: Key '{key}' not found in JSON. Skipping.")

    return any_modified


def main() -> None:
    """Parse command-line arguments and launch the JSON update pipeline."""
    if len(sys.argv) < 2:
        print("Usage: python json_updater.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]

    config_data = load_json(file_path)
    modified = False

    print("\n--- Running OrthoFinder Command Update Mode ---")
    modified = handle_orthofinder_cmd_update(config_data)

    if modified:
        save_json(file_path, config_data)
    else:
        print("\nOperation finished, but no changes were applied.")


if __name__ == "__main__":
    main()
