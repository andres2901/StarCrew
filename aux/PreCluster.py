import sys
import collections
import argparse
import os

def read_input_file(input_path):
    """
    Reads the two-column input file (Node A, Node B). 
    
    Args:
        input_path (str): Path to the input file.
        
    Returns:
        tuple: (unweighted_pairs, weighted_edges)
               unweighted_pairs: list of (node1, node2) for clustering.
               weighted_edges: list of (node1, node2, 1.0) for network output.
    """
    print(f"STEP 1: Reading input file from '{input_path}'...")
    all_raw_weighted_edges = []
    unweighted_pairs_for_wccs = []
    default_weight = 1.0 # Default weight for Cytoscape output

    try:
        with open(input_path, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith('#'):
                    # Skip empty lines or lines starting with a comment hash
                    continue

                # Try splitting by tab, then by space
                parts = line.split('\t')
                if len(parts) < 2:
                    parts = line.split() 
                
                if len(parts) >= 2:
                    node1, node2 = parts[0], parts[1]
                    
                    # Store for network output (with default weight)
                    all_raw_weighted_edges.append((node1, node2, default_weight))
                    
                    # Store for WCC calculation
                    unweighted_pairs_for_wccs.append((node1, node2))
                else:
                    print(f"Warning: Skipping line {i+1} due to insufficient columns (expected 2): {line}", file=sys.stderr)
                    continue
    except FileNotFoundError:
        print(f"Error: Input file '{input_path}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file '{input_path}': {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Finished reading. Found {len(unweighted_pairs_for_wccs)} edges.")
    return unweighted_pairs_for_wccs, all_raw_weighted_edges

def find_clusters(pairs):
    """
    Finds weakly connected components (clusters) from a list of UNWEIGHTED pairs (edges) 
    using an iterative Depth First Search (DFS).
    
    Args:
        pairs (list of tuples): List of (node1, node2) edges.
        
    Returns:
        list of sets: A list where each set represents a cluster (WCC).
    """
    # 1. Build the adjacency list and collect all unique nodes
    graph = collections.defaultdict(set)
    all_nodes = set()
    for u, v in pairs:
        # Add edges in both directions since this is an undirected graph (WCCs)
        graph[u].add(v)
        graph[v].add(u)
        all_nodes.add(u)
        all_nodes.add(v)

    visited = set()
    clusters = []

    # 2. Iterate through all nodes to find unvisited ones (start of a new cluster)
    for node in all_nodes:
        if node not in visited:
            current_cluster = set()
            stack = [node]  # Use stack for iterative DFS
            
            while stack:
                vertex = stack.pop()
                if vertex not in visited:
                    visited.add(vertex)
                    current_cluster.add(vertex)
                    
                    # Push neighbors onto the stack
                    for neighbor in graph.get(vertex, []):
                        if neighbor not in visited:
                            stack.append(neighbor)
                            
            clusters.append(current_cluster)
            
    return clusters

def write_clusters_to_file(filename, clusters_data):
    """
    Writes clusters to a specified file.
    Clusters are sorted by size (largest first) and assigned a sequential ID.
    
    Args:
        filename (str): Path to the output file.
        clusters_data (list of sets): List of sets, each being a cluster.
        
    Returns:
        dict: A map of original cluster index to the formatted Cluster ID string (based on sorted order).
    """
    print(f"Writing clusters to {filename}...")
    with open(filename, 'w') as f_out:
        # Sort clusters by size in descending order
        sorted_clusters = sorted(clusters_data, key=len, reverse=True)
        cluster_id_map = {}
        
        for i, cluster_set in enumerate(sorted_clusters):
            # Format cluster ID: e.g., Cluster0001, Cluster0002
            cluster_id = f"Cluster{i+1:04d}"
            
            # Map the index in the *sorted* list to the formatted ID
            cluster_id_map[f"original_{i}"] = cluster_id 
            
            elements = sorted(list(cluster_set))
            f_out.write(f"{cluster_id}\t{' '.join(elements)}\n")
            
        return cluster_id_map

def process_clustering(input_path, output_dir):
    """
    Orchestrates the entire clustering process: reads data, finds clusters, 
    and writes all output files.
    """
    
    # 1. Read input data
    unweighted_pairs_for_wccs, all_raw_weighted_edges = read_input_file(input_path)
    
    if not unweighted_pairs_for_wccs:
        print("Aborting: No valid edges found in the input file.", file=sys.stderr)
        return

    print("STEP 2: Finding Connected Components...")
    # 2. Perform clustering
    main_clusters = find_clusters(unweighted_pairs_for_wccs)
    print(f"Found {len(main_clusters)} connected components (clusters).")

    # 3. Write output files
    print("STEP 3: Writing output files...")

    # Write clusters file
    main_clusters_output_file = os.path.join(output_dir, "Clusters.txt")
    main_cluster_id_map = write_clusters_to_file(main_clusters_output_file, main_clusters)
    
    print("\nProcessing complete.")


def main():
    """Sets up command-line arguments and calls the main processing function."""
    parser = argparse.ArgumentParser(description='Simple community detection by finding connected components from a two-column (unweighted) input file.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input-file', required=True,
                        help="Path to the input file containing two columns (Node A, Node B).")
    parser.add_argument('-o', '--output-dir', default='./',
                        help="The directory to store all output files. Default is the current directory.")

    args = parser.parse_args()

    # Check and create the output directory
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
        
    process_clustering(args.input_file, args.output_dir)

if __name__ == "__main__":
    main()
