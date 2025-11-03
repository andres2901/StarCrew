import sys
import collections
import networkx as nx
import argparse
import os

def find_clusters(pairs):
    """Finds weakly connected components (clusters) from a list of UNWEIGHTED pairs (edges)."""
    graph = collections.defaultdict(set)
    all_nodes = set()
    for u, v in pairs:
        graph[u].add(v)
        graph[v].add(u)
        all_nodes.add(u)
        all_nodes.add(v)

    visited = set()
    clusters = []

    for node in all_nodes:
        if node not in visited:
            current_cluster = set()
            stack = [node]
            
            while stack:
                vertex = stack.pop()
                if vertex not in visited:
                    visited.add(vertex)
                    current_cluster.add(vertex)
                    for neighbor in graph.get(vertex, []):
                        if neighbor not in visited:
                            stack.append(neighbor)
            clusters.append(current_cluster)
            
    return clusters

def write_clusters_to_file(filename, clusters_data):
    """
    Writes clusters to a specified file with proper formatting.
    This version handles the simple list of clusters.
    """
    with open(filename, 'w') as f_out:
        sorted_clusters = sorted(clusters_data, key=len, reverse=True)
        cluster_id_map = {}
        for i, cluster_set in enumerate(sorted_clusters):
            cluster_id = f"Cluster{i+1:04d}"
            cluster_id_map[f"original_{i}"] = cluster_id
            elements = sorted(list(cluster_set))
            f_out.write(f"{cluster_id}\t{' '.join(elements)}\n")
        return cluster_id_map

def write_network_edges(filename, edges):
    """
    Writes the network edges to a file suitable for Cytoscape.
    """
    with open(filename, 'w') as f_out:
        f_out.write("Source\tTarget\tWeight\n")
        for n1, n2, weight in edges:
            f_out.write(f"{n1}\t{n2}\t{weight}\n")
    print(f"Network edges written to '{filename}' for Cytoscape visualization.")

def write_cytoscape_node_attributes(filename, clusters_data, cluster_id_map):
    """
    Writes node attributes (cluster IDs) to a file for Cytoscape.
    """
    node_attributes = {}
    for i, cluster_set in enumerate(clusters_data):
        cluster_id = cluster_id_map[f"original_{i}"]
        for node in cluster_set:
            node_attributes[node] = {'ClusterID': cluster_id}
            
    with open(filename, 'w') as f_out:
        f_out.write("NodeID\tClusterID\n")
        for node_id in sorted(node_attributes.keys()):
            attrs = node_attributes[node_id]
            f_out.write(f"{node_id}\t{attrs['ClusterID']}\n")
    print(f"Node attributes written to '{filename}' for Cytoscape visualization.")

def main():
    """
    Main function to read input, find clusters, and print formatted output.
    """
    # 1. Setup argparse for command-line flags
    parser = argparse.ArgumentParser(description='Simple community detection by finding connected components.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input-file', required=True,
                        help="Path to the input file containing graph edges.")
    parser.add_argument('-o', '--output-dir', default='./',
                        help="The directory to store all output files. Default is the current directory.")

    args = parser.parse_args()

    # 2. Check and create the output directory
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    # 3. Reading the input file with header check and tab delimiter
    all_raw_weighted_edges = []
    unweighted_pairs_for_wccs = []

    try:
        with open(args.input_file, 'r') as f:
            lines = f.readlines()
            start_line = 0

            # Check if the first line is a header by trying to parse the weight column.
            if len(lines) > 0:
                first_line_parts = lines[0].strip().split('\t')
                try:
                    # If parsing the third column as a float fails, it's likely a header.
                    _ = float(first_line_parts[2])
                except (ValueError, IndexError):
                    print("Info: Skipping what appears to be a header row.", file=sys.stderr)
                    start_line = 1

            for line in lines[start_line:]:
                parts = line.strip().split('\t')
                # We expect at least 3 columns for node1, node2, and weight
                if len(parts) >= 3:
                    node1, node2 = parts[0], parts[1]
                    try:
                        weight = float(parts[2])
                        all_raw_weighted_edges.append((node1, node2, weight))
                        unweighted_pairs_for_wccs.append((node1, node2))
                    except ValueError:
                        print(f"Warning: Could not parse weight on line: {line.strip()}. Skipping line.", file=sys.stderr)
                        continue
                else:
                    print(f"Warning: Skipping line due to insufficient columns (expected at least 3, got {len(parts)}): {line.strip()}", file=sys.stderr)
                    continue
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file '{args.input_file}': {e}", file=sys.stderr)
        sys.exit(1)

    # 4. Perform clustering and write output files to the specified directory
    main_clusters = find_clusters(unweighted_pairs_for_wccs)

    main_clusters_output_file = os.path.join(output_dir, "Clusters.txt")
    main_cluster_id_map = write_clusters_to_file(main_clusters_output_file, main_clusters)
    print(f"Main clusters written to '{main_clusters_output_file}'")

    network_edges_output_file = os.path.join(output_dir, "network_edges.txt")
    write_network_edges(network_edges_output_file, all_raw_weighted_edges)

    node_attributes_output_file = os.path.join(output_dir, "node_attributes.txt")
    write_cytoscape_node_attributes(node_attributes_output_file, main_clusters, main_cluster_id_map)
    
if __name__ == "__main__":
    main()

