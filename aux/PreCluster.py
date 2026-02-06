#!/usr/bin/env python3
"""
Network Clustering Tool.

Finds Weakly Connected Components (WCCs) in an undirected graph.
"""

import sys
import collections
import argparse
import os

def read_input_file(input_path):
    """Parses a two-column file into an unweighted edge list."""
    print(f"STEP 1: Reading input from '{input_path}'...")
    edges = []

    try:
        with open(input_path, 'r') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                # Support both Tab and Space separators
                parts = line.split('\t') if '\t' in line else line.split()
                
                if len(parts) >= 2:
                    edges.append((parts[0], parts[1]))
                else:
                    print(f"Warning: Skipping line {i+1} (needs 2 columns).", file=sys.stderr)
    except FileNotFoundError:
        sys.exit(f"Error: File '{input_path}' not found.")
    except Exception as e:
        sys.exit(f"Error reading file: {e}")

    print(f"Finished. Total edges: {len(edges)}")
    return edges

def find_clusters(pairs):
    """
    Groups nodes into Weakly Connected Components using iterative DFS.
    """
    # 1. Build Adjacency List
    graph = collections.defaultdict(set)
    all_nodes = set()
    for u, v in pairs:
        graph[u].add(v)
        graph[v].add(u)
        all_nodes.add(u)
        all_nodes.add(v)

    visited = set()
    clusters = []

    # 2. DFS Traversal
    for node in all_nodes:
        if node not in visited:
            current_cluster = set()
            stack = [node]
            
            while stack:
                vertex = stack.pop()
                if vertex not in visited:
                    visited.add(vertex)
                    current_cluster.add(vertex)
                    # Add unvisited neighbors to stack
                    for neighbor in graph.get(vertex, []):
                        if neighbor not in visited:
                            stack.append(neighbor)
            
            clusters.append(current_cluster)
            
    return clusters

def write_results(output_dir, clusters):
    """Writes clusters to a file, sorted by size (largest first)."""
    output_path = os.path.join(output_dir, "Clusters.txt")
    print(f"STEP 3: Writing clusters to {output_path}...")

    # Sort: Largest components first, then alphabetically by first node
    sorted_clusters = sorted(clusters, key=lambda x: (-len(x), sorted(list(x))[0]))

    with open(output_path, 'w') as f_out:
        for i, cluster_set in enumerate(sorted_clusters):
            cluster_id = f"Cluster{i+1:04d}"
            elements = sorted(list(cluster_set))
            f_out.write(f"{cluster_id}\t{' '.join(elements)}\n")
            
    print(f"Successfully identified {len(clusters)} clusters.")

def main():
    """CLI Entry Point."""
    parser = argparse.ArgumentParser(
        description='Community detection via Connected Components.'
    )
    parser.add_argument('-i', '--input-file', required=True, help="Input (NodeA NodeB)")
    parser.add_argument('-o', '--output-dir', default='./', help="Output directory")

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Execution logic
    edge_list = read_input_file(args.input_file)
    if not edge_list:
        sys.exit("No data to process.")

    print("STEP 2: Finding Connected Components...")
    clusters = find_clusters(edge_list)
    
    write_results(args.output_dir, clusters)

if __name__ == "__main__":
    main()