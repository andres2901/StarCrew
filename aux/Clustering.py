import sys
import collections
import networkx as nx
import numpy as np
from sklearn.cluster import SpectralClustering
import argparse
import os
import re # Added for regular expression parsing

def find_clusters(pairs):
    """Finds weakly connected components (clusters) from a list of WEIGHTED pairs (edges)."""
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

def find_sub_clusters_spectral(all_raw_weighted_edges, main_clusters, 
                                 min_final_cluster_size=2,
                                 min_nodes_for_meaningful_spectral_split=4, 
                                 modularity_split_threshold=0.05):
    """Applies an iterative spectral clustering approach."""
    sub_clusters_with_parent_info = []
    
    G_full = nx.Graph() 
    for n1, n2, weight in all_raw_weighted_edges:
        if weight > 0:
            G_full.add_edge(n1, n2, weight=weight)

    for original_cluster_idx, main_cluster_nodes in enumerate(main_clusters):
        subgraph = G_full.subgraph(main_cluster_nodes)

        if (subgraph.number_of_nodes() < min_nodes_for_meaningful_spectral_split or 
            subgraph.number_of_edges() == 0):
            sub_clusters_with_parent_info.append((original_cluster_idx, main_cluster_nodes, False))
            continue

        adj_matrix = nx.to_numpy_array(subgraph, weight='weight')
        node_list = list(subgraph.nodes())
        
        best_partition = None
        best_modularity_score = -1.0
        
        k = 2
        while True:
            if k > len(node_list):
                break

            try:
                sc = SpectralClustering(n_clusters=k, affinity='precomputed', n_init=10, assign_labels='kmeans', random_state=0)
                labels = sc.fit_predict(adj_matrix)
                
                proposed_partition = collections.defaultdict(set)
                for i, node_label in enumerate(labels):
                    proposed_partition[node_label].add(node_list[i])
                
                proposed_partition_list = list(proposed_partition.values())

                modularity_score = nx.community.modularity(subgraph, proposed_partition_list, weight='weight')
                
                if modularity_score > best_modularity_score and modularity_score > modularity_split_threshold:
                    best_modularity_score = modularity_score
                    best_partition = proposed_partition_list
                    k += 1
                else:
                    break
            except Exception as e:
                print(f"Warning: Spectral clustering failed for cluster {original_cluster_idx} with k={k}. Error: {e}", file=sys.stderr)
                break
                
        if best_partition:
            current_communities = collections.defaultdict(set)
            for i, community_set in enumerate(best_partition):
                current_communities[i].update(community_set)
            
            something_merged_in_iteration = True
            while something_merged_in_iteration:
                something_merged_in_iteration = False
                small_cluster_ids = [cid for cid, c_nodes in current_communities.items() if len(c_nodes) < min_final_cluster_size]
                if not small_cluster_ids or (len(small_cluster_ids) == 1 and len(current_communities) == 1):
                    break

                small_cluster_ids.sort(key=lambda cid: len(current_communities[cid]))
                cluster_to_merge_id = small_cluster_ids[0]
                cluster_to_merge_nodes = current_communities[cluster_to_merge_id]
                best_merge_target_id = None
                max_aggregate_weight = -1 
                for target_cluster_id, target_cluster_nodes in current_communities.items():
                    if target_cluster_id == cluster_to_merge_id:
                        continue 
                    current_aggregate_weight = 0
                    for node_in_merge_cluster in cluster_to_merge_nodes:
                        for neighbor in subgraph.neighbors(node_in_merge_cluster):
                            if neighbor in target_cluster_nodes: 
                                current_aggregate_weight += subgraph[node_in_merge_cluster][neighbor]['weight']
                    if best_merge_target_id is None or current_aggregate_weight > max_aggregate_weight:
                        max_aggregate_weight = current_aggregate_weight
                        best_merge_target_id = target_cluster_id
                
                if best_merge_target_id is not None and max_aggregate_weight > 0:
                    current_communities[best_merge_target_id].update(cluster_to_merge_nodes)
                    del current_communities[cluster_to_merge_id] 
                    something_merged_in_iteration = True 
            
            final_sub_cluster_map = collections.defaultdict(set)
            for cid, c_nodes in current_communities.items():
                if c_nodes: 
                    final_sub_cluster_map[cid].update(c_nodes)

            if len(final_sub_cluster_map) > 1:
                for community_set in final_sub_cluster_map.values():
                    if community_set: 
                        sub_clusters_with_parent_info.append((original_cluster_idx, community_set, True))
            else:
                sub_clusters_with_parent_info.append((original_cluster_idx, main_cluster_nodes, False))
        else:
            sub_clusters_with_parent_info.append((original_cluster_idx, main_cluster_nodes, False))

    return sub_clusters_with_parent_info

def write_clusters_to_file(filename, clusters_data, main_cluster_id_map=None):
    """Writes clusters to a specified file with proper formatting."""
    with open(filename, 'w') as f_out:
        if main_cluster_id_map:
            processed_sub_clusters = []
            # Sorting sub-clusters first by parent index, then by the smallest element in the set for stable IDs
            sorted_raw_sub_clusters = sorted(clusters_data, key=lambda item: (item[0], sorted(list(item[1]))[0] if item[1] else ''))
            temp_final_parent_id_sub_cluster_counts = collections.defaultdict(int)
            for original_parent_idx, cluster_set, was_split in sorted_raw_sub_clusters:
                parent_cluster_id = main_cluster_id_map[original_parent_idx]
                if was_split:
                    temp_final_parent_id_sub_cluster_counts[parent_cluster_id] += 1
                    sub_id_num = temp_final_parent_id_sub_cluster_counts[parent_cluster_id]
                    sub_id_str = f"{sub_id_num:02d}"
                    full_cluster_id = f"{parent_cluster_id}.{sub_id_str}"
                else:
                    sub_id_num = 0 
                    full_cluster_id = parent_cluster_id 
                processed_sub_clusters.append((full_cluster_id, sub_id_num, cluster_set))
            
            # Final sorting by main cluster ID and then sub cluster number
            sorted_final_sub_clusters = sorted(processed_sub_clusters, key=lambda item: (item[0].split('.')[0], item[1]))
            for cluster_id, _, cluster_set in sorted_final_sub_clusters:
                elements = sorted(list(cluster_set))
                f_out.write(f"{cluster_id}\t{' '.join(elements)}\n")
        else:
            sorted_main_clusters = sorted(clusters_data, key=len, reverse=True)
            current_main_cluster_id_map = {} 
            for i, cluster_set in enumerate(sorted_main_clusters):
                original_idx = -1
                # Find the original index of the current cluster set (required for mapping in sub-clustering)
                for idx, s in enumerate(clusters_data):
                    # NOTE: comparing sets is slow, but necessary here for mapping back to the original index
                    if s == cluster_set: 
                        original_idx = idx
                        break
                if original_idx == -1: 
                    raise ValueError(f"Cluster set {cluster_set} not found in original main_clusters list. "
                                     "This indicates an unexpected issue in matching clusters.")
                cluster_id = f"Cluster{i+1:04d}" 
                current_main_cluster_id_map[original_idx] = cluster_id
                elements = sorted(list(cluster_set))
                f_out.write(f"{cluster_id}\t{' '.join(elements)}\n")
            return current_main_cluster_id_map

def write_network_edges(filename, edges):
    """Writes the network edges to a file suitable for Cytoscape."""
    with open(filename, 'w') as f_out:
        f_out.write("Source\tTarget\tWeight\n")
        for n1, n2, weight in edges:
            f_out.write(f"{n1}\t{n2}\t{weight}\n")
    print(f"Network edges written to '{filename}' for Cytoscape visualization.")

def write_cytoscape_node_attributes(filename, sub_clusters_info, main_cluster_id_map, node_metadata, metadata_keys):
    """Writes node attributes (cluster IDs, split status, and dynamic metadata) to a file for Cytoscape."""
    node_attributes = {}
    
    # 1. Map sub-cluster IDs to nodes
    sorted_raw_sub_clusters = sorted(sub_clusters_info, key=lambda item: (item[0], sorted(list(item[1]))[0] if item[1] else ''))
    temp_final_parent_id_sub_cluster_counts = collections.defaultdict(int)
    for original_parent_idx, cluster_set, was_split in sorted_raw_sub_clusters:
        parent_cluster_id = main_cluster_id_map[original_parent_idx]
        if was_split:
            temp_final_parent_id_sub_cluster_counts[parent_cluster_id] += 1
            sub_id_num = temp_final_parent_id_sub_cluster_counts[parent_cluster_id]
            sub_id_str = f"{sub_id_num:02d}"
            full_cluster_id = f"{parent_cluster_id}.{sub_id_str}"
        else:
            full_cluster_id = parent_cluster_id
            
        for node in cluster_set:
            if node not in node_attributes:
                node_attributes[node] = {
                    'SubClusterID': full_cluster_id,
                    'MainClusterID': parent_cluster_id,
                    'WasSplit': was_split
                }

    # 2. Determine which headers to use dynamically
    default_headers = ["NodeID", "SubClusterID", "MainClusterID", "WasSplit"]
    
    # Check for derived keys like 'genus' if 'species' was present in the metadata
    final_metadata_keys = []
    if 'species' in metadata_keys:
        final_metadata_keys.append('genus') # Genus is derived from species
    final_metadata_keys.extend(metadata_keys)
        
    header = '\t'.join(default_headers + final_metadata_keys) + '\n'

    # 3. Write data
    with open(filename, 'w') as f_out:
        f_out.write(header)
        for node_id in sorted(node_attributes.keys()):
            attrs = node_attributes[node_id]
            line_parts = [
                node_id, 
                attrs['SubClusterID'], 
                attrs['MainClusterID'], 
                str(attrs['WasSplit'])
            ]
            
            # Add dynamic metadata
            metadata = node_metadata.get(node_id, {})
            for key in final_metadata_keys:
                # Provide 'N/A' fallback for missing data
                line_parts.append(str(metadata.get(key, 'N/A')))
            
            f_out.write('\t'.join(line_parts) + '\n')
            
    print(f"Node attributes written to '{filename}' for Cytoscape visualization.")

def write_cluster_stats(filename, main_clusters, sub_clusters_with_parent_info, main_cluster_id_map):
    """Generates a statistics file for each main cluster."""
    cluster_stats = collections.defaultdict(lambda: {'Size': 0, 'NumberSubCluster': 1})
    
    # Populate initial stats from main clusters
    for i, cluster_set in enumerate(main_clusters):
        original_idx = -1
        # Find the original index of the current cluster set
        for idx, s in enumerate(main_clusters):
            if s == cluster_set:
                original_idx = idx
                break
        if original_idx != -1:
            cluster_id = main_cluster_id_map[original_idx]
            cluster_stats[cluster_id]['Size'] = len(cluster_set)

    # Update stats based on sub-cluster splits
    sub_cluster_counts = collections.defaultdict(int)
    for original_parent_idx, _, was_split in sub_clusters_with_parent_info:
        if was_split:
            parent_cluster_id = main_cluster_id_map[original_parent_idx]
            sub_cluster_counts[parent_cluster_id] += 1
    
    # Set the final sub-cluster count
    for cluster_id, count in sub_cluster_counts.items():
        cluster_stats[cluster_id]['NumberSubCluster'] = count

    # Write the statistics to the file
    with open(filename, 'w') as f_out:
        f_out.write("ClusterID\tSize\tNumberSubCluster\n")
        sorted_cluster_ids = sorted(cluster_stats.keys())
        for cluster_id in sorted_cluster_ids:
            stats = cluster_stats[cluster_id]
            f_out.write(f"{cluster_id}\t{stats['Size']}\t{stats['NumberSubCluster']}\n")
            
    print(f"Cluster statistics written to '{filename}'.")

def process_clustering(input_file, output_dir, min_final_cluster_size, min_nodes_for_split, modularity_split_threshold):
    """function to read input, find clusters, sort them, and print formatted output."""

    # 1. Initialize data structures
    all_raw_weighted_edges = []
    unweighted_pairs_for_wccs = []
    node_metadata = {}
    metadata_keys = []  # Stores the attribute names (e.g., 'species', 'length')
    metadata_col_map = {} # Stores {'key': {'node1_idx': i, 'node2_idx': i+1}}

    # 2. Check and create the output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
        
    # 3. Reading the input file and storing metadata with automatic column detection
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
            start_line = 0

            # --- Header detection and parsing (NEW LOGIC) ---
            if len(lines) > 0:
                first_line_parts = lines[0].strip().split('\t')
                
                # Assume a header if the third column is not parsable as a float (weight)
                try:
                    _ = float(first_line_parts[2])
                except (ValueError, IndexError):
                    print("Info: Skipping what appears to be a header row.", file=sys.stderr)
                    start_line = 1
                    
                    # Dynamic Header Parsing
                    header_parts = lines[0].strip().split('\t')
                    
                    # Metadata starts at index 5 (after 5th column, which are indices 0-4)
                    if len(header_parts) > 5:
                        for i in range(5, len(header_parts), 2):
                            col1_header = header_parts[i]
                            col2_header = header_parts[i+1] if i + 1 < len(header_parts) else None
                            
                            # Check for the required paired pattern: element01_KEY and element02_KEY
                            match1 = re.match(r"element01_(.*)", col1_header)
                            match2 = re.match(r"element02_(.*)", col2_header) if col2_header else None
                            
                            if match1 and match2 and match1.group(1) == match2.group(1):
                                key = match1.group(1)
                                metadata_keys.append(key)
                                metadata_col_map[key] = {'node1_idx': i, 'node2_idx': i + 1}
                            else:
                                # Stop metadata parsing if the paired pattern is broken
                                if len(header_parts) > i:
                                    print(f"Warning: Metadata column pattern broken at index {i}. Expected 'element01_KEY'/'element02_KEY' pair. Stopping metadata parsing.", file=sys.stderr)
                                break
                    
            has_metadata = bool(metadata_keys)
            
            # --- Data line reading and processing ---
            for line in lines[start_line:]:
                parts = line.strip().split('\t')
                
                # Check for minimum required columns (Node1, Node2, Weight)
                if len(parts) < 3:
                    print(f"Warning: Skipping line due to insufficient columns (expected at least 3, got {len(parts)}): {line.strip()}", file=sys.stderr)
                    continue

                node1, node2 = parts[0], parts[1]
                try:
                    weight = float(parts[2])
                    all_raw_weighted_edges.append((node1, node2, weight))
                    unweighted_pairs_for_wccs.append((node1, node2))

                    # Dynamic Metadata Extraction
                    if has_metadata:
                        # Initialize temporary storage for this line's metadata
                        line_metadata_node1 = {}
                        line_metadata_node2 = {}
                        
                        for key in metadata_keys:
                            indices = metadata_col_map[key]
                            
                            # Check if data line has enough columns for this metadata pair
                            if indices['node2_idx'] < len(parts):
                                val1 = parts[indices['node1_idx']]
                                val2 = parts[indices['node2_idx']]
                                
                                # Special handling for 'genus' if 'species' is a key (derived field)
                                if key.lower() == 'species':
                                    line_metadata_node1['genus'] = val1.split()[0] if val1 and ' ' in val1 else val1
                                    line_metadata_node2['genus'] = val2.split()[0] if val2 and ' ' in val2 else val2
                                
                                # Store the primary metadata key
                                line_metadata_node1[key] = val1
                                line_metadata_node2[key] = val2
                            else:
                                print(f"Warning: Data line is shorter than expected for metadata key '{key}'. Using 'N/A'.", file=sys.stderr)
                                
                        # Store gathered metadata for both nodes (only if not already stored)
                        if line_metadata_node1 and node1 not in node_metadata:
                            node_metadata[node1] = line_metadata_node1
                        
                        if line_metadata_node2 and node2 not in node_metadata:
                            node_metadata[node2] = line_metadata_node2
                            
                except ValueError:
                    print(f"Warning: Could not parse weight on line: {line.strip()}. Skipping line.", file=sys.stderr)
                    continue

    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file '{input_file}': {e}", file=sys.stderr)
        sys.exit(1)

    # 4. Perform clustering and write output files to the specified directory
    main_clusters = find_clusters(unweighted_pairs_for_wccs)

    sub_clusters_with_parent_info = find_sub_clusters_spectral(
        all_raw_weighted_edges,
        main_clusters,
        min_final_cluster_size=min_final_cluster_size,
        min_nodes_for_meaningful_spectral_split=min_nodes_for_split,
        modularity_split_threshold=modularity_split_threshold
    )

    main_clusters_output_file = os.path.join(output_dir, "main_clusters.txt")
    main_cluster_id_map = write_clusters_to_file(main_clusters_output_file, main_clusters)
    print(f"Main clusters written to '{main_clusters_output_file}'")

    sub_clusters_output_file = os.path.join(output_dir, "sub_clusters.txt")
    write_clusters_to_file(sub_clusters_output_file, sub_clusters_with_parent_info, main_cluster_id_map=main_cluster_id_map)
    print(f"Sub-clusters written to '{sub_clusters_output_file}'")

    network_edges_output_file = os.path.join(output_dir, "network_edges.txt")
    write_network_edges(network_edges_output_file, all_raw_weighted_edges)

    node_attributes_output_file = os.path.join(output_dir, "node_attributes.txt")
    # Pass the dynamically determined metadata_keys list
    write_cytoscape_node_attributes(node_attributes_output_file, sub_clusters_with_parent_info, main_cluster_id_map, node_metadata, metadata_keys)
    
    cluster_stats_output_file = os.path.join(output_dir, "cluster_stats.txt")
    write_cluster_stats(cluster_stats_output_file, main_clusters, sub_clusters_with_parent_info, main_cluster_id_map)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Community detection on a weighted graph with iterative spectral clustering. Supports dynamic metadata columns.')
    parser.add_argument('-i', '--input-file', required=True,
                        help="Path to the input file containing graph edges and optional metadata.")
    parser.add_argument('-o', '--output-dir', default='./',
                        help="The directory to store all output files. Default is the current directory.")
    parser.add_argument('-m', '--min-size', type=int, default=1, dest='min_final_cluster_size',
                        help="The minimum desired size for any final sub-cluster. (default: 1).")
    parser.add_argument('-n', '--min-nodes', type=int, default=4, dest='min_nodes_for_split',
                        help="Minimum number of nodes in a cluster for spectral clustering to be attempted (default: 4).")
    parser.add_argument('-t', '--threshold', type=float, default=0.05, dest='modularity_split_threshold',
                        help="The minimum modularity score for a split to be accepted (default: 0.05) [range: -0.5, 1.0].")

    args = parser.parse_args()
    process_clustering(args.input_file, args.output_dir, args.min_final_cluster_size, args.min_nodes_for_split, args.modularity_split_threshold)