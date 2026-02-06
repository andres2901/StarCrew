#!/usr/bin/env python3
"""
Network Clustering and Spectral Analysis Script.

This script identifies weakly connected components in a network and
applies adaptive spectral clustering to find sub-clusters based on
modularity and connectivity.
"""

import sys
import collections
import os
import re
import argparse
import networkx as nx
import numpy as np
from sklearn.cluster import SpectralClustering


def find_clusters(pairs):
    """
    Find weakly connected components (clusters) from a list of edges.

    Uses a stack-based Depth-First Search (DFS) to identify all nodes
    reachable from an unvisited starting point.
    """
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
                               min_final_cluster_size,
                               min_nodes_for_meaningful_spectral_split,
                               modularity_split_threshold):
    """
    Apply modularity-based analysis followed by adaptive spectral clustering.

    Ensures all final subclusters are formally connected components and
    merges clusters smaller than the specified minimum size.
    """
    sub_clusters_with_parent_info = []
    g_full = nx.Graph()
    for n1, n2, weight in all_raw_weighted_edges:
        if weight > 0:
            g_full.add_edge(n1, n2, weight=weight)

    for original_cluster_idx, main_cluster_nodes in enumerate(main_clusters):
        subgraph = g_full.subgraph(main_cluster_nodes)

        # Skip if too small for splitting
        if (subgraph.number_of_nodes() < min_nodes_for_meaningful_spectral_split or
                subgraph.number_of_edges() == 0):
            sub_clusters_with_parent_info.append((original_cluster_idx, main_cluster_nodes, False))
            continue

        # Step 1: Modularity-based Pre-analysis
        mod_communities = nx.community.louvain_communities(subgraph, weight='weight')
        mod_k = len(mod_communities)

        if mod_k <= 1:
            sub_clusters_with_parent_info.append((original_cluster_idx, main_cluster_nodes, False))
            continue

        # Step 2: Adaptive Spectral Clustering
        adj_matrix = nx.to_numpy_array(subgraph, weight='weight')
        node_list = list(subgraph.nodes())
        best_partition = None
        best_modularity_score = -1.0
        k = max(2, mod_k - 2)

        while k <= len(node_list):
            try:
                sc = SpectralClustering(
                    n_clusters=k, affinity='precomputed', n_init=10, assign_labels='kmeans'
                )
                labels = sc.fit_predict(adj_matrix)

                proposed_partition = collections.defaultdict(set)
                for i, node_label in enumerate(labels):
                    proposed_partition[node_label].add(node_list[i])

                proposed_partition_list = list(proposed_partition.values())
                mod_score = nx.community.modularity(subgraph, proposed_partition_list, weight='weight')

                if mod_score > best_modularity_score and mod_score > modularity_split_threshold:
                    best_modularity_score = mod_score
                    best_partition = proposed_partition_list
                    k += 1
                else:
                    break
            except Exception as e:
                print(f"Warning: Spectral failed for cluster {original_cluster_idx} at k={k}: {e}", file=sys.stderr)
                break

        if best_partition:
            # Ensure connectivity within proposed spectral clusters
            connected_partition = []
            for part_nodes in best_partition:
                part_subgraph = subgraph.subgraph(part_nodes)
                for component in nx.connected_components(part_subgraph):
                    connected_partition.append(component)

            current_communities = {i: set(nodes) for i, nodes in enumerate(connected_partition)}

            # Merge small clusters
            while True:
                small_ids = [cid for cid, nodes in current_communities.items() if len(nodes) < min_final_cluster_size]
                if not small_ids or len(current_communities) <= 1:
                    break

                small_ids.sort(key=lambda cid: len(current_communities[cid]))
                to_merge_id = small_ids[0]
                to_merge_nodes = current_communities[to_merge_id]

                best_target = None
                max_weight = -1
                for target_id, target_nodes in current_communities.items():
                    if target_id == to_merge_id:
                        continue
                    conn_weight = 0
                    for n in to_merge_nodes:
                        for nbr in subgraph.neighbors(n):
                            if nbr in target_nodes:
                                conn_weight += subgraph[n][nbr]['weight']
                    if best_target is None or conn_weight > max_weight:
                        max_weight = conn_weight
                        best_target = target_id

                if best_target is not None and max_weight > 0:
                    current_communities[best_target].update(to_merge_nodes)
                    del current_communities[to_merge_id]
                else:
                    break

            if len(current_communities) > 1:
                for community_set in current_communities.values():
                    sub_clusters_with_parent_info.append((original_cluster_idx, community_set, True))
            else:
                sub_clusters_with_parent_info.append((original_cluster_idx, main_cluster_nodes, False))
        else:
            sub_clusters_with_parent_info.append((original_cluster_idx, main_cluster_nodes, False))

    return sub_clusters_with_parent_info


def write_clusters_to_file(filename, clusters_data, main_cluster_id_map=None):
    """Write clustering results to a text file."""
    with open(filename, 'w') as f_out:
        if main_cluster_id_map:
            # Writing sub-clusters
            sorted_raw = sorted(clusters_data, key=lambda x: (x[0], sorted(list(x[1]))[0] if x[1] else ''))
            counts = collections.defaultdict(int)
            processed = []
            for p_idx, c_set, was_split in sorted_raw:
                p_id = main_cluster_id_map[p_idx]
                if was_split:
                    counts[p_id] += 1
                    full_id = f"{p_id}.{counts[p_id]:02d}"
                    processed.append((full_id, counts[p_id], c_set))
                else:
                    processed.append((p_id, 0, c_set))

            final_sorted = sorted(processed, key=lambda x: (x[0].split('.')[0], x[1]))
            for cid, _, c_set in final_sorted:
                f_out.write(f"{cid}\t{' '.join(sorted(list(c_set)))}\n")
        else:
            # Writing main clusters
            sorted_main = sorted(clusters_data, key=len, reverse=True)
            new_map = {}
            for i, c_set in enumerate(sorted_main):
                orig_idx = -1
                for idx, s in enumerate(clusters_data):
                    if s == c_set:
                        orig_idx = idx
                        break
                cid = f"Cluster{i+1:04d}"
                new_map[orig_idx] = cid
                f_out.write(f"{cid}\t{' '.join(sorted(list(c_set)))}\n")
            return new_map


def write_network_edges(filename, edges):
    """Write network edge list to a TSV file."""
    with open(filename, 'w') as f:
        f.write("Source\tTarget\tWeight\n")
        for n1, n2, w in edges:
            f.write(f"{n1}\t{n2}\t{w}\n")


def write_cytoscape_node_attributes(filename, sub_info, m_map, node_meta, meta_keys):
    """Export node attributes for Cytoscape visualization."""
    node_attrs = {}
    sorted_raw = sorted(sub_info, key=lambda x: (x[0], sorted(list(x[1]))[0] if x[1] else ''))
    counts = collections.defaultdict(int)
    for p_idx, c_set, was_split in sorted_raw:
        p_id = m_map[p_idx]
        if was_split:
            counts[p_id] += 1
            full_id = f"{p_id}.{counts[p_id]:02d}"
        else:
            full_id = p_id
        for node in c_set:
            node_attrs[node] = {'SubID': full_id, 'MainID': p_id, 'Split': was_split}

    headers = ["NodeID", "SubClusterID", "MainClusterID", "WasSplit"] + meta_keys
    with open(filename, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for nid in sorted(node_attrs.keys()):
            attrs = node_attrs[nid]
            row = [nid, attrs['SubID'], attrs['MainID'], str(attrs['Split'])]
            meta = node_meta.get(nid, {})
            for k in meta_keys:
                row.append(str(meta.get(k, 'N/A')))
            f.write('\t'.join(row) + '\n')


def write_cluster_stats(filename, main_c, sub_info, m_map):
    """Write summary statistics for each cluster."""
    stats = {m_map[i]: {'Size': len(c), 'Subs': 1} for i, c in enumerate(main_c)}
    sub_counts = collections.defaultdict(int)
    for p_idx, _, was_split in sub_info:
        if was_split:
            sub_counts[m_map[p_idx]] += 1
    for cid, count in sub_counts.items():
        stats[cid]['Subs'] = count
    with open(filename, 'w') as f:
        f.write("ClusterID\tSize\tNumberSubCluster\n")
        for cid in sorted(stats.keys()):
            stat = stats[cid]
            f.write(f"{cid}\t{stat['Size']}\t{stat['Subs']}\n")


def process_clustering(input_file, output_dir, min_size, min_nodes, threshold):
    """Orchestrate the clustering pipeline."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    edges, pairs, node_meta, meta_keys, meta_map = [], [], {}, [], {}

    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
            if not lines:
                return

            # Header detection and metadata extraction
            start_row = 0
            header_parts = lines[0].strip().split('\t')
            try:
                float(header_parts[2])
            except (ValueError, IndexError):
                start_row = 1
                for i in range(5, len(header_parts), 2):
                    m1 = re.match(r"element01_(.*)", header_parts[i])
                    m2 = re.match(r"element02_(.*)", header_parts[i+1]) if i+1 < len(header_parts) else None
                    if m1 and m2 and m1.group(1) == m2.group(1):
                        key = m1.group(1)
                        meta_keys.append(key)
                        meta_map[key] = (i, i+1)
                    else:
                        break

            for line in lines[start_row:]:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                n1, n2 = parts[0], parts[1]
                try:
                    weight = float(parts[2])
                    edges.append((n1, n2, weight))
                    pairs.append((n1, n2))
                    if meta_keys:
                        if n1 not in node_meta:
                            node_meta[n1] = {k: parts[meta_map[k][0]] for k in meta_keys if meta_map[k][0] < len(parts)}
                        if n2 not in node_meta:
                            node_meta[n2] = {k: parts[meta_map[k][1]] for k in meta_keys if meta_map[k][1] < len(parts)}
                except ValueError:
                    continue
    except Exception as e:
        sys.exit(f"Error during file processing: {e}")

    main_c = find_clusters(pairs)
    sub_info = find_sub_clusters_spectral(edges, main_c, min_size, min_nodes, threshold)

    m_map = write_clusters_to_file(os.path.join(output_dir, "main_clusters.txt"), main_c)
    write_clusters_to_file(os.path.join(output_dir, "sub_clusters.txt"), sub_info, main_cluster_id_map=m_map)
    write_network_edges(os.path.join(output_dir, "network_edges.txt"), edges)
    write_cytoscape_node_attributes(os.path.join(output_dir, "node_attributes.txt"), sub_info, m_map, node_meta, meta_keys)
    write_cluster_stats(os.path.join(output_dir, "cluster_stats.txt"), main_c, sub_info, m_map)


def main():
    """Main CLI execution block."""
    parser = argparse.ArgumentParser(description="Network clustering with Spectral and Modularity analysis.")
    parser.add_argument('-i', '--input-file', required=True, help="Input TSV file with edges.")
    parser.add_argument('-o', '--output-dir', default='./', help="Directory for output files.")
    parser.add_argument('-m', '--min-size', type=int, default=1, help="Minimum final cluster size.")
    parser.add_argument('-n', '--min-nodes', type=int, default=4, help="Min nodes to trigger spectral split.")
    parser.add_argument('-t', '--threshold', type=float, default=0.05, help="Modularity split threshold.")

    args = parser.parse_args()
    process_clustering(args.input_file, args.output_dir, args.min_size, args.min_nodes, args.threshold)


if __name__ == "__main__":
    main()