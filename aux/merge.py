import gffutils
import sys
import os
import argparse
import collections

def get_maintained_gene_ids(gff_file_path):
    """Identifies and returns the IDs of the genes to be maintained"""
    dbfn = 'temp_gff.db'
    
    # --- 1. Create or connect to the gffutils database ---
    try:
        # Check if the database already exists and delete it if forced
        if os.path.exists(dbfn):
             os.remove(dbfn)
             
        db = gffutils.create_db(
            gff_file_path,
            dbfn=dbfn,
            force=True, # Will now re-create the DB if it exists
            keep_order=True,
            # Add an index for 'gene' features by start/end for faster querying
            id_spec='ID',
            merge_strategy='replace', # Default, but good to be explicit
            disable_infer_transcripts=True, # Optional, often useful if only dealing with genes
            disable_infer_genes=True # Optional, often useful if only dealing with genes
        )
    except Exception as e:
        print(f"Error creating gffutils database: {e}", file=sys.stderr)
        return []

    # --- 2. Initial Setup ---
    maintained_gene_ids = []
    
    # Fetch all gene features from the DB
    all_genes = list(db.features_of_type('gene'))
    # Use a set of IDs for genes that still need to be processed
    genes_to_process = set(gene.id for gene in all_genes)
    
    # --- 3. Iterative Overlap Group Processing (Cluster Finding) ---
    while genes_to_process:
        # Start a new potential cluster with an un-processed gene
        gene1_id = genes_to_process.pop()
        
        # Use a queue for a Breadth-First Search (BFS) style traversal
        # to find all connected genes in the cluster
        queue = collections.deque([gene1_id])
        current_cluster = [db[gene1_id]]
        
        # Flag to check if we found any overlap at all
        found_overlap = False 
        
        # Process the queue until all connected genes are found
        while queue:
            current_id = queue.popleft()
            current_gene = db[current_id]
            
            # Find *all* genes that overlap the current_gene
            # The limit argument in features_of_type speeds up the search
            overlapping_potential = list(db.features_of_type(
                'gene', 
                limit=(current_gene.seqid, current_gene.start, current_gene.end)
            ))
            
            for gene2 in overlapping_potential:
                # Check for overlap, same strand, and if it's a new gene to the cluster
                # gffutils overlap: (a.start <= b.end) and (a.end >= b.start)
                # The 'limit' ensures they are on the same seqid and broadly overlap
                
                # Check if it's a gene we haven't added to the cluster yet
                is_new_to_cluster = gene2.id != current_gene.id and gene2.id in genes_to_process
                
                # Further check: same strand and actual overlap (limit is only a hint)
                # gffutils features_of_type with limit is usually sufficient for actual overlap
                
                if is_new_to_cluster and current_gene.strand == gene2.strand:
                    # Found an overlapping gene connected to the cluster
                    found_overlap = True
                    
                    # Add to the cluster list
                    current_cluster.append(gene2)
                    
                    # Add to the queue for checking its own overlaps later
                    queue.append(gene2.id)
                    
                    # Remove from the main set to ensure it's not started as a new cluster later
                    genes_to_process.remove(gene2.id)

        # --- 4. Select the Longest Gene from the Cluster ---
        
        # If found_overlap is True, it means current_cluster has more than one gene
        # The cluster contains all overlapping genes (if any)
        
        longest_gene_in_group = max(current_cluster, key=lambda x: x.end - x.start)
        maintained_gene_ids.append(longest_gene_in_group.id)


    # --- 5. Cleanup ---
    if os.path.exists(dbfn):
        os.remove(dbfn)

    # Return the maintained IDs, sorted for consistent output
    return sorted(maintained_gene_ids)

def main():
    parser = argparse.ArgumentParser(description="Find the longest gene within overlapping groups and write their IDs to a file.")
    parser.add_argument("input_gff", help="Path to the input GFF file.")
    parser.add_argument("output_file", help="Path to the output file for maintained gene IDs.")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_gff):
        print(f"Error: Input file '{args.input_gff}' not found.", file=sys.stderr)
        sys.exit(1)
    
    maintained_ids = get_maintained_gene_ids(args.input_gff)
    
    if maintained_ids:
        with open(args.output_file, 'w') as out_file:
            for gene_id in maintained_ids:
                out_file.write(f"{gene_id}\n")
        print(f"Maintained gene IDs have been written to {args.output_file}")
    else:
        print("No maintained gene IDs found to write.")

if __name__ == "__main__":
    main()