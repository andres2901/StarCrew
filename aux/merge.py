#!/usr/bin/env python3
"""
GFF Overlap Resolver.

Reduces GFF redundancy by identifying overlapping genes on the same 
strand and selecting only the longest gene.
"""

import gffutils
import sys
import os
import argparse
import collections

def get_maintained_gene_ids(gff_file_path):
    """
    Identifies and returns the IDs of the longest genes in overlapping clusters.
    """
    dbfn = f'temp_{os.path.basename(gff_file_path)}.db'
    
    try:
        if os.path.exists(dbfn):
             os.remove(dbfn)
             
        db = gffutils.create_db(
            gff_file_path,
            dbfn=dbfn,
            force=True,
            keep_order=True,
            id_spec='ID',
            merge_strategy='replace',
            disable_infer_transcripts=True,
            disable_infer_genes=True
        )
    except Exception as e:
        print(f"Error creating gffutils database: {e}", file=sys.stderr)
        return []

    maintained_gene_ids = []
    
    all_genes = list(db.features_of_type('gene'))
    genes_to_process = set(gene.id for gene in all_genes)
    
    while genes_to_process:
        root_gene_id = genes_to_process.pop()
        
        queue = collections.deque([root_gene_id])
        current_cluster = [db[root_gene_id]]
        
        while queue:
            current_id = queue.popleft()
            current_gene = db[current_id]
            
            overlapping_potential = db.features_of_type(
                'gene', 
                limit=(current_gene.seqid, current_gene.start, current_gene.end)
            )
            
            for gene2 in overlapping_potential:
                if (gene2.id in genes_to_process and 
                    gene2.id != current_gene.id and 
                    gene2.strand == current_gene.strand):
                    
                    current_cluster.append(gene2)
                    queue.append(gene2.id)
                    genes_to_process.remove(gene2.id)

        longest_gene = max(current_cluster, key=lambda x: x.end - x.start)
        maintained_gene_ids.append(longest_gene.id)

    if os.path.exists(dbfn):
        os.remove(dbfn)

    return sorted(maintained_gene_ids)

def main():
    """CLI Entry Point."""
    parser = argparse.ArgumentParser(
        description="Pick the longest gene from same-strand overlapping clusters."
    )
    parser.add_argument("input_gff", help="Path to the input GFF file.")
    parser.add_argument("output_file", help="Path to write the list of maintained IDs.")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_gff):
        print(f"Error: File '{args.input_gff}' not found.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Processing overlaps in {args.input_gff}...")
    maintained_ids = get_maintained_gene_ids(args.input_gff)
    
    if maintained_ids:
        with open(args.output_file, 'w') as out_file:
            for gene_id in maintained_ids:
                out_file.write(f"{gene_id}\n")
        print(f"Done. Kept {len(maintained_ids)} genes. IDs saved to: {args.output_file}")
    else:
        print("No gene features found to process.")

if __name__ == "__main__":
    main()