import gffutils
import sys
import os
import argparse

def get_maintained_gene_ids(gff_file_path):
    """
    Identifies and returns the IDs of the genes to be maintained:
    - All non-overlapping genes.
    - The longest gene from each group of overlapping genes.
    """
    try:
        dbfn = 'temp_gff.db'
        db = gffutils.create_db(
            gff_file_path,
            dbfn=dbfn,
            force=True,
            keep_order=True
        )
    except Exception as e:
        print(f"Error creating gffutils database: {e}", file=sys.stderr)
        return []

    maintained_gene_ids = []
    
    all_genes = list(db.features_of_type('gene', order_by='start'))
    genes_to_process = set(gene.id for gene in all_genes)
    
    while genes_to_process:
        gene1_id = genes_to_process.pop()
        gene1 = db[gene1_id]
        
        overlapping_genes = [
            gene2 for gene2 in db.features_of_type('gene', limit=(gene1.seqid, gene1.start, gene1.end))
            if gene1.id != gene2.id and gene1.strand == gene2.strand and gene2.id in genes_to_process
        ]
        
        current_group = [gene1]
        
        if overlapping_genes:
            for gene2 in overlapping_genes:
                current_group.append(gene2)
                if gene2.id in genes_to_process:
                    genes_to_process.remove(gene2.id)
            
            longest_gene_in_group = max(current_group, key=lambda x: x.end - x.start)
            maintained_gene_ids.append(longest_gene_in_group.id)
        else:
            maintained_gene_ids.append(gene1.id)

    if os.path.exists(dbfn):
        os.remove(dbfn)

    return sorted(maintained_gene_ids)

# -------------------------------------------------------------
# Main execution block
# -------------------------------------------------------------

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