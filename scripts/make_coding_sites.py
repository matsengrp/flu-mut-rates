#!/usr/bin/env python3
# filepath: /fh/fast/matsen_e/hhaddox/2025/flu-syn-rates/scripts/make_coding_sites.py
import argparse
import math
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import sys
import lzma

def main():
    parser = argparse.ArgumentParser(description="Create coding sites file from reference sequence.")
    parser.add_argument('--ref_fasta', required=True, help='Path to reference FASTA file')
    parser.add_argument('--output', required=True, help='Path to output CSV file')
    parser.add_argument('--gene_name', required=True, help='Name of the gene')
    args = parser.parse_args()

    # Check if the file is compressed with xz
    is_xz = args.ref_fasta.endswith('.xz')
    
    # Open the file with appropriate method
    try:
        if is_xz:
            fasta_file = lzma.open(args.ref_fasta, 'rt')  # Open in text mode for SeqIO
            print(f"Reading compressed .xz file: {args.ref_fasta}")
        else:
            fasta_file = open(args.ref_fasta, 'r')
            
        # Read first sequence from the FASTA file
        for record in SeqIO.parse(fasta_file, 'fasta'):
            ref_seq = str(record.seq)
            ref_id = record.id
            break  # Only use the first sequence
        
        fasta_file.close()
        
    except Exception as e:
        sys.stderr.write(f"Error opening or reading FASTA file: {e}\n")
        sys.exit(1)
    
    # Check if sequence length is a multiple of 3
    if len(ref_seq) % 3 != 0:
        sys.stderr.write(f"Error: Reference sequence length ({len(ref_seq)}) is not a multiple of 3.\n")
        sys.exit(1)
    
    # Create a dataframe to store coding sites information
    coding_sites_dict = defaultdict(list)
    
    # Process each site in the reference sequence
    for site in range(1, len(ref_seq) + 1):
        coding_sites_dict["site"].append(site)
        
        # Calculate position within codon (1, 2, or 3)
        offset = site  # Since we assume the sequence starts at position 1
        rem = offset % 3
        coding_sites_dict["codon_position"].append(3 if rem == 0 else rem)
        
        # Calculate which codon this site belongs to
        coding_sites_dict["codon_site"].append(math.ceil(offset / 3))
        
        # Assign gene name from command-line argument
        coding_sites_dict["gene"].append(args.gene_name)
    
    # Convert to DataFrame and save to CSV
    coding_sites_df = pd.DataFrame(coding_sites_dict)
    coding_sites_df.to_csv(args.output, index=False)
    print(f"Created coding sites file: {args.output}")
    print(f"Processed sequence '{ref_id}': {len(ref_seq)} nucleotides ({len(ref_seq)/3} codons)")
    print(f"Assigned gene name: {args.gene_name}")

if __name__ == "__main__":
    main()