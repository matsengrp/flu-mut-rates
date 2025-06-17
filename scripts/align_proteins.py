import argparse
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description="Align protein sequences across subtypes using MUSCLE.")
    parser.add_argument('--output_dir', required=True, help='Base output directory containing subtype/segment subdirectories')
    parser.add_argument('--segment', required=True, help='Gene segment (e.g., HA, NA)')
    parser.add_argument('--subtypes', nargs='+', required=True, help='List of subtypes to process')
    parser.add_argument('--muscle_path', default='muscle', help='Path to MUSCLE executable (default: muscle)')
    args = parser.parse_args()

    # Dictionary to store genes found in each subtype
    subtype_genes = {}
    
    # Check each subtype's protein_sequences directory
    for subtype in args.subtypes:
        protein_dir = os.path.join(args.output_dir, subtype, args.segment, "protein_sequences")
        
        if not os.path.exists(protein_dir):
            sys.stderr.write(f"Error: Directory {protein_dir} does not exist\n")
            sys.exit(1)
        
        # Get list of genes (FASTA files) in this directory
        fasta_files = [f for f in os.listdir(protein_dir) if f.endswith('.fasta')]
        genes = [os.path.splitext(f)[0] for f in fasta_files]  # Remove .fasta extension
        
        subtype_genes[subtype] = set(genes)
        print(f"Found {len(genes)} genes in {subtype} {args.segment}: {sorted(genes)}")
    
    # Check that all subtypes have the same set of genes
    all_gene_sets = list(subtype_genes.values())
    common_genes = set.intersection(*all_gene_sets) if all_gene_sets else set()
    
    if len(common_genes) == 0:
        sys.stderr.write("Error: No common genes found across all subtypes\n")
        sys.exit(1)
    
    # Check if any subtypes are missing genes
    for subtype, genes in subtype_genes.items():
        missing_genes = common_genes - genes
        extra_genes = genes - common_genes
        
        if missing_genes:
            print(f"Warning: {subtype} is missing genes: {sorted(missing_genes)}")
        if extra_genes:
            print(f"Warning: {subtype} has extra genes: {sorted(extra_genes)}")
    
    print(f"Processing {len(common_genes)} common genes: {sorted(common_genes)}")
    
    # Create output directory for aligned sequences
    aligned_dir = os.path.join(args.output_dir, "aligned_proteins", args.segment)
    if not os.path.exists(aligned_dir):
        os.makedirs(aligned_dir)
    
    # Process each gene
    for gene in sorted(common_genes):
        print(f"Processing gene: {gene}")
        
        # Collect sequences for this gene from all subtypes
        sequences = []
        
        for subtype in args.subtypes:
            protein_file = os.path.join(args.output_dir, subtype, args.segment, "protein_sequences", f"{gene}.fasta")
            
            if os.path.exists(protein_file):
                try:
                    # Read the sequence
                    record = SeqIO.read(protein_file, "fasta")
                    
                    # Create new record with updated ID (just use subtype as ID for clarity)
                    new_record = SeqRecord(
                        record.seq,
                        id=subtype,
                        description=f"{subtype} {args.segment} {gene}"
                    )
                    sequences.append(new_record)
                    
                except Exception as e:
                    print(f"Warning: Error reading {protein_file}: {e}")
            else:
                print(f"Warning: File {protein_file} not found")
        
        if len(sequences) < 2:
            print(f"Warning: Only found {len(sequences)} sequences for gene {gene}, skipping alignment")
            continue
        
        # Write unaligned sequences to temporary file
        temp_input = os.path.join(aligned_dir, f"{gene}_temp_input.fasta")
        SeqIO.write(sequences, temp_input, "fasta")
        
        # Define output file for aligned sequences
        aligned_output = os.path.join(aligned_dir, f"{gene}_aligned.fasta")
        
        # Run MUSCLE alignment
        try:
            cmd = [args.muscle_path, "-align", temp_input, "-output", aligned_output]
            print(f"Running: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Clean up temporary file
            os.remove(temp_input)
            
            print(f"Successfully aligned {len(sequences)} sequences for gene {gene}")
            print(f"Aligned sequences saved to: {aligned_output}")
            
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"Error running MUSCLE for gene {gene}: {e}\n")
            sys.stderr.write(f"MUSCLE stderr: {e.stderr}\n")
            # Clean up temporary file even if MUSCLE fails
            if os.path.exists(temp_input):
                os.remove(temp_input)
            continue
            
        except FileNotFoundError:
            sys.stderr.write(f"Error: MUSCLE executable not found at {args.muscle_path}\n")
            sys.stderr.write("Please install MUSCLE or specify the correct path with --muscle_path\n")
            # Clean up temporary file
            if os.path.exists(temp_input):
                os.remove(temp_input)
            sys.exit(1)
    
    print(f"Alignment complete. Aligned sequences saved in: {aligned_dir}")

if __name__ == "__main__":
    main()