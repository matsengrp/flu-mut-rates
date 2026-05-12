import argparse
import json
import math
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import sys
import gffutils
import os

def main():
    parser = argparse.ArgumentParser(description="Create coding sites file from reference sequence.")
    parser.add_argument('--ref_fasta', required=True, help='Path to reference FASTA file')
    parser.add_argument('--output_directory', required=True, help='Directory for output files')
    parser.add_argument('--gff_file', required=True, help='Path to GFF annotation file')
    parser.add_argument('--subtype', required=True, help='Influenza subtype (e.g., H1N1, H3N2)')
    parser.add_argument('--segment', required=True, help='Gene segment (e.g., HA, NA)')
    parser.add_argument('--ignore_genes', nargs='*', default=[],
                        help='Gene names to exclude from coding site mapping (e.g. PB1-F2)')
    parser.add_argument('--additional_orfs', default='[]',
                        help='JSON-encoded list of {name, start, end} dicts for ORFs '
                             'not present in the GFF (e.g. PA X-ORF). Treated as '
                             'single-CDS overlapping genes.')
    args = parser.parse_args()
    additional_orfs = json.loads(args.additional_orfs)

    # Read in the reference sequence from the FASTA file, assuming it is the first sequence
    with open(args.ref_fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            ref_seq = str(record.seq)
            ref_id = record.id
            break

    # Process GFF file to extract gene and CDS information
    try:
        db = gffutils.create_db(args.gff_file, ':memory:', force=True, merge_strategy='create_unique')
    except Exception as e:
        sys.stderr.write(f"Error processing GFF file: {e}\n")
        sys.exit(1)
    
    # Get all genes and their associated CDSs
    genes_info = {}
    
    for gene_feature in db.features_of_type('gene'):
        gene_id = gene_feature.attributes.get('gene', [gene_feature.id])[0]
        gene_name = gene_feature.attributes.get('Name', [gene_id])[0]
        
        # Find all CDSs for this gene
        cdss = []
        try:
            # Look for CDSs that are children of this gene
            for cds in db.children(gene_feature, featuretype='CDS'):
                cdss.append((cds.start, cds.end, cds.strand))
        except:
            # If no parent-child relationship, find CDSs by gene attribute
            for cds in db.features_of_type('CDS'):
                cds_gene = cds.attributes.get('gene', [''])[0]
                if cds_gene == gene_id:
                    cdss.append((cds.start, cds.end, cds.strand))
        
        if cdss:
            # Sort CDSs by start position
            cdss.sort(key=lambda x: x[0])
            genes_info[gene_name] = {
                'cdss': cdss,
                'strand': cdss[0][2]  # Assume all CDSs have same strand
            }
    
    if args.ignore_genes:
        genes_info = {name: data for name, data in genes_info.items() if name not in args.ignore_genes}
        print(f"Ignoring genes: {args.ignore_genes}")

    for orf in additional_orfs:
        genes_info[orf['name']] = {
            'cdss': [(int(orf['start']), int(orf['end']), '+')],
            'strand': '+',
        }
    if additional_orfs:
        print(f"Adding {len(additional_orfs)} non-GFF ORFs: "
              f"{[orf['name'] for orf in additional_orfs]}")

    print(f"Found {len(genes_info)} genes: {list(genes_info.keys())}")
    
    # Create output directory and protein sequences subdirectory if they don't exist
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    
    protein_seq_dir = os.path.join(args.output_directory, "protein_sequences")
    if not os.path.exists(protein_seq_dir):
        os.makedirs(protein_seq_dir)
    
    # Translate each gene and save as FASTA
    for gene_name, gene_info in genes_info.items():
        cdss = gene_info['cdss']
        strand = gene_info['strand']
        
        # Extract and concatenate all CDS sequences for this gene
        gene_nucleotide_seq = ""
        for cds_start, cds_end, cds_strand in cdss:
            # Extract sequence (convert to 0-based indexing)
            cds_seq = ref_seq[cds_start-1:cds_end]
            gene_nucleotide_seq += cds_seq
        
        # Convert to Biopython Seq object
        bio_seq = Seq(gene_nucleotide_seq)
        
        # Reverse complement if on negative strand
        if strand == '-':
            bio_seq = bio_seq.reverse_complement()
        
        # Check if length is divisible by 3
        if len(bio_seq) % 3 != 0:
            print(f"Warning: Gene {gene_name} nucleotide sequence length ({len(bio_seq)}) is not divisible by 3")
            # Trim to nearest multiple of 3
            bio_seq = bio_seq[:len(bio_seq) - (len(bio_seq) % 3)]
        
        # Translate to amino acids
        try:
            amino_acid_seq = bio_seq.translate()
            
            # Create SeqRecord with custom ID format and blank description
            record = SeqRecord(
                amino_acid_seq,
                id=f"{args.subtype} {args.segment} {gene_name.upper()}",
                description=""
            )
            
            # Save to FASTA file in protein_sequences subdirectory
            output_fasta = os.path.join(protein_seq_dir, f"{gene_name.upper()}.fasta")
            SeqIO.write(record, output_fasta, "fasta")
            print(f"Saved translated sequence for {gene_name.upper()} to {output_fasta} ({len(amino_acid_seq)} amino acids)")
            
        except Exception as e:
            print(f"Error translating gene {gene_name}: {e}")
    
    # Create a dataframe to store coding sites information
    coding_sites_data = []
    
    # Process each site in the reference sequence
    for site in range(1, len(ref_seq) + 1):
        
        # Check which genes overlap this site
        overlapping_genes = []
        
        for gene_name, gene_info in genes_info.items():
            cdss = gene_info['cdss']
            
            # Check if site falls within any CDS of this gene
            for cds_start, cds_end, strand in cdss:
                if cds_start <= site <= cds_end:
                    overlapping_genes.append((gene_name, gene_info))
                    break
        
        # If site overlaps with genes, calculate codon position and site
        for gene_name, gene_info in overlapping_genes:
            cdss = gene_info['cdss']
            
            # Calculate cumulative position within the gene's coding sequence
            cumulative_pos = 0
            found = False
            
            for cds_start, cds_end, strand in cdss:
                if cds_start <= site <= cds_end:
                    # Position within this CDS
                    pos_in_cds = site - cds_start + 1
                    cumulative_pos += pos_in_cds
                    found = True
                    break
                else:
                    # Add length of this CDS if it comes before our site
                    if cds_end < site:
                        cumulative_pos += cds_end - cds_start + 1
            
            if found:
                # Calculate codon position (1, 2, or 3)
                rem = cumulative_pos % 3
                codon_pos = 3 if rem == 0 else rem
                
                # Calculate codon site number
                codon_site_num = math.ceil(cumulative_pos / 3)
                
                # Create entry for this gene
                coding_sites_data.append({
                    'site': site,
                    'codon_position': codon_pos,
                    'codon_site': codon_site_num,
                    'gene': gene_name.upper()
                })
        
        # If no overlapping genes, add noncoding entry
        if not overlapping_genes:
            coding_sites_data.append({
                'site': site,
                'codon_position': 'noncoding',
                'codon_site': 'noncoding',
                'gene': 'noncoding'
            })
    
    # Convert to DataFrame
    coding_sites_df = pd.DataFrame(coding_sites_data)
    
    # Compress multiple rows per site using groupby
    final_df = coding_sites_df.groupby('site').agg({
        'codon_position': lambda x: ';'.join(str(val) for val in x),
        'codon_site': lambda x: ';'.join(str(val) for val in x),
        'gene': lambda x: ';'.join(str(val) for val in x)
    }).reset_index()
    
    # Save coding sites CSV in output directory
    coding_sites_output = os.path.join(args.output_directory, "coding_sites.csv")
    final_df.to_csv(coding_sites_output, index=False)
    
    print(f"Created coding sites file: {coding_sites_output}")
    print(f"Processed sequence '{ref_id}': {len(ref_seq)} nucleotides")
    
    # Print summary statistics (using original uncompressed data for accurate counts)
    gene_counts = coding_sites_df['gene'].value_counts()
    print("Gene coverage summary:")
    for gene, count in gene_counts.items():
        if gene != 'noncoding':
            codons = count / 3
            print(f"  {gene}: {count} nucleotides ({codons:.1f} codons)")
        else:
            print(f"  {gene}: {count} nucleotides")
    
    # Check for genes with lengths not divisible by 3
    for gene_name in genes_info.keys():
        gene_sites = coding_sites_df[coding_sites_df['gene'] == gene_name]
        if len(gene_sites) % 3 != 0:
            print(f"Warning: Gene {gene_name} has {len(gene_sites)} nucleotides (not divisible by 3)")

if __name__ == "__main__":
    main()