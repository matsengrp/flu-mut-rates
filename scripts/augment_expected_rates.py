"""
Augment expected_rates_by_predictor.csv with CG and GC mutation types.

This script adds missing CG and GC mutation types to the full neutral model
expected rates, using rates from segment_wide_rates.csv, and outputs a
complete table to a new file.
"""

import argparse
import pandas as pd
import sys


def load_segment_rates(filepath):
    """Load segment-wide rates and filter to CG/GC synonymous mutations."""
    df = pd.read_csv(filepath, keep_default_na=False)
    filtered = df.query("mut_type in ['CG', 'GC'] and mut_class == 'synonymous' and host == 'all'")
    return filtered[['mut_type', 'segment', 'rate']]


def generate_motifs(mut_type):
    """
    Generate all 16 possible 3-mer motifs for a mutation type.
    Center nucleotide must match the wild-type (first letter of mut_type).
    """
    letters = ['A', 'C', 'G', 'T']
    center = mut_type[0]  # Wild-type nucleotide
    motifs = [f"{left}{center}{right}" for left in letters for right in letters]
    return sorted(motifs)


def augment_full_model(df_existing, segment_rates):
    """Augment full model with CG and GC entries (16 motifs × 8 segments each)."""
    # Segment order matching rates_model.py
    segments = ['PB1', 'NS', 'NA', 'HA', 'PB2', 'MP', 'PA', 'NP']

    new_rows = []
    for mut_type in ['CG', 'GC']:
        # Get segment-specific rates
        mut_rates = segment_rates[segment_rates['mut_type'] == mut_type]

        # Generate motifs
        motifs = generate_motifs(mut_type)

        # Create row for each segment × motif combination
        for segment in segments:
            segment_rate = mut_rates[mut_rates['segment'] == segment]['rate'].values
            if len(segment_rate) == 0:
                print(f"Warning: No rate found for {mut_type} in {segment}, skipping segment")
                continue

            rate = segment_rate[0]

            for motif in motifs:
                new_rows.append({
                    'mut_type': mut_type,
                    'segment': segment,
                    'motif': motif,
                    'predicted_rate': rate,
                    'tau_squared': 'n.a.'
                })

    df_augmented = pd.concat([df_existing, pd.DataFrame(new_rows)], ignore_index=True)
    df_augmented = df_augmented.sort_values(['mut_type', 'segment', 'motif']).reset_index(drop=True)
    return df_augmented


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description='Augment full model expected rates with CG and GC mutation types'
    )
    parser.add_argument('--segment_wide_rates', required=True,
                       help='Path to segment_wide_rates.csv')
    parser.add_argument('--input_file', required=True,
                       help='Path to full model expected_rates_by_predictor.csv')
    parser.add_argument('--output_file', required=True,
                       help='Path to output expected_rates.csv')

    args = parser.parse_args()

    # Load segment-wide rates
    print(f"Loading segment-wide rates from {args.segment_wide_rates}")
    try:
        segment_rates = load_segment_rates(args.segment_wide_rates)
    except FileNotFoundError:
        print(f"Error: File not found: {args.segment_wide_rates}")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading segment_wide_rates: {e}")
        sys.exit(1)

    # Validate we have CG and GC data
    mut_types = set(segment_rates['mut_type'].unique())
    if 'CG' not in mut_types or 'GC' not in mut_types:
        print(f"Error: Missing CG or GC in segment_wide_rates. Found: {mut_types}")
        sys.exit(1)

    # Validate we have all 8 segments
    segments = set(segment_rates['segment'].unique())
    expected_segments = {'PB1', 'NS', 'NA', 'HA', 'PB2', 'MP', 'PA', 'NP'}
    if not expected_segments.issubset(segments):
        missing = expected_segments - segments
        print(f"Error: Missing segments in data. Missing: {missing}")
        sys.exit(1)

    print(f"Found CG and GC rates for all 8 segments")

    # Load and augment full model
    print(f"\nLoading full model from: {args.input_file}")
    try:
        df_full = pd.read_csv(args.input_file, keep_default_na=False)
        print(f"  Original rows: {len(df_full)}")
        df_full_aug = augment_full_model(df_full, segment_rates)
        print(f"  Augmented rows: {len(df_full_aug)}")
        df_full_aug.to_csv(args.output_file, index=False)
        print(f"  Saved to {args.output_file}")
    except Exception as e:
        print(f"Error augmenting full model: {e}")
        sys.exit(1)

    # Validation: Check each mutation type has same number of entries
    print("\n" + "="*60)
    print("Validating entry counts...")
    print("="*60)

    # Full model: each mut_type should have exactly 128 rows (16 motifs × 8 segments)
    full_counts = df_full_aug['mut_type'].value_counts()
    if not all(count == 128 for count in full_counts.values):
        print("ERROR: Mutation types have unequal entry counts:")
        print(full_counts)
        sys.exit(1)
    print(f"✓ All {len(full_counts)} mutation types have 128 entries each")

    # Full model: each mut_type × segment combo should have exactly 16 rows (one per motif)
    for mut_type in df_full_aug['mut_type'].unique():
        for segment in ['PB1', 'NS', 'NA', 'HA', 'PB2', 'MP', 'PA', 'NP']:
            count = len(df_full_aug[(df_full_aug['mut_type'] == mut_type) &
                                     (df_full_aug['segment'] == segment)])
            if count != 16:
                print(f"ERROR: {mut_type} in {segment} has {count} entries, expected 16")
                sys.exit(1)
    print(f"✓ All mutation type × segment combinations have 16 motifs each")

    print("\n" + "="*60)
    print("Augmentation complete!")
    print("="*60)
    print(f"  Rows: {len(df_full)} -> {len(df_full_aug)} (added {len(df_full_aug) - len(df_full)})")
    print("="*60)


if __name__ == "__main__":
    main()
