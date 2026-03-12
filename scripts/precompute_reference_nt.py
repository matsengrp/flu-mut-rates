"""Pre-compute reference nucleotide sequences for all segment/subtype combinations.

Outputs results/reference_nt.json with structure:
  { "SEGMENT:SUBTYPE": { "1": "A", "2": "C", ... }, ... }

Keys use string site positions (JSON keys must be strings).
"""

import json
import yaml
from Bio import SeqIO

with open("config.yaml") as f:
    config = yaml.safe_load(f)

DATA_DIR = config["data_dir"]

# segment -> list of subtypes to iterate over
SEGMENT_SUBTYPES = {
    "HA": config["ha_subtypes"],
    "NA": config["na_subtypes"],
    "NP": ["all"], "PA": ["all"], "PB1": ["all"], "PB2": ["all"],
    "MP": ["all"], "NS": ["all"],
}


def load_reference_nt(segment, subtype):
    fasta_path = f"{DATA_DIR}/{segment}/{subtype}/curated_reference.fasta"
    rec = next(SeqIO.parse(fasta_path, "fasta"))
    sites = {str(i + 1): nt for i, nt in enumerate(str(rec.seq))}
    return rec.id, sites


result = {}
accessions = {}
for segment, subtypes in SEGMENT_SUBTYPES.items():
    for subtype in subtypes:
        key = f"{segment}:{subtype}"
        print(f"Processing {key}...")
        acc, nt = load_reference_nt(segment, subtype)
        result[key] = nt
        accessions[key] = acc

result["_accessions"] = accessions

with open("results/reference_nt.json", "w") as f:
    json.dump(result, f)

print(f"Done. Wrote {len(result) - 1} entries to results/reference_nt.json")
