"""Pre-compute reference amino acid sequences for all gene/subtype combinations.

Outputs results/reference_aa.json with structure:
  { "GENE:SUBTYPE": { "1": "M", "2": "A", ... }, ... }

Keys use string codon sites (JSON keys must be strings).
"""

import json
import yaml
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq

with open("config.yaml") as f:
    config = yaml.safe_load(f)

DATA_DIR = config["data_dir"]

GENE_TO_SEGMENT = {
    "HA": "HA", "NA": "NA", "NP": "NP",
    "PA": "PA", "PB1": "PB1", "PB2": "PB2",
    "M1": "MP", "M2": "MP",
    "NS1": "NS", "NEP": "NS",
}

# gene -> list of subtypes to iterate over
GENE_SUBTYPES = {
    "HA": config["ha_subtypes"],
    "NA": config["na_subtypes"],
    "NP": ["all"], "PA": ["all"], "PB1": ["all"], "PB2": ["all"],
    "M1": ["all"], "M2": ["all"],
    "NS1": ["all"], "NEP": ["all"],
}


def load_reference_aa(segment, subtype, gene):
    base = f"{DATA_DIR}/{segment}/{subtype}"
    db = gffutils.create_db(
        f"{base}/curated_reference.gff",
        ":memory:",
        force=True,
        merge_strategy="create_unique",
    )
    rec = next(SeqIO.parse(f"{base}/curated_reference.fasta", "fasta"))
    ref_seq = str(rec.seq)

    cdss = []
    for cds in db.features_of_type("CDS"):
        if cds.attributes.get("gene", [""])[0].upper() == gene.upper():
            cdss.append((cds.start, cds.end))
    cdss.sort()

    if not cdss:
        return None, {}

    cds_nt = "".join(ref_seq[s - 1 : e] for s, e in cdss)
    cds_nt = cds_nt[: len(cds_nt) - len(cds_nt) % 3]
    protein = str(Seq(cds_nt).translate(to_stop=True))

    return rec.id, {str(i + 1): aa for i, aa in enumerate(protein)}


result = {}
accessions = {}
for gene, subtypes in GENE_SUBTYPES.items():
    segment = GENE_TO_SEGMENT[gene]
    for subtype in subtypes:
        key = f"{gene}:{subtype}"
        print(f"Processing {key}...")
        acc, aa = load_reference_aa(segment, subtype, gene)
        if aa:
            result[key] = aa
            accessions[key] = acc
        else:
            print(f"  WARNING: no CDS found for {key}")

result["_accessions"] = accessions

with open("results/reference_aa.json", "w") as f:
    json.dump(result, f)

print(f"Done. Wrote {len(result) - 1} entries to results/reference_aa.json")
