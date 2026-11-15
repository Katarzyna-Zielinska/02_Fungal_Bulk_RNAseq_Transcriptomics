#!/usr/bin/env python3

import csv
import re
from pathlib import Path

target_file = Path(
    "results/go_enrichment/biological_interpretation/"
    "annotation_audit/genes_requiring_further_annotation.csv"
)

fasta_file = Path("data/reference/proteins/candida_auris_proteins.faa")

eggnog_file = Path(
    "results/annotation/eggnog/candida_auris.emapper.annotations"
)

output_dir = Path(
    "results/go_enrichment/biological_interpretation/"
    "annotation_audit/eggnog_refined"
)

output_dir.mkdir(parents=True, exist_ok=True)

def normalize_id(text):
    match = re.search(r"B9J08_\d+", text)
    return match.group(0) if match else None

# --------------------------------------------------
# 1. LOAD TARGET GENES
# --------------------------------------------------

with target_file.open(newline="") as handle:
    targets = list(csv.DictReader(handle))

# --------------------------------------------------
# 2. LOAD PROTEIN FASTA
# --------------------------------------------------

proteins = {}
header = None
sequence = []

with fasta_file.open() as handle:
    for line in handle:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if header is not None:
                proteins[normalize_id(header)] = {
                    "fasta_header": header,
                    "protein_sequence": "".join(sequence)
                }

            header = line[1:]
            sequence = []
        else:
            sequence.append(line)

if header is not None:
    proteins[normalize_id(header)] = {
        "fasta_header": header,
        "protein_sequence": "".join(sequence)
    }

# --------------------------------------------------
# 3. LOAD EGGNOG ANNOTATIONS
# --------------------------------------------------

eggnog = {}

with eggnog_file.open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 21:
            continue

        gene_id = normalize_id(fields[0])

        if gene_id:
            eggnog[gene_id] = {
                "eggnog_query": fields[0],
                "seed_ortholog": fields[1],
                "evalue": fields[2],
                "score": fields[3],
                "eggNOG_OGs": fields[4],
                "max_annot_lvl": fields[5],
                "eggNOG_COG_category": fields[6],
                "eggNOG_description": fields[7],
                "eggNOG_Preferred_name": fields[8],
                "eggNOG_GOs": fields[9],
                "eggNOG_EC": fields[10],
                "eggNOG_KEGG_ko": fields[11],
                "eggNOG_KEGG_Pathway": fields[12],
                "eggNOG_KEGG_Module": fields[13],
                "eggNOG_KEGG_Reaction": fields[14],
                "eggNOG_KEGG_rclass": fields[15],
                "eggNOG_BRITE": fields[16],
                "eggNOG_KEGG_TC": fields[17],
                "eggNOG_CAZy": fields[18],
                "eggNOG_BiGG_Reaction": fields[19],
                "eggNOG_PFAMs": fields[20]
            }

# --------------------------------------------------
# 4. MERGE DATA
# --------------------------------------------------

results = []

for target in targets:
    gene_id = target["gene_id"]
    normalized = normalize_id(gene_id)

    protein = proteins.get(normalized)
    annotation = eggnog.get(normalized)

    if protein is None or annotation is None:
        raise ValueError(
            f"Missing protein or eggNOG annotation: {gene_id}"
        )

    record = dict(target)

    record["normalized_gene_id"] = normalized
    record["fasta_header"] = protein["fasta_header"]
    record["protein_length_aa"] = len(
        protein["protein_sequence"]
    )

    record.update(annotation)

    results.append(record)

# --------------------------------------------------
# 5. WRITE CSV
# --------------------------------------------------

output_file = output_dir / "18_core_genes_refined_eggnog_annotation.csv"

fieldnames = list(results[0].keys())

with output_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames
    )
    writer.writeheader()
    writer.writerows(results)

# --------------------------------------------------
# 6. VALIDATION
# --------------------------------------------------

print("\n========== REFINED ANNOTATION SUMMARY ==========")

print("Target genes:", len(targets))
print("Merged genes:", len(results))
print(
    "Genes with PFAM:",
    sum(
        bool(r["eggNOG_PFAMs"]) and r["eggNOG_PFAMs"] != "-"
        for r in results
    )
)
print(
    "Genes with GO:",
    sum(
        bool(r["eggNOG_GOs"]) and r["eggNOG_GOs"] != "-"
        for r in results
    )
)
print(
    "Genes with preferred name:",
    sum(
        bool(r["eggNOG_Preferred_name"])
        and r["eggNOG_Preferred_name"] != "-"
        for r in results
    )
)

print("\nOutput:")
print(output_file)

if len(results) != 18:
    raise SystemExit("FAIL: Expected 18 genes.")

if len({r["gene_id"] for r in results}) != 18:
    raise SystemExit("FAIL: Duplicate gene IDs detected.")

print("\nPASS: Refined eggNOG annotation table generated.")
