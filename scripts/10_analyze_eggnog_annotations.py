
#!/usr/bin/env python3

"""
Functional annotation analysis of Candida auris proteins
using eggNOG-mapper and NCBI gene annotation.

Input:
    - candida_auris.emapper.annotations
    - protein_gene_mapping.csv

Output:
    - candida_auris_gene_functional_annotation.csv
    - eggnog_annotation_summary.tsv
    - eggnog_annotation_summary.json
    - eggnog_unannotated_proteins.csv
"""

from pathlib import Path
import json

import pandas as pd
import numpy as np


# ============================================================
# 1. PROJECT PATHS
# ============================================================

PROJECT = Path(__file__).resolve().parents[1]

ANNOTATION_FILE = (
    PROJECT
    / "results/annotation/eggnog/"
    / "candida_auris.emapper.annotations"
)

MAPPING_FILE = (
    PROJECT
    / "data/annotation/"
    / "protein_gene_mapping.csv"
)

OUTPUT_DIR = (
    PROJECT
    / "results/annotation/eggnog"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. REQUIRED COLUMNS
# ============================================================

MAPPING_COLUMNS = [
    "fasta_header",
    "locus_tag_fasta",
    "gene_id",
    "gene_name",
    "gene_biotype",
    "locus_tag",
    "protein_id",
    "product",
    "dbxref",
]

EGGNOG_COLUMNS = [
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]

FUNCTIONAL_FIELDS = [
    "eggNOG_OGs",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]


# ============================================================
# 3. VALIDATE INPUT FILES
# ============================================================

print("\nChecking input files...")

if not ANNOTATION_FILE.exists():
    raise FileNotFoundError(
        f"eggNOG annotation file not found: {ANNOTATION_FILE}"
    )

if not MAPPING_FILE.exists():
    raise FileNotFoundError(
        f"Protein mapping file not found: {MAPPING_FILE}"
    )

print("Input files found.")


# ============================================================
# 4. LOAD PROTEIN MAPPING
# ============================================================

print("\nLoading protein-to-gene mapping...")

mapping = pd.read_csv(
    MAPPING_FILE,
    dtype=str
)

missing_mapping_columns = (
    set(MAPPING_COLUMNS) - set(mapping.columns)
)

if missing_mapping_columns:
    raise ValueError(
        f"Missing mapping columns: {missing_mapping_columns}"
    )

if mapping["fasta_header"].isna().any():
    raise ValueError(
        "Missing FASTA headers in protein mapping."
    )

if mapping["fasta_header"].duplicated().any():
    raise ValueError(
        "Duplicate FASTA headers in protein mapping."
    )

print(f"Proteins in mapping: {len(mapping)}")
print(
    "Unique gene IDs:",
    mapping["gene_id"].nunique()
)


# ============================================================
# 5. READ EGGNOG HEADER
# ============================================================

print("\nReading eggNOG annotation header...")

header_line = None

with open(ANNOTATION_FILE, "r", encoding="utf-8") as handle:

    for line in handle:

        if line.startswith("#query"):
            header_line = line
            break

if header_line is None:
    raise ValueError(
        "Could not find eggNOG header starting with #query."
    )

column_names = (
    header_line
    .lstrip("#")
    .strip()
    .split("\t")
)

if column_names != EGGNOG_COLUMNS:
    raise ValueError(
        "Unexpected eggNOG columns.\n"
        f"Observed: {column_names}\n"
        f"Expected: {EGGNOG_COLUMNS}"
    )

print(f"eggNOG columns identified: {len(column_names)}")


# ============================================================
# 6. LOAD EGGNOG ANNOTATIONS
# ============================================================

print("\nLoading eggNOG annotations...")

annotations = pd.read_csv(
    ANNOTATION_FILE,
    sep="\t",
    comment="#",
    names=column_names,
    header=None,
    dtype=str,
    keep_default_na=False
)

print(f"eggNOG annotation records: {len(annotations)}")

if annotations["query"].duplicated().any():
    raise ValueError(
        "Duplicate query IDs in eggNOG annotation file."
    )

if annotations["query"].isna().any():
    raise ValueError(
        "Missing query IDs in eggNOG annotations."
    )


# ============================================================
# 7. NORMALIZE EMPTY VALUES
# ============================================================

print("\nNormalizing missing annotation values...")

missing_values = ["", "-", "--", "NA", "NaN", "nan"]

for column in FUNCTIONAL_FIELDS:
    annotations[column] = annotations[column].replace(
        missing_values,
        pd.NA
    )


# ============================================================
# 8. MERGE MAPPING WITH EGGNOG
# ============================================================

print("\nMerging gene mapping with eggNOG annotations...")

merged = mapping.merge(
    annotations,
    left_on="fasta_header",
    right_on="query",
    how="left",
    validate="one_to_one",
    indicator=True
)

merged["annotation_record"] = (
    merged["_merge"] == "both"
)

merged.drop(
    columns=["_merge"],
    inplace=True
)

if len(merged) != len(mapping):
    raise ValueError(
        "The merge changed the number of protein records."
    )

print(f"Proteins after merge: {len(merged)}")


# ============================================================
# 9. IDENTIFY PROTEINS WITHOUT EGGNOG RECORDS
# ============================================================

unannotated = merged.loc[
    ~merged["annotation_record"]
].copy()

print(
    "Proteins without eggNOG record:",
    len(unannotated)
)


# ============================================================
# 10. CALCULATE ANNOTATION COVERAGE
# ============================================================

print("\nCalculating annotation coverage...")

total_proteins = len(merged)

annotated_records = int(
    merged["annotation_record"].sum()
)

summary = {
    "total_proteins": total_proteins,

    "proteins_with_eggnog_record": annotated_records,

    "proteins_without_eggnog_record": (
        total_proteins - annotated_records
    ),

    "eggnog_record_coverage_percent": round(
        annotated_records / total_proteins * 100,
        2
    ),
}


# ============================================================
# 11. FUNCTIONAL ANNOTATION COVERAGE
# ============================================================

for column in FUNCTIONAL_FIELDS:

    valid_values = merged[column].notna()

    count = int(valid_values.sum())

    percentage = round(
        count / total_proteins * 100,
        2
    )

    summary[f"{column}_count"] = count

    summary[f"{column}_coverage_percent"] = percentage


# ============================================================
# 12. ADDITIONAL ANNOTATION STATISTICS
# ============================================================

# Number of annotated proteins with at least one
# functional annotation field populated.

merged["has_functional_annotation"] = (
    merged[FUNCTIONAL_FIELDS]
    .notna()
    .any(axis=1)
)

proteins_with_function = int(
    merged["has_functional_annotation"].sum()
)

summary["proteins_with_at_least_one_functional_field"] = (
    proteins_with_function
)

summary["functional_annotation_coverage_percent"] = round(
    proteins_with_function / total_proteins * 100,
    2
)


# Number of eggNOG records without any functional field.

records_without_function = int(
    (
        merged["annotation_record"]
        & ~merged["has_functional_annotation"]
    ).sum()
)

summary["records_without_functional_fields"] = (
    records_without_function
)


# ============================================================
# 13. CHECK PROTEIN AND GENE IDENTIFIERS
# ============================================================

summary["unique_gene_ids"] = int(
    merged["gene_id"].nunique()
)

summary["unique_locus_tags"] = int(
    merged["locus_tag"].nunique()
)

summary["unique_protein_ids"] = int(
    merged["protein_id"].nunique()
)


# ============================================================
# 14. SAVE COMPLETE FUNCTIONAL ANNOTATION TABLE
# ============================================================

print("\nSaving complete annotation table...")

OUTPUT_TABLE = (
    OUTPUT_DIR
    / "candida_auris_gene_functional_annotation.csv"
)

merged.to_csv(
    OUTPUT_TABLE,
    index=False
)


# ============================================================
# 15. SAVE SUMMARY TABLE
# ============================================================

print("Saving annotation summary...")

summary_df = pd.DataFrame(
    [
        {
            "feature": key,
            "value": value
        }
        for key, value in summary.items()
    ]
)

SUMMARY_TSV = (
    OUTPUT_DIR
    / "eggnog_annotation_summary.tsv"
)

summary_df.to_csv(
    SUMMARY_TSV,
    sep="\t",
    index=False
)


# ============================================================
# 16. SAVE SUMMARY JSON
# ============================================================

SUMMARY_JSON = (
    OUTPUT_DIR
    / "eggnog_annotation_summary.json"
)

with open(
    SUMMARY_JSON,
    "w",
    encoding="utf-8"
) as handle:

    json.dump(
        summary,
        handle,
        indent=4,
        ensure_ascii=False
    )


# ============================================================
# 17. SAVE UNANNOTATED PROTEINS
# ============================================================

UNANNOTATED_CSV = (
    OUTPUT_DIR
    / "eggnog_unannotated_proteins.csv"
)

unannotated.to_csv(
    UNANNOTATED_CSV,
    index=False
)


# ============================================================
# 18. SAVE FUNCTIONAL COVERAGE TABLE
# ============================================================

coverage_rows = []

for column in FUNCTIONAL_FIELDS:

    coverage_rows.append(
        {
            "annotation_type": column,

            "annotated_proteins": summary[
                f"{column}_count"
            ],

            "total_proteins": total_proteins,

            "coverage_percent": summary[
                f"{column}_coverage_percent"
            ],
        }
    )

coverage_df = pd.DataFrame(coverage_rows)

COVERAGE_TSV = (
    OUTPUT_DIR
    / "eggnog_functional_coverage.tsv"
)

coverage_df.to_csv(
    COVERAGE_TSV,
    sep="\t",
    index=False
)


# ============================================================
# 19. FINAL VALIDATION
# ============================================================

assert len(merged) == len(mapping)

assert merged["gene_id"].nunique() == len(mapping)

assert (
    int(merged["annotation_record"].sum())
    + len(unannotated)
    == total_proteins
)

assert np.isfinite(
    np.array(
        [
            summary["eggnog_record_coverage_percent"],
            summary["functional_annotation_coverage_percent"]
        ],
        dtype=float
    )
).all()


# ============================================================
# 20. PRINT FINAL REPORT
# ============================================================

print("\n")
print("=" * 65)
print("CANDIDA AURIS - EGGNOG FUNCTIONAL ANNOTATION REPORT")
print("=" * 65)

print("\nGENERAL STATISTICS")
print("-" * 65)

print(f"Total proteins: {total_proteins}")

print(f"Unique genes: {summary['unique_gene_ids']}")

print(
    "Proteins with eggNOG records:",
    annotated_records
)

print(
    "Proteins without eggNOG records:",
    len(unannotated)
)

print(
    "eggNOG record coverage:",
    f"{summary['eggnog_record_coverage_percent']}%"
)

print(
    "Proteins with at least one functional field:",
    proteins_with_function
)

print(
    "Functional annotation coverage:",
    f"{summary['functional_annotation_coverage_percent']}%"
)


print("\nFUNCTIONAL ANNOTATION COVERAGE")
print("-" * 65)

for _, row in coverage_df.iterrows():

    print(
        f"{row['annotation_type']:<20}"
        f"{row['annotated_proteins']:>8}"
        f"{row['coverage_percent']:>10.2f}%"
    )


print("\nOUTPUT FILES")
print("-" * 65)

for path in [
    OUTPUT_TABLE,
    SUMMARY_TSV,
    SUMMARY_JSON,
    UNANNOTATED_CSV,
    COVERAGE_TSV
]:

    print(path)


print("\nValidation: PASSED")
print("=" * 65)
