#!/usr/bin/env python3

"""
Integrate edgeR and limma-voom differential expression results
with eggNOG functional annotations for Candida auris.

Genes not tested in differential expression are classified
separately as 'not_tested'.
"""

from pathlib import Path
import json

import numpy as np
import pandas as pd


# ============================================================
# 1. CONFIGURATION
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[1]

EDGER_FILE = (
    BASE_DIR / "results/differential_expression/edgeR_all_genes.csv"
)

LIMMA_FILE = (
    BASE_DIR / "results/differential_expression/limma_voom_all_genes.csv"
)

ANNOTATION_FILE = (
    BASE_DIR
    / "results/annotation/eggnog/"
    / "candida_auris_gene_functional_annotation.csv"
)

OUTPUT_DIR = BASE_DIR / "results/integrated_analysis"

FDR_THRESHOLD = 0.05

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

def check_file(path):
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")


def check_unique_gene_ids(df, name):
    if "gene_id" not in df.columns:
        raise ValueError(f"Missing gene_id column in {name}.")

    if df["gene_id"].isna().any():
        raise ValueError(f"Missing gene_id values in {name}.")

    if df["gene_id"].duplicated().any():
        raise ValueError(f"Duplicated gene_id values in {name}.")


def check_columns(df, columns, name):
    missing = [column for column in columns if column not in df.columns]

    if missing:
        raise ValueError(f"Missing columns in {name}: {missing}")


def convert_numeric(df, columns, name):
    for column in columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

        if df[column].isna().any():
            raise ValueError(
                f"Missing or non-numeric values in {name}.{column}"
            )


# ============================================================
# 3. CHECK INPUT FILES
# ============================================================

print("\n=== CHECKING INPUT FILES ===")

for path in [EDGER_FILE, LIMMA_FILE, ANNOTATION_FILE]:
    check_file(path)
    print(f"Found: {path}")


# ============================================================
# 4. LOAD DATA
# ============================================================

print("\n=== LOADING DATA ===")

edger = pd.read_csv(EDGER_FILE)
limma = pd.read_csv(LIMMA_FILE)
annotation = pd.read_csv(ANNOTATION_FILE)

print(f"edgeR genes: {len(edger)}")
print(f"limma-voom genes: {len(limma)}")
print(f"Annotated genes: {len(annotation)}")


# ============================================================
# 5. VALIDATE INPUT TABLES
# ============================================================

print("\n=== VALIDATING INPUT TABLES ===")

check_unique_gene_ids(edger, "edgeR")
check_unique_gene_ids(limma, "limma-voom")
check_unique_gene_ids(annotation, "annotation")

check_columns(
    edger,
    ["gene_id", "logFC", "logCPM", "F", "PValue", "FDR"],
    "edgeR"
)

check_columns(
    limma,
    ["gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"],
    "limma-voom"
)

check_columns(
    annotation,
    [
        "gene_id",
        "gene_name",
        "product",
        "Description",
        "Preferred_name",
        "GOs",
        "EC",
        "KEGG_ko",
        "KEGG_Pathway",
        "PFAMs",
        "annotation_record",
        "has_functional_annotation"
    ],
    "annotation"
)

convert_numeric(
    edger,
    ["logFC", "logCPM", "F", "PValue", "FDR"],
    "edgeR"
)

convert_numeric(
    limma,
    ["logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"],
    "limma-voom"
)

print("Input validation passed.")


# ============================================================
# 6. RENAME DE COLUMNS
# ============================================================

print("\n=== PREPARING DE TABLES ===")

edger = edger.rename(
    columns={
        "logFC": "edgeR_logFC",
        "logCPM": "edgeR_logCPM",
        "F": "edgeR_F",
        "PValue": "edgeR_PValue",
        "FDR": "edgeR_FDR"
    }
)

limma = limma.rename(
    columns={
        "logFC": "limma_logFC",
        "AveExpr": "limma_AveExpr",
        "t": "limma_t",
        "P.Value": "limma_PValue",
        "adj.P.Val": "limma_FDR",
        "B": "limma_B"
    }
)


# ============================================================
# 7. MERGE EDGE R AND LIMMA RESULTS
# ============================================================

de = pd.merge(
    edger,
    limma,
    on="gene_id",
    how="outer",
    indicator=True,
    validate="one_to_one"
)

print(f"Genes in merged DE table: {len(de)}")

if not de["_merge"].eq("both").all():
    raise ValueError(
        "edgeR and limma-voom contain different gene sets."
    )

de = de.drop(columns=["_merge"])

de["DE_results_available"] = True

# Explicit Boolean dtype.
de["edgeR_significant"] = (
    de["edgeR_FDR"].lt(FDR_THRESHOLD).astype(bool)
)

de["limma_significant"] = (
    de["limma_FDR"].lt(FDR_THRESHOLD).astype(bool)
)


# ============================================================
# 8. CLASSIFY DE RESULTS
# ============================================================

print("\n=== CLASSIFYING DIFFERENTIAL EXPRESSION ===")

de["DE_class"] = "not_significant"

common = (
    de["edgeR_significant"]
    & de["limma_significant"]
)

edger_only = (
    de["edgeR_significant"]
    & ~de["limma_significant"]
)

limma_only = (
    ~de["edgeR_significant"]
    & de["limma_significant"]
)

de.loc[common, "DE_class"] = "common_significant"
de.loc[edger_only, "DE_class"] = "edgeR_only"
de.loc[limma_only, "DE_class"] = "limma_only"

de["direction_concordant"] = (
    np.sign(de["edgeR_logFC"])
    == np.sign(de["limma_logFC"])
)

de["consensus_logFC"] = (
    de["edgeR_logFC"] + de["limma_logFC"]
) / 2

de["consensus_direction"] = np.select(
    [
        de["consensus_logFC"] > 0,
        de["consensus_logFC"] < 0
    ],
    [
        "up_in_resistant",
        "down_in_resistant"
    ],
    default="no_change"
)

de.to_csv(
    OUTPUT_DIR / "integrated_DE_methods.csv",
    index=False
)


# ============================================================
# 9. INTEGRATE WITH FUNCTIONAL ANNOTATION
# ============================================================

print("\n=== INTEGRATING FUNCTIONAL ANNOTATION ===")

integrated = pd.merge(
    annotation,
    de,
    on="gene_id",
    how="left",
    validate="one_to_one"
)

# Fill missing DE availability for genes absent from DE tables.
integrated["DE_results_available"] = (
    integrated["DE_results_available"]
    .fillna(False)
    .astype(bool)
)

# CRITICAL: convert significance columns to real Boolean dtype
# BEFORE applying logical negation (~).
integrated["edgeR_significant"] = (
    integrated["edgeR_significant"]
    .fillna(False)
    .astype(bool)
)

integrated["limma_significant"] = (
    integrated["limma_significant"]
    .fillna(False)
    .astype(bool)
)

integrated["direction_concordant"] = (
    integrated["direction_concordant"]
    .fillna(False)
    .astype(bool)
)

tested = integrated["DE_results_available"]

edgeR_sig = integrated["edgeR_significant"]
limma_sig = integrated["limma_significant"]

# Initialize all genes as not tested.
integrated["DE_class"] = "not_tested"

# Classify only genes tested by both methods.
integrated.loc[tested, "DE_class"] = "not_significant"

integrated.loc[
    tested & edgeR_sig & limma_sig,
    "DE_class"
] = "common_significant"

integrated.loc[
    tested & edgeR_sig & ~limma_sig,
    "DE_class"
] = "edgeR_only"

integrated.loc[
    tested & ~edgeR_sig & limma_sig,
    "DE_class"
] = "limma_only"


# ============================================================
# 10. SAVE INTEGRATED TABLE
# ============================================================

print("\n=== SAVING INTEGRATED TABLE ===")

integrated = integrated.sort_values(
    by=["DE_class", "gene_id"]
)

integrated.to_csv(
    OUTPUT_DIR / "integrated_DE_functional_annotation.csv",
    index=False
)

print("Integrated table saved.")


# ============================================================
# 11. EXPORT SUBSETS
# ============================================================

print("\n=== EXPORTING GENE SUBSETS ===")

subsets = {
    "common_significant_genes.csv": "common_significant",
    "edgeR_only_significant_genes.csv": "edgeR_only",
    "limma_only_significant_genes.csv": "limma_only",
    "non_significant_genes.csv": "not_significant",
    "genes_not_tested.csv": "not_tested"
}

for filename, category in subsets.items():

    subset = integrated[
        integrated["DE_class"] == category
    ]

    subset.to_csv(
        OUTPUT_DIR / filename,
        index=False
    )

    print(f"{filename}: {len(subset)} genes")


# ============================================================
# 12. SUMMARY STATISTICS
# ============================================================

print("\n=== GENERATING SUMMARY ===")

total_genes = len(integrated)

tested_count = int(tested.sum())
not_tested_count = int((~tested).sum())

common_count = int(
    (integrated["DE_class"] == "common_significant").sum()
)

edger_only_count = int(
    (integrated["DE_class"] == "edgeR_only").sum()
)

limma_only_count = int(
    (integrated["DE_class"] == "limma_only").sum()
)

not_significant_count = int(
    (integrated["DE_class"] == "not_significant").sum()
)

edger_significant_count = int(
    (tested & edgeR_sig).sum()
)

limma_significant_count = int(
    (tested & limma_sig).sum()
)

direction_concordant_count = int(
    (tested & integrated["direction_concordant"]).sum()
)

direction_discordant_count = (
    tested_count - direction_concordant_count
)

class_counts = {
    str(key): int(value)
    for key, value in integrated["DE_class"].value_counts().items()
}

summary = {
    "total_annotated_genes": int(len(annotation)),
    "edgeR_tested_genes": int(len(edger)),
    "limma_voom_tested_genes": int(len(limma)),
    "genes_tested_by_both_methods": tested_count,
    "genes_not_tested": not_tested_count,
    "FDR_threshold": FDR_THRESHOLD,
    "edgeR_significant_genes": edger_significant_count,
    "limma_voom_significant_genes": limma_significant_count,
    "common_significant_genes": common_count,
    "edgeR_only_significant_genes": edger_only_count,
    "limma_only_significant_genes": limma_only_count,
    "not_significant_in_either_method": not_significant_count,
    "direction_concordant_genes": direction_concordant_count,
    "direction_discordant_genes": direction_discordant_count,
    "DE_class_counts": class_counts
}

summary_rows = [
    {"metric": key, "value": value}
    for key, value in summary.items()
    if key != "DE_class_counts"
]

pd.DataFrame(summary_rows).to_csv(
    OUTPUT_DIR / "integration_summary.tsv",
    sep="\t",
    index=False
)

with open(
    OUTPUT_DIR / "integration_summary.json",
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
# 13. FINAL VALIDATION
# ============================================================

print("\n=== FINAL VALIDATION ===")

assert integrated["gene_id"].is_unique
assert integrated["DE_class"].notna().all()

assert total_genes == len(annotation)

assert tested_count == len(de)

assert common_count == 1210
assert edger_only_count == 115
assert limma_only_count == 4
assert not_significant_count == 4047
assert not_tested_count == 43

assert edger_significant_count == (
    common_count + edger_only_count
)

assert limma_significant_count == (
    common_count + limma_only_count
)

assert tested_count + not_tested_count == total_genes

assert (
    integrated["DE_class"].value_counts().sum()
    == total_genes
)

print("All validation checks passed.")


# ============================================================
# 14. FINAL REPORT
# ============================================================

print("\n========================================")
print("INTEGRATION COMPLETED SUCCESSFULLY")
print("========================================")

print(f"Total annotated genes: {total_genes}")
print(f"Genes tested by both methods: {tested_count}")
print(f"Genes not tested: {not_tested_count}")

print("\nDifferential expression classification:")

for category in [
    "common_significant",
    "edgeR_only",
    "limma_only",
    "not_significant",
    "not_tested"
]:
    count = int(
        (integrated["DE_class"] == category).sum()
    )

    print(f"{category}: {count}")

print("\nDirection concordance:")
print(f"Concordant: {direction_concordant_count}")
print(f"Discordant: {direction_discordant_count}")

print("\nOutput directory:")
print(OUTPUT_DIR)

print("\nAll files were generated and validated.")

