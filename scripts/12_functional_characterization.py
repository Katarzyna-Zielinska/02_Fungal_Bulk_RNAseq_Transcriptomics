#!/usr/bin/env python3

"""
Functional characterization of common significant genes
in Candida auris resistant versus sensitive isolates.

Input:
    integrated_DE_functional_annotation.csv

Output:
    Functional characterization tables in:
    results/functional_characterization/
"""

from pathlib import Path
import json
import re

import pandas as pd


# ============================================================
# 1. CONFIGURATION
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    BASE_DIR
    / "results/integrated_analysis/"
    / "integrated_DE_functional_annotation.csv"
)

OUTPUT_DIR = BASE_DIR / "results/functional_characterization"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FDR_THRESHOLD = 0.05


# ============================================================
# 2. LOAD DATA
# ============================================================

print("\n=== LOADING INTEGRATED DATA ===")

if not INPUT_FILE.exists():
    raise FileNotFoundError(INPUT_FILE)

df = pd.read_csv(INPUT_FILE)

required_columns = [
    "gene_id",
    "gene_name",
    "product",
    "Description",
    "Preferred_name",
    "COG_category",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "PFAMs",
    "edgeR_logFC",
    "edgeR_FDR",
    "limma_logFC",
    "limma_FDR",
    "DE_class",
    "consensus_direction"
]

missing = [
    col for col in required_columns
    if col not in df.columns
]

if missing:
    raise ValueError(f"Missing columns: {missing}")

print("Loaded genes:", len(df))


# ============================================================
# 3. SELECT COMMON SIGNIFICANT GENES
# ============================================================

print("\n=== SELECTING COMMON SIGNIFICANT GENES ===")

common = df[
    df["DE_class"] == "common_significant"
].copy()

if len(common) != 1210:
    raise ValueError(
        f"Expected 1210 common significant genes, found {len(common)}"
    )

common["edgeR_FDR"] = pd.to_numeric(
    common["edgeR_FDR"],
    errors="coerce"
)

common["edgeR_logFC"] = pd.to_numeric(
    common["edgeR_logFC"],
    errors="coerce"
)

common["limma_FDR"] = pd.to_numeric(
    common["limma_FDR"],
    errors="coerce"
)

common["limma_logFC"] = pd.to_numeric(
    common["limma_logFC"],
    errors="coerce"
)

if common[
    ["edgeR_FDR", "edgeR_logFC", "limma_FDR", "limma_logFC"]
].isna().any().any():
    raise ValueError("Missing numerical DE values in common significant genes.")

common = common.sort_values("edgeR_FDR")

print("Common significant genes:", len(common))


# ============================================================
# 4. EXPORT COMPLETE COMMON GENE LIST
# ============================================================

print("\n=== EXPORTING COMMON GENE LIST ===")

common.to_csv(
    OUTPUT_DIR / "common_significant_functionally_annotated.csv",
    index=False
)


# ============================================================
# 5. SPLIT BY EXPRESSION DIRECTION
# ============================================================

print("\n=== SPLITTING BY EXPRESSION DIRECTION ===")

up = common[
    common["edgeR_logFC"] > 0
].copy()

down = common[
    common["edgeR_logFC"] < 0
].copy()

up.to_csv(
    OUTPUT_DIR / "genes_upregulated_in_resistant.csv",
    index=False
)

down.to_csv(
    OUTPUT_DIR / "genes_downregulated_in_resistant.csv",
    index=False
)

print("Upregulated in resistant:", len(up))
print("Downregulated in resistant:", len(down))


# ============================================================
# 6. FUNCTIONAL ANNOTATION COVERAGE
# ============================================================

print("\n=== FUNCTIONAL ANNOTATION COVERAGE ===")

annotation_columns = [
    "product",
    "Description",
    "Preferred_name",
    "COG_category",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "PFAMs"
]

coverage_rows = []

for column in annotation_columns:

    present = (
        common[column].notna()
        & common[column].astype(str).str.strip().ne("")
    )

    coverage_rows.append({
        "annotation_field": column,
        "genes_with_annotation": int(present.sum()),
        "total_genes": len(common),
        "coverage_percent": round(
            present.mean() * 100,
            2
        )
    })

coverage = pd.DataFrame(coverage_rows)

coverage.to_csv(
    OUTPUT_DIR / "functional_annotation_coverage.tsv",
    sep="\t",
    index=False
)

print(coverage.to_string(index=False))


# ============================================================
# 7. FUNCTIONAL DESCRIPTION SUMMARY
# ============================================================

print("\n=== FUNCTIONAL DESCRIPTION SUMMARY ===")

description_summary = (
    common["Description"]
    .dropna()
    .astype(str)
    .value_counts()
    .rename_axis("functional_description")
    .reset_index(name="gene_count")
)

description_summary.to_csv(
    OUTPUT_DIR / "functional_description_summary.csv",
    index=False
)


# ============================================================
# 8. COG CATEGORY SUMMARY
# ============================================================

print("\n=== COG CATEGORY SUMMARY ===")

cog = common[
    common["COG_category"].notna()
].copy()

cog["COG_category"] = cog["COG_category"].astype(str)

cog_rows = []

for _, row in cog.iterrows():

    categories = re.findall(
        r"[A-Z]",
        row["COG_category"]
    )

    for category in set(categories):

        cog_rows.append({
            "gene_id": row["gene_id"],
            "COG_category": category,
            "edgeR_logFC": row["edgeR_logFC"],
            "edgeR_FDR": row["edgeR_FDR"]
        })

cog_long = pd.DataFrame(cog_rows)

if not cog_long.empty:

    cog_summary = (
        cog_long.groupby("COG_category")
        .agg(
            gene_count=("gene_id", "nunique"),
            mean_logFC=("edgeR_logFC", "mean")
        )
        .reset_index()
        .sort_values("gene_count", ascending=False)
    )

    cog_summary.to_csv(
        OUTPUT_DIR / "COG_category_summary.csv",
        index=False
    )


# ============================================================
# 9. PFAM DOMAIN SUMMARY
# ============================================================

print("\n=== PFAM DOMAIN SUMMARY ===")

pfam_rows = []

for _, row in common.iterrows():

    if pd.isna(row["PFAMs"]):
        continue

    domains = [
        item.strip()
        for item in str(row["PFAMs"]).split(",")
        if item.strip()
    ]

    for domain in set(domains):

        pfam_rows.append({
            "gene_id": row["gene_id"],
            "PFAM_domain": domain,
            "edgeR_logFC": row["edgeR_logFC"],
            "edgeR_FDR": row["edgeR_FDR"]
        })

pfam_long = pd.DataFrame(pfam_rows)

if not pfam_long.empty:

    pfam_summary = (
        pfam_long.groupby("PFAM_domain")
        .agg(
            gene_count=("gene_id", "nunique"),
            mean_logFC=("edgeR_logFC", "mean")
        )
        .reset_index()
        .sort_values("gene_count", ascending=False)
    )

    pfam_summary.to_csv(
        OUTPUT_DIR / "PFAM_domain_summary.csv",
        index=False
    )


# ============================================================
# 10. KEGG PATHWAY MEMBERSHIP
# ============================================================

print("\n=== KEGG PATHWAY MEMBERSHIP ===")

pathway_rows = []

for _, row in common.iterrows():

    if pd.isna(row["KEGG_Pathway"]):
        continue

    pathways = [
        item.strip()
        for item in str(row["KEGG_Pathway"]).split(",")
        if item.strip()
    ]

    for pathway in set(pathways):

        pathway_rows.append({
            "gene_id": row["gene_id"],
            "KEGG_Pathway": pathway,
            "edgeR_logFC": row["edgeR_logFC"],
            "edgeR_FDR": row["edgeR_FDR"]
        })

pathway_long = pd.DataFrame(pathway_rows)

if not pathway_long.empty:

    pathway_summary = (
        pathway_long.groupby("KEGG_Pathway")
        .agg(
            gene_count=("gene_id", "nunique"),
            mean_logFC=("edgeR_logFC", "mean"),
            minimum_FDR=("edgeR_FDR", "min")
        )
        .reset_index()
        .sort_values("gene_count", ascending=False)
    )

    pathway_summary.to_csv(
        OUTPUT_DIR / "KEGG_pathway_membership.csv",
        index=False
    )

    pathway_long.to_csv(
        OUTPUT_DIR / "KEGG_pathway_gene_membership.csv",
        index=False
    )

    print("KEGG pathways:", pathway_summary.shape[0])

else:
    print("No KEGG pathway annotations found.")


# ============================================================
# 11. KEGG ORTHOLOGY SUMMARY
# ============================================================

print("\n=== KEGG ORTHOLOGY SUMMARY ===")

ko_rows = []

for _, row in common.iterrows():

    if pd.isna(row["KEGG_ko"]):
        continue

    kos = [
        item.strip()
        for item in str(row["KEGG_ko"]).split(",")
        if item.strip()
    ]

    for ko in set(kos):

        ko_rows.append({
            "gene_id": row["gene_id"],
            "KEGG_ko": ko,
            "edgeR_logFC": row["edgeR_logFC"],
            "edgeR_FDR": row["edgeR_FDR"]
        })

ko_long = pd.DataFrame(ko_rows)

if not ko_long.empty:

    ko_summary = (
        ko_long.groupby("KEGG_ko")
        .agg(
            gene_count=("gene_id", "nunique"),
            mean_logFC=("edgeR_logFC", "mean")
        )
        .reset_index()
        .sort_values("gene_count", ascending=False)
    )

    ko_summary.to_csv(
        OUTPUT_DIR / "KEGG_orthology_summary.csv",
        index=False
    )


# ============================================================
# 12. PRELIMINARY BIOLOGICAL CATEGORIES
# ============================================================

print("\n=== PRELIMINARY BIOLOGICAL CATEGORIES ===")

categories = {
    "transport_and_membrane": [
        "transporter",
        "transport",
        "permease",
        "major facilitator",
        "membrane",
        "efflux",
        "arn family"
    ],
    "iron_homeostasis": [
        "iron",
        "ferric",
        "ferrous",
        "siderophore",
        "ftr1"
    ],
    "oxidative_stress_and_redox": [
        "oxidoreductase",
        "redoxin",
        "oxidative",
        "glutathione",
        "peroxidase",
        "catalase",
        "superoxide",
        "flr1"
    ],
    "metabolism": [
        "metabolic",
        "kinase",
        "dehydrogenase",
        "metabolism",
        "biosynthesis",
        "synthase"
    ],
    "protein_homeostasis": [
        "ubiquitin",
        "proteasome",
        "chaperone",
        "heat shock",
        "protein folding"
    ],
    "cell_wall_and_membrane_biology": [
        "cell wall",
        "glucan",
        "chitin",
        "ergosterol",
        "sterol",
        "membrane"
    ]
}

# Combine available descriptive annotation fields.
common["search_text"] = (
    common[
        [
            "product",
            "Description",
            "Preferred_name",
            "PFAMs"
        ]
    ]
    .fillna("")
    .astype(str)
    .agg(" ".join, axis=1)
    .str.lower()
)

category_summary = []
category_gene_rows = []

for category, keywords in categories.items():

    mask = common["search_text"].apply(
        lambda text: any(
            keyword in text
            for keyword in keywords
        )
    )

    selected = common[mask].copy()

    category_summary.append({
        "category": category,
        "gene_count": len(selected),
        "genes_upregulated": int(
            (selected["edgeR_logFC"] > 0).sum()
        ),
        "genes_downregulated": int(
            (selected["edgeR_logFC"] < 0).sum()
        )
    })

    for _, row in selected.iterrows():

        category_gene_rows.append({
            "category": category,
            "gene_id": row["gene_id"],
            "gene_name": row["gene_name"],
            "product": row["product"],
            "Description": row["Description"],
            "Preferred_name": row["Preferred_name"],
            "edgeR_logFC": row["edgeR_logFC"],
            "edgeR_FDR": row["edgeR_FDR"],
            "limma_FDR": row["limma_FDR"]
        })

category_summary_df = pd.DataFrame(category_summary)

category_summary_df.to_csv(
    OUTPUT_DIR / "preliminary_biological_category_summary.csv",
    index=False
)

category_genes_df = pd.DataFrame(category_gene_rows)

category_genes_df.to_csv(
    OUTPUT_DIR / "preliminary_biological_category_genes.csv",
    index=False
)

print(category_summary_df.to_string(index=False))


# ============================================================
# 13. TOP CANDIDATE GENES
# ============================================================

print("\n=== TOP CANDIDATE GENES ===")

top_candidates = common.sort_values(
    "edgeR_FDR"
).head(50)

top_candidates.to_csv(
    OUTPUT_DIR / "top_50_common_significant_candidates.csv",
    index=False
)


# ============================================================
# 14. SUMMARY REPORT
# ============================================================

summary = {
    "total_integrated_genes": int(len(df)),
    "common_significant_genes": int(len(common)),
    "upregulated_in_resistant": int(len(up)),
    "downregulated_in_resistant": int(len(down)),
    "fdr_threshold": FDR_THRESHOLD,
    "genes_with_GO": int(
        common["GOs"].notna().sum()
    ),
    "genes_with_KEGG_ko": int(
        common["KEGG_ko"].notna().sum()
    ),
    "genes_with_KEGG_pathway": int(
        common["KEGG_Pathway"].notna().sum()
    ),
    "genes_with_PFAM": int(
        common["PFAMs"].notna().sum()
    ),
    "genes_with_COG": int(
        common["COG_category"].notna().sum()
    )
}

with open(
    OUTPUT_DIR / "functional_characterization_summary.json",
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
# 15. FINAL VALIDATION
# ============================================================

print("\n=== FINAL VALIDATION ===")

assert len(common) == 1210

assert len(up) + len(down) == len(common)

assert common["gene_id"].is_unique

assert not common["edgeR_FDR"].isna().any()

assert (
    common["edgeR_FDR"] < FDR_THRESHOLD
).all()

assert (
    common["limma_FDR"] < FDR_THRESHOLD
).all()

assert len(top_candidates) == 50

print("All validation checks passed.")

print("\n========================================")
print("FUNCTIONAL CHARACTERIZATION COMPLETED")
print("========================================")

print("Common significant genes:", len(common))
print("Upregulated in resistant:", len(up))
print("Downregulated in resistant:", len(down))
print("Output directory:", OUTPUT_DIR)

