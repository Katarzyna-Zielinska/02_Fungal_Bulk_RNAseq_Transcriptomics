#!/usr/bin/env python3

import csv
import json
from collections import Counter
from pathlib import Path


# ============================================================
# CONFIGURATION
# ============================================================

BASE_DIR = Path(".")

RMATS_DIR = BASE_DIR / "results/rmats/alternative_splicing"

ANNOTATION_FILE = (
    BASE_DIR
    / "results/integrated_analysis/"
    / "integrated_DE_functional_annotation.csv"
)

OUTPUT_DIR = BASE_DIR / "results/integrated_analysis/splicing"

EVENT_TYPES = ["SE", "MXE", "A5SS", "A3SS", "RI"]
EVIDENCE_TYPES = ["JC", "JCEC"]

FDR_THRESHOLD = 0.05


# ============================================================
# HELPERS
# ============================================================

def read_csv(path):
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    if not rows:
        print("No rows to write:", path)
        return

    fieldnames = list(dict.fromkeys(
        key
        for row in rows
        for key in row.keys()
    ))

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def write_tsv(path, rows, fieldnames):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def read_rmats_file(path):
    """
    Read rMATS files while preserving duplicated ID columns.

    The first ID is the event ID.
    The second ID is preserved as secondary_id.
    """

    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")

        header = next(reader)

        normalized_header = []
        id_count = 0

        for column in header:
            if column == "ID":
                id_count += 1

                if id_count == 1:
                    normalized_header.append("event_id")
                else:
                    normalized_header.append("secondary_id")
            else:
                normalized_header.append(column)

        rows = []

        for values in reader:
            if not values:
                continue

            if len(values) != len(normalized_header):
                raise ValueError(
                    "Column count mismatch in "
                    f"{path}: expected {len(normalized_header)}, "
                    f"received {len(values)}"
                )

            row = dict(zip(normalized_header, values))
            rows.append(row)

    return rows


def is_significant(row):
    value = row.get("FDR", "")

    if value in ("", "NA", "nan", "NaN"):
        return False

    return float(value) < FDR_THRESHOLD


def safe_float(value):
    if value in ("", "NA", "nan", "NaN", None):
        return None

    return float(value)


# ============================================================
# VALIDATE INPUTS
# ============================================================

print("=" * 70)
print("PROJECT 2 - INTEGRATION OF SPLICING AND GENE EXPRESSION")
print("=" * 70)

if not ANNOTATION_FILE.is_file():
    raise FileNotFoundError(ANNOTATION_FILE)

for event_type in EVENT_TYPES:
    for evidence_type in EVIDENCE_TYPES:
        path = (
            RMATS_DIR
            / f"{event_type}.MATS.{evidence_type}.txt"
        )

        if not path.is_file():
            raise FileNotFoundError(path)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# LOAD FUNCTIONAL ANNOTATION AND DIFFERENTIAL EXPRESSION
# ============================================================

print("\nLoading integrated gene annotation...")

annotation_rows = read_csv(ANNOTATION_FILE)

annotation_by_gene = {}

for row in annotation_rows:
    gene_id = row.get("gene_id", "")

    if not gene_id:
        continue

    if gene_id in annotation_by_gene:
        raise ValueError(
            f"Duplicated gene_id in annotation: {gene_id}"
        )

    annotation_by_gene[gene_id] = row

print("Annotation records:", len(annotation_rows))
print("Unique annotated gene IDs:", len(annotation_by_gene))


# ============================================================
# LOAD ALL RMATS RESULTS
# ============================================================

print("\nLoading rMATS results...")

all_events = []

event_counts = Counter()
significant_counts = Counter()

for event_type in EVENT_TYPES:
    for evidence_type in EVIDENCE_TYPES:

        filename = (
            f"{event_type}.MATS.{evidence_type}.txt"
        )

        path = RMATS_DIR / filename

        rows = read_rmats_file(path)

        print(
            f"{filename}: {len(rows)} events"
        )

        for row in rows:

            row["event_type"] = event_type
            row["evidence_type"] = evidence_type

            row["gene_id"] = row.get("GeneID", "")

            row["significant_FDR05"] = (
                is_significant(row)
            )

            row["FDR_threshold"] = FDR_THRESHOLD

            row["inc_level_difference_numeric"] = (
                safe_float(
                    row.get("IncLevelDifference", "")
                )
            )

            row["PValue_numeric"] = (
                safe_float(row.get("PValue", ""))
            )

            row["FDR_numeric"] = (
                safe_float(row.get("FDR", ""))
            )

            event_counts[
                (event_type, evidence_type)
            ] += 1

            if row["significant_FDR05"]:
                significant_counts[
                    (event_type, evidence_type)
                ] += 1

            annotation = annotation_by_gene.get(
                row["gene_id"]
            )

            if annotation is not None:

                row["annotation_match"] = True

                for key, value in annotation.items():

                    if key == "gene_id":
                        continue

                    row[f"annotation_{key}"] = value

            else:
                row["annotation_match"] = False

            all_events.append(row)


# ============================================================
# IDENTIFY GENES WITH SIGNIFICANT SPLICING
# ============================================================

significant_events = [
    row
    for row in all_events
    if row["significant_FDR05"]
]

significant_gene_ids = {
    row["gene_id"]
    for row in significant_events
    if row["gene_id"]
}


# ============================================================
# INTEGRATE EXPRESSION AND SPLICING
# ============================================================

print("\nIntegrating expression and splicing...")

for row in all_events:

    annotation = annotation_by_gene.get(
        row["gene_id"], {}
    )

    row["edgeR_significant"] = (
        annotation.get("edgeR_significant", "")
    )

    row["limma_significant"] = (
        annotation.get("limma_significant", "")
    )

    row["DE_class"] = (
        annotation.get("DE_class", "")
    )

    row["edgeR_logFC"] = (
        annotation.get("edgeR_logFC", "")
    )

    row["edgeR_FDR"] = (
        annotation.get("edgeR_FDR", "")
    )

    row["limma_logFC"] = (
        annotation.get("limma_logFC", "")
    )

    row["limma_FDR"] = (
        annotation.get("limma_FDR", "")
    )

    row["significant_splicing_gene"] = (
        row["gene_id"] in significant_gene_ids
    )

    de_significant = (
        str(row["edgeR_significant"]).lower() == "true"
        and str(row["limma_significant"]).lower() == "true"
    )

    row["common_significant_DE_and_splicing"] = (
        row["significant_FDR05"] and de_significant
    )


# ============================================================
# UNIQUE GENE-LEVEL SUMMARY
# ============================================================

gene_summary = []

for gene_id in sorted(significant_gene_ids):

    gene_events = [
        row
        for row in significant_events
        if row["gene_id"] == gene_id
    ]

    annotation = annotation_by_gene.get(
        gene_id, {}
    )

    event_types = sorted({
        row["event_type"]
        for row in gene_events
    })

    evidence_types = sorted({
        row["evidence_type"]
        for row in gene_events
    })

    jc_events = [
        row
        for row in gene_events
        if row["evidence_type"] == "JC"
    ]

    jcec_events = [
        row
        for row in gene_events
        if row["evidence_type"] == "JCEC"
    ]

    gene_summary.append({

        "gene_id": gene_id,

        "locus_tag": annotation.get(
            "locus_tag", ""
        ),

        "product": annotation.get(
            "product", ""
        ),

        "functional_description": annotation.get(
            "Description", ""
        ),

        "KEGG_ko": annotation.get(
            "KEGG_ko", ""
        ),

        "edgeR_logFC": annotation.get(
            "edgeR_logFC", ""
        ),

        "edgeR_FDR": annotation.get(
            "edgeR_FDR", ""
        ),

        "limma_logFC": annotation.get(
            "limma_logFC", ""
        ),

        "limma_FDR": annotation.get(
            "limma_FDR", ""
        ),

        "DE_class": annotation.get(
            "DE_class", ""
        ),

        "number_of_significant_event_records": len(
            gene_events
        ),

        "number_of_significant_JC_records": len(
            jc_events
        ),

        "number_of_significant_JCEC_records": len(
            jcec_events
        ),

        "event_types": ";".join(event_types),

        "evidence_types": ";".join(evidence_types),

        "common_significant_DE": (
            annotation.get("DE_class", "")
            == "common_significant"
        )
    })


# ============================================================
# SAVE OUTPUT TABLES
# ============================================================

print("\nSaving output tables...")

write_csv(
    OUTPUT_DIR / "rmats_all_events_with_annotation.csv",
    all_events
)

write_csv(
    OUTPUT_DIR / "rmats_significant_events_FDR05.csv",
    significant_events
)

write_csv(
    OUTPUT_DIR / "splicing_expression_integrated.csv",
    [
        row
        for row in significant_events
        if row["common_significant_DE_and_splicing"]
    ]
)

write_csv(
    OUTPUT_DIR / "splicing_gene_level_summary.csv",
    gene_summary
)


# ============================================================
# SUMMARY TABLES
# ============================================================

summary_rows = []

for event_type in EVENT_TYPES:

    for evidence_type in EVIDENCE_TYPES:

        total = event_counts[
            (event_type, evidence_type)
        ]

        significant = significant_counts[
            (event_type, evidence_type)
        ]

        summary_rows.append({

            "event_type": event_type,

            "evidence_type": evidence_type,

            "total_events": total,

            "significant_FDR05": significant,

            "not_significant": total - significant,

            "FDR_threshold": FDR_THRESHOLD
        })


write_tsv(
    OUTPUT_DIR / "rmats_event_summary.tsv",
    summary_rows,
    [
        "event_type",
        "evidence_type",
        "total_events",
        "significant_FDR05",
        "not_significant",
        "FDR_threshold"
    ]
)


# ============================================================
# VALIDATION
# ============================================================

print("\nValidation...")

print(
    "Total event records:",
    len(all_events)
)

print(
    "Significant event records:",
    len(significant_events)
)

print(
    "Unique genes with significant splicing:",
    len(significant_gene_ids)
)

print(
    "Genes with significant DE and splicing:",
    len({
        row["gene_id"]
        for row in significant_events
        if row["common_significant_DE_and_splicing"]
    })
)

print(
    "Significant records with annotation:",
    sum(
        row["annotation_match"]
        for row in significant_events
    )
)

print(
    "Significant records without annotation:",
    sum(
        not row["annotation_match"]
        for row in significant_events
    )
)


# ============================================================
# JSON SUMMARY
# ============================================================

json_summary = {

    "project": "Project 2 - Candida auris",

    "analysis": "rMATS and differential expression integration",

    "rmats_version": "4.4.0",

    "FDR_threshold": FDR_THRESHOLD,

    "event_types": EVENT_TYPES,

    "evidence_types": EVIDENCE_TYPES,

    "total_event_records": len(all_events),

    "significant_event_records": len(significant_events),

    "unique_genes_with_significant_splicing": len(
        significant_gene_ids
    ),

    "genes_with_significant_DE_and_splicing": len({
        row["gene_id"]
        for row in significant_events
        if row["common_significant_DE_and_splicing"]
    }),

    "event_summary": summary_rows,

    "unique_splicing_genes": sorted(
        significant_gene_ids
    )
}


with (
    OUTPUT_DIR / "splicing_expression_summary.json"
).open("w", encoding="utf-8") as handle:

    json.dump(
        json_summary,
        handle,
        indent=4,
        ensure_ascii=False
    )


# ============================================================
# FINISH
# ============================================================

print("\n" + "=" * 70)
print("INTEGRATION COMPLETED")
print("=" * 70)

print("Output directory:", OUTPUT_DIR)

for path in sorted(OUTPUT_DIR.iterdir()):

    if path.is_file():

        print(
            path.name,
            "|",
            path.stat().st_size,
            "bytes"
        )
