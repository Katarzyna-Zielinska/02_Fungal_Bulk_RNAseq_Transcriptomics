import csv
from collections import defaultdict
from pathlib import Path

BASE = Path("results/integrated_analysis/splicing")

INPUT = BASE / "splicing_expression_integrated.csv"
OUTPUT = BASE / "unique_splicing_events_summary.csv"

events = defaultdict(list)

with INPUT.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle):
        key = (
            row["gene_id"],
            row["event_type"],
            row["event_id"]
        )
        events[key].append(row)

output_rows = []

for (gene_id, event_type, event_id), records in sorted(events.items()):

    jc = next(
        (r for r in records if r["evidence_type"] == "JC"),
        None
    )

    jcec = next(
        (r for r in records if r["evidence_type"] == "JCEC"),
        None
    )

    reference = jc if jc is not None else jcec

    output_rows.append({
        "gene_id": gene_id,
        "event_type": event_type,
        "event_id": event_id,
        "JC_present": jc is not None,
        "JCEC_present": jcec is not None,
        "JC_FDR": jc["FDR"] if jc else "",
        "JCEC_FDR": jcec["FDR"] if jcec else "",
        "JC_delta_PSI": jc["IncLevelDifference"] if jc else "",
        "JCEC_delta_PSI": jcec["IncLevelDifference"] if jcec else "",
        "edgeR_logFC": reference["edgeR_logFC"],
        "edgeR_FDR": reference["edgeR_FDR"],
        "limma_logFC": reference["limma_logFC"],
        "limma_FDR": reference["limma_FDR"],
        "DE_class": reference["DE_class"],
        "annotation_product": reference["annotation_product"],
        "functional_description": reference["annotation_Description"],
        "annotation_match": reference["annotation_match"],
        "common_significant_DE_and_splicing":
            reference["common_significant_DE_and_splicing"],
        "JC_JCEC_complete": jc is not None and jcec is not None
    })

fields = list(output_rows[0].keys()) if output_rows else []

with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(output_rows)

print("Zapisano:", OUTPUT)
print("Liczba rekordów wejściowych:", sum(map(len, events.values())))
print("Liczba unikalnych zdarzeń:", len(output_rows))
print("Liczba genów:", len(set(r["gene_id"] for r in output_rows)))
print("Kompletne pary JC/JCEC:", sum(r["JC_JCEC_complete"] for r in output_rows))

print("\nPODSUMOWANIE:")
for row in output_rows:
    print(
        row["gene_id"],
        row["event_type"],
        row["event_id"],
        "| JC FDR:", row["JC_FDR"],
        "| JCEC FDR:", row["JCEC_FDR"],
        "| kompletna para:", row["JC_JCEC_complete"]
    )

print("\n========== DEDUPLICATION REPORT COMPLETE ==========")
