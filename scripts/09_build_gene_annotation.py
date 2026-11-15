import re
import pandas as pd
from pathlib import Path

project = Path(
    "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"
)

gff = Path(
    "/mnt/e/Fungal_Bioinformatics_Portfolio/01_Fungal_NGS_Pipeline/"
    "data/reference/ncbi_dataset/data/GCA_002759435.2/genomic.gff"
)

output_dir = project / "data" / "annotation"
output_dir.mkdir(parents=True, exist_ok=True)

def parse_attributes(attribute_string):
    attributes = {}
    for item in attribute_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes

genes = []
cds_records = []

with gff.open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue

        feature = fields[2]
        attributes = parse_attributes(fields[8])

        if feature == "gene":
            gene_id = attributes.get("ID", "")
            genes.append({
                "gene_id": gene_id,
                "gene_name": attributes.get("Name", ""),
                "gene_biotype": attributes.get("gene_biotype", ""),
                "locus_tag": attributes.get("locus_tag", "")
            })

        elif feature == "CDS":
            cds_records.append({
                "locus_tag": attributes.get("locus_tag", ""),
                "protein_id": attributes.get("protein_id", ""),
                "product": attributes.get("product", ""),
                "dbxref": attributes.get("Dbxref", "")
            })

gene_df = pd.DataFrame(genes)
cds_df = pd.DataFrame(cds_records)

# Remove duplicate CDS annotations for the same locus.
# Retain all distinct descriptions in a semicolon-separated field.
cds_df = (
    cds_df.groupby("locus_tag", as_index=False)
    .agg({
        "protein_id": lambda x: ";".join(sorted(set(v for v in x if v))),
        "product": lambda x: ";".join(sorted(set(v for v in x if v))),
        "dbxref": lambda x: ";".join(sorted(set(v for v in x if v)))
    })
)

annotation = gene_df.merge(
    cds_df,
    on="locus_tag",
    how="left"
)

annotation["product"] = annotation["product"].fillna("")

annotation.to_csv(
    output_dir / "candida_auris_gene_annotation.csv",
    index=False
)

print("=== GENE ANNOTATION SUMMARY ===")
print("Genes:", len(gene_df))
print("Genes with CDS product:", (annotation["product"] != "").sum())
print(
    "Hypothetical proteins:",
    annotation["product"].str.contains(
        "hypothetical protein", case=False
    ).sum()
)
print(
    "Genes with other product descriptions:",
    (
        (annotation["product"] != "") &
        (~annotation["product"].str.contains(
            "hypothetical protein", case=False
        ))
    ).sum()
)
print("Saved:", output_dir / "candida_auris_gene_annotation.csv")
print("\nExample annotations:")
print(annotation.head(10).to_string(index=False))
