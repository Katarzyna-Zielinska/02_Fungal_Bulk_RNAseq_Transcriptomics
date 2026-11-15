import pandas as pd
import re
from pathlib import Path

project = Path(
    "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"
)

fasta = project / "data/reference/proteins/candida_auris_proteins.faa"
annotation_file = project / "data/annotation/candida_auris_gene_annotation.csv"

annotation = pd.read_csv(annotation_file)

fasta_records = []

with fasta.open() as handle:
    for line in handle:
        if line.startswith(">"):
            header = line[1:].strip()

            match = re.search(r"(B9J08_\d+)", header)

            fasta_records.append({
                "fasta_header": header,
                "locus_tag_fasta": match.group(1) if match else None
            })

proteins = pd.DataFrame(fasta_records)

proteins["gene_id"] = (
    "gene-" + proteins["locus_tag_fasta"].fillna("")
)

merged = proteins.merge(
    annotation,
    on="gene_id",
    how="left",
    suffixes=("_fasta", "_annotation")
)

mapped = merged["locus_tag"].notna()

print("=== PROTEIN-GENE MAPPING VALIDATION ===")
print("FASTA sequences:", len(proteins))
print("Unique gene IDs:", proteins["gene_id"].nunique())
print("Successfully mapped:", mapped.sum())
print("Unmapped:", (~mapped).sum())
print("Duplicated gene IDs:", proteins["gene_id"].duplicated().sum())

merged.to_csv(
    project / "data/annotation/protein_gene_mapping.csv",
    index=False
)

print("\nSaved: data/annotation/protein_gene_mapping.csv")
