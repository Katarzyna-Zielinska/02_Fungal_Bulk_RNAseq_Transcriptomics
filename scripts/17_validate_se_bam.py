import csv
import pysam
from pathlib import Path
from collections import defaultdict

BASE = Path(".")
RMATS = BASE / "results/rmats/alternative_splicing"
BAM_DIR = BASE / "results/star"

GENE = "gene-B9J08_004599"
EVENT_ID = "8"
CHROM = "PEKT02000009.1"

SAMPLES = [
    ("SRR10059687", "Sensitive", 0),
    ("SRR10059688", "Sensitive", 1),
    ("SRR10059689", "Sensitive", 2),
    ("SRR10059690", "Resistant", 3),
    ("SRR10059691", "Resistant", 4),
    ("SRR10059692", "Resistant", 5),
]

JUNCTIONS = {
    "upstream_to_SE": (320498, 320552),
    "SE_to_downstream": (320699, 320763),
    "exon_skipping": (320498, 320763),
}

def parse_values(value):
    return [x.strip() for x in value.split(",")]

def load_mats(mode):
    path = RMATS / f"SE.MATS.{mode}.txt"
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [
            row for row in reader
            if row["GeneID"].strip('"') == GENE
            and row["ID"] == EVENT_ID
        ]

    if len(rows) != 1:
        raise ValueError(
            f"{path}: oczekiwano 1 rekordu, znaleziono {len(rows)}"
        )
    return rows[0]

def load_individual_counts():
    path = RMATS / "individualCounts.SE.txt"

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [row for row in reader if row["ID"] == EVENT_ID]

    if len(rows) != 1:
        raise ValueError(
            f"individualCounts: oczekiwano 1 rekordu ID {EVENT_ID}, "
            f"znaleziono {len(rows)}"
        )

    return rows[0]

def audit_bam(sample):
    path = BAM_DIR / sample / f"{sample}_Aligned.sortedByCoord.out.bam"

    counts = defaultdict(lambda: {
        "all": 0,
        "NH1": 0,
        "names": set()
    })

    with pysam.AlignmentFile(path, "rb") as bam:
        if CHROM not in bam.references:
            raise ValueError(f"{sample}: brak kontigu {CHROM}")

        for read in bam.fetch(CHROM, 320400, 320850):
            if read.is_unmapped or read.cigartuples is None:
                continue

            ref_pos = read.reference_start

            for operation, length in read.cigartuples:
                if operation == 3:
                    junction = (ref_pos + 1, ref_pos + length)

                    if junction in JUNCTIONS.values():
                        item = counts[junction]
                        item["all"] += 1

                        if read.has_tag("NH") and read.get_tag("NH") == 1:
                            item["NH1"] += 1
                            item["names"].add(read.query_name)

                    ref_pos += length

                elif operation in (0, 2, 7, 8):
                    ref_pos += length

    result = {}

    for name, junction in JUNCTIONS.items():
        item = counts[junction]
        result[name] = {
            "BAM_all": item["all"],
            "BAM_NH1": item["NH1"],
            "BAM_unique_names": len(item["names"])
        }

    return result

jc = load_mats("JC")
jcec = load_mats("JCEC")
individual = load_individual_counts()

count_fields = [
    "upstream_to_target_count",
    "target_to_downstream_count",
    "target_count",
    "upstream_to_downstream_count"
]

parsed_counts = {
    field: parse_values(individual[field])
    for field in count_fields
}

inc_jc_1 = parse_values(jc["IncLevel1"])
inc_jc_2 = parse_values(jc["IncLevel2"])
inc_jcec_1 = parse_values(jcec["IncLevel1"])
inc_jcec_2 = parse_values(jcec["IncLevel2"])

rows = []

for sample, group, index in SAMPLES:
    bam = audit_bam(sample)

    up = bam["upstream_to_SE"]["BAM_unique_names"]
    down = bam["SE_to_downstream"]["BAM_unique_names"]
    skip = bam["exon_skipping"]["BAM_unique_names"]

    mean_inclusion = (up + down) / 2
    psi_bam = (
        mean_inclusion / (mean_inclusion + skip)
        if mean_inclusion + skip > 0 else None
    )

    group_index = index if group == "Sensitive" else index - 3

    rows.append({
        "sample": sample,
        "group": group,
        "gene_id": GENE,
        "event_id": EVENT_ID,
        "BAM_upstream_all": bam["upstream_to_SE"]["BAM_all"],
        "BAM_upstream_NH1": bam["upstream_to_SE"]["BAM_NH1"],
        "BAM_upstream_unique_names": up,
        "BAM_downstream_all": bam["SE_to_downstream"]["BAM_all"],
        "BAM_downstream_NH1": bam["SE_to_downstream"]["BAM_NH1"],
        "BAM_downstream_unique_names": down,
        "BAM_skipping_all": bam["exon_skipping"]["BAM_all"],
        "BAM_skipping_NH1": bam["exon_skipping"]["BAM_NH1"],
        "BAM_skipping_unique_names": skip,
        "rMATS_upstream_count": parsed_counts[count_fields[0]][index],
        "rMATS_downstream_count": parsed_counts[count_fields[1]][index],
        "rMATS_target_count": parsed_counts[count_fields[2]][index],
        "rMATS_skipping_count": parsed_counts[count_fields[3]][index],
        "rMATS_PSI_JC": (
            inc_jc_1[group_index] if group == "Sensitive"
            else inc_jc_2[group_index]
        ),
        "rMATS_PSI_JCEC": (
            inc_jcec_1[group_index] if group == "Sensitive"
            else inc_jcec_2[group_index]
        ),
        "BAM_PSI_proxy": round(psi_bam, 6) if psi_bam is not None else None,
        "BAM_PSI_method": "mean(unique upstream/downstream junction names) / (mean + unique skipping names)"
    })

output = BASE / "results/rmats/se_bam_validation_gene-B9J08_004599.csv"
output.parent.mkdir(parents=True, exist_ok=True)

with output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print("Zapisano:", output)
print("Liczba próbek:", len(rows))
print("Gen:", GENE)
print("Zdarzenie SE ID:", EVENT_ID)
print("JC FDR:", jc["FDR"])
print("JCEC FDR:", jcec["FDR"])
print("JC ΔPSI:", jc["IncLevelDifference"])
print("JCEC ΔPSI:", jcec["IncLevelDifference"])
print("\n========== WALIDACJA ZAKOŃCZONA ==========")
