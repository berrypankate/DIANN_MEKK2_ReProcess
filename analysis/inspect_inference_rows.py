import pandas as pd
from pathlib import Path

# ============================================================
# PATHS
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

FILES = {
    "Baseline": PROJECT_DIR / "MEKK2_final.parquet",
    "PF": PROJECT_DIR / "results/MEKK_PF/final/MEKK_PF_final.parquet",
    "NPI": PROJECT_DIR / "results/MEKK2_NPI/final/MEKK2_NPI_final.parquet",
}

MAP3K2 = "Q9Y2U5"
MAP3K3 = "Q99759"

# Shared peptide we already identified
SHARED_PEPTIDE = "QVQFDPDSPETSK"


# ============================================================
# HELPERS
# ============================================================

def condition_from_run(run):
    run = str(run)

    if "_G_B24_" in run:
        return "DMSO"

    if "_H_B24_" in run:
        return "DTAGv1"

    return "Unknown"


def classify(protein_ids):
    proteins = set(str(protein_ids).split(";"))

    has_2 = MAP3K2 in proteins
    has_3 = MAP3K3 in proteins

    if has_2 and has_3:
        return "Shared MAP3K2/MAP3K3"

    if has_2:
        return "MAP3K2-specific"

    if has_3:
        return "MAP3K3-specific"

    return "Other"


# ============================================================
# LOAD ROWS FROM ALL THREE MODES
# ============================================================

all_target_rows = []

for mode, path in FILES.items():

    print()
    print("=" * 100)
    print(f"LOADING {mode}")
    print(path)
    print("=" * 100)

    df = pd.read_parquet(path)

    print(f"Total rows: {len(df):,}")

    # Rows mentioning MAP3K2 or MAP3K3 anywhere in Protein.Ids
    target = df[
        df["Protein.Ids"]
        .fillna("")
        .apply(
            lambda x:
                MAP3K2 in str(x).split(";")
                or MAP3K3 in str(x).split(";")
        )
    ].copy()

    target["Mode"] = mode
    target["Condition"] = target["Run"].apply(condition_from_run)
    target["Peptide_Class"] = target["Protein.Ids"].apply(classify)

    all_target_rows.append(target)

    print(f"MAP3K2/MAP3K3-related rows: {len(target):,}")

rows = pd.concat(all_target_rows, ignore_index=True)


# ============================================================
# 1. SHARED PEPTIDE — EXACT ROWS
# ============================================================

shared = rows[
    rows["Stripped.Sequence"] == SHARED_PEPTIDE
].copy()

shared_columns = [
    "Mode",
    "Run",
    "Condition",
    "Stripped.Sequence",
    "Protein.Ids",
    "Protein.Group",
    "Protein.Names",
    "Genes",
    "Proteotypic",
    "Precursor.Quantity",
    "PG.MaxLFQ",
    "Genes.MaxLFQ",
    "Genes.MaxLFQ.Unique",
    "Q.Value",
    "PG.Q.Value",
]

shared_columns = [
    c for c in shared_columns
    if c in shared.columns
]

print()
print("=" * 120)
print("SHARED PEPTIDE: EXACT DIA-NN ROWS")
print("=" * 120)

print(
    shared[shared_columns]
    .sort_values(["Mode", "Run"])
    .to_string(index=False)
)

shared[shared_columns].to_csv(
    SCRIPT_DIR / "shared_peptide_rows_by_mode.csv",
    index=False
)


# ============================================================
# 2. DOES PROTEIN.GROUP CHANGE BETWEEN MODES?
# ============================================================

grouping = (
    rows[
        [
            "Mode",
            "Stripped.Sequence",
            "Protein.Ids",
            "Protein.Group",
            "Proteotypic",
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            "Stripped.Sequence",
            "Mode",
            "Protein.Group"
        ]
    )
)

print()
print("=" * 120)
print("UNIQUE PEPTIDE -> PROTEIN MAPPING BY MODE")
print("=" * 120)

print(grouping.to_string(index=False))

grouping.to_csv(
    SCRIPT_DIR / "peptide_protein_mapping_by_mode.csv",
    index=False
)


# ============================================================
# 3. FIND PEPTIDES WHOSE GROUPING CHANGES
# ============================================================

group_sets = (
    grouping
    .groupby(["Stripped.Sequence", "Mode"])["Protein.Group"]
    .apply(
        lambda x: ";".join(
            sorted(set(map(str, x.dropna())))
        )
    )
    .unstack("Mode")
)

group_sets["Grouping_Changed"] = (
    group_sets.nunique(axis=1, dropna=False) > 1
)

changed = group_sets[
    group_sets["Grouping_Changed"]
].copy()

print()
print("=" * 120)
print("PEPTIDES WHOSE PROTEIN.GROUP CHANGES BETWEEN MODES")
print("=" * 120)

if changed.empty:
    print("No Protein.Group differences found.")
else:
    print(changed.to_string())

changed.to_csv(
    SCRIPT_DIR / "peptides_with_changed_grouping.csv"
)


# ============================================================
# 4. MAP3K2 / MAP3K3 PROTEIN-GROUP ROW COUNTS
# ============================================================

protein_group_counts = (
    rows[
        rows["Protein.Group"].isin([MAP3K2, MAP3K3])
    ]
    .groupby(
        [
            "Mode",
            "Condition",
            "Protein.Group"
        ]
    )
    .size()
    .reset_index(name="Row_Count")
)

print()
print("=" * 120)
print("PROTEIN-GROUP ROW COUNTS")
print("=" * 120)

print(protein_group_counts.to_string(index=False))

protein_group_counts.to_csv(
    SCRIPT_DIR / "protein_group_row_counts.csv",
    index=False
)


# ============================================================
# 5. PG.MAXLFQ — EXACT RUN-LEVEL VALUES
# ============================================================

pg_rows = rows[
    rows["Protein.Group"].isin([MAP3K2, MAP3K3])
].copy()

pg_run = (
    pg_rows
    .groupby(
        [
            "Mode",
            "Condition",
            "Run",
            "Protein.Group"
        ],
        as_index=False
    )
    .agg(
        PG_MaxLFQ_First=("PG.MaxLFQ", "first"),
        PG_MaxLFQ_Min=("PG.MaxLFQ", "min"),
        PG_MaxLFQ_Max=("PG.MaxLFQ", "max"),
        Underlying_Rows=("PG.MaxLFQ", "size"),
        Unique_Peptides=("Stripped.Sequence", "nunique"),
    )
)

print()
print("=" * 120)
print("RUN-LEVEL PG.MaxLFQ VALUES")
print("=" * 120)

print(pg_run.to_string(index=False))

pg_run.to_csv(
    SCRIPT_DIR / "pg_maxlfq_run_level_by_mode.csv",
    index=False
)


# ============================================================
# 6. MEDIAN PG.MAXLFQ BY MODE / CONDITION
# ============================================================

pg_median = (
    pg_run
    .groupby(
        [
            "Mode",
            "Condition",
            "Protein.Group"
        ],
        as_index=False
    )
    .agg(
        Median_PG_MaxLFQ=("PG_MaxLFQ_First", "median"),
        Median_Unique_Peptides=("Unique_Peptides", "median"),
    )
)

print()
print("=" * 120)
print("MEDIAN PG.MaxLFQ BY MODE")
print("=" * 120)

print(pg_median.to_string(index=False))

pg_median.to_csv(
    SCRIPT_DIR / "pg_maxlfq_medians_by_mode.csv",
    index=False
)


# ============================================================
# 7. EXACT PEPTIDE SET DIFFERENCES
# ============================================================

peptide_sets = {}

for mode in FILES:

    subset = rows[
        (rows["Mode"] == mode)
        &
        (rows["Q.Value"] <= 0.01)
        &
        (rows["PG.Q.Value"] <= 0.01)
    ]

    peptide_sets[mode] = set(
        subset["Stripped.Sequence"].dropna()
    )


baseline = peptide_sets["Baseline"]
pf = peptide_sets["PF"]
npi = peptide_sets["NPI"]

print()
print("=" * 120)
print("PEPTIDE SET DIFFERENCES")
print("=" * 120)

print("PF only vs Baseline:")
print(sorted(pf - baseline))

print()
print("Baseline only vs PF:")
print(sorted(baseline - pf))

print()
print("NPI only vs Baseline:")
print(sorted(npi - baseline))

print()
print("Baseline only vs NPI:")
print(sorted(baseline - npi))


# ============================================================
# COMPLETE
# ============================================================

print()
print("=" * 120)
print("SAVED OUTPUTS")
print("=" * 120)

for filename in [
    "shared_peptide_rows_by_mode.csv",
    "peptide_protein_mapping_by_mode.csv",
    "peptides_with_changed_grouping.csv",
    "protein_group_row_counts.csv",
    "pg_maxlfq_run_level_by_mode.csv",
    "pg_maxlfq_medians_by_mode.csv",
]:
    print(SCRIPT_DIR / filename)