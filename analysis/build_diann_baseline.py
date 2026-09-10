import pandas as pd

INPUT_FILE = "../MEKK2_final.parquet"

MAP3K2_ID = "Q9Y2U5"
MAP3K3_ID = "Q99759"

df = pd.read_parquet(INPUT_FILE)

# -----------------------------
# Condition mapping
# -----------------------------

def get_condition(run):
    if "_G_B24_" in run:
        return "DMSO"
    elif "_H_B24_" in run:
        return "DTAGv1"
    return "Unknown"


# -----------------------------
# Peptide-level evidence
# -----------------------------

def classify_peptide(protein_ids):
    proteins = set(str(protein_ids).split(";"))

    has_map3k2 = MAP3K2_ID in proteins
    has_map3k3 = MAP3K3_ID in proteins

    if has_map3k2 and has_map3k3:
        return "Shared"
    elif has_map3k2:
        return "MAP3K2-specific"
    elif has_map3k3:
        return "MAP3K3-specific"

    return "Other"


peptides = df[
    df["Protein.Ids"].fillna("").apply(
        lambda s: (
            MAP3K2_ID in str(s).split(";")
            or MAP3K3_ID in str(s).split(";")
        )
    )
].copy()

# Confidence filter
peptides = peptides[
    (peptides["Q.Value"] <= 0.01)
    & (peptides["PG.Q.Value"] <= 0.01)
].copy()

peptides["Condition"] = peptides["Run"].apply(get_condition)
peptides["Peptide_Class"] = peptides["Protein.Ids"].apply(classify_peptide)

# Count unique peptide sequences per run
peptide_counts = (
    peptides
    .groupby(["Run", "Condition", "Peptide_Class"])["Stripped.Sequence"]
    .nunique()
    .unstack(fill_value=0)
    .reset_index()
)

for col in ["MAP3K2-specific", "MAP3K3-specific", "Shared"]:
    if col not in peptide_counts.columns:
        peptide_counts[col] = 0


# -----------------------------
# Protein-level evidence
# -----------------------------

proteins = df[
    df["Protein.Group"].isin([MAP3K2_ID, MAP3K3_ID])
    & (df["PG.Q.Value"] <= 0.01)
].copy()

proteins["Condition"] = proteins["Run"].apply(get_condition)

# PG.MaxLFQ is repeated on peptide rows, so take one value per protein/run
protein_run = (
    proteins
    .groupby(
        ["Run", "Condition", "Protein.Group"],
        as_index=False
    )
    .agg(
        PG_MaxLFQ=("PG.MaxLFQ", "first")
    )
)

protein_wide = protein_run.pivot_table(
    index=["Run", "Condition"],
    columns="Protein.Group",
    values="PG_MaxLFQ"
).reset_index()

protein_wide = protein_wide.rename(
    columns={
        MAP3K2_ID: "MAP3K2_PG_MaxLFQ",
        MAP3K3_ID: "MAP3K3_PG_MaxLFQ"
    }
)

if "MAP3K2_PG_MaxLFQ" not in protein_wide.columns:
    protein_wide["MAP3K2_PG_MaxLFQ"] = pd.NA

if "MAP3K3_PG_MaxLFQ" not in protein_wide.columns:
    protein_wide["MAP3K3_PG_MaxLFQ"] = pd.NA

protein_wide["MAP3K2_Reported"] = (
    protein_wide["MAP3K2_PG_MaxLFQ"].notna()
)

protein_wide["MAP3K3_Reported"] = (
    protein_wide["MAP3K3_PG_MaxLFQ"].notna()
)


# -----------------------------
# Merge into benchmark table
# -----------------------------

baseline = peptide_counts.merge(
    protein_wide,
    on=["Run", "Condition"],
    how="outer"
)

baseline = baseline[
    [
        "Run",
        "Condition",
        "MAP3K2-specific",
        "MAP3K3-specific",
        "Shared",
        "MAP3K2_Reported",
        "MAP3K3_Reported",
        "MAP3K2_PG_MaxLFQ",
        "MAP3K3_PG_MaxLFQ"
    ]
].sort_values(["Condition", "Run"])

print()
print("DIA-NN BASELINE")
print("=" * 120)
print(baseline.to_string(index=False))

baseline.to_csv(
    "diann_baseline_per_run.csv",
    index=False
)

print()
print("Saved: diann_baseline_per_run.csv")