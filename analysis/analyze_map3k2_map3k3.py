import pandas as pd
import numpy as np

# ============================================================
# SETTINGS
# ============================================================

INPUT_FILE = "../MEKK2_final.parquet"

MAP3K2_ID = "Q9Y2U5"
MAP3K3_ID = "Q99759"

# Working assumption:
# DMSO = control
# DTAGv1 = MAP3K2 depletion
# This still needs supervisor confirmation.

# ============================================================
# LOAD DIA-NN RESULT
# ============================================================

df = pd.read_parquet(INPUT_FILE)

print(f"Total DIA-NN rows: {len(df):,}")

# ============================================================
# KEEP MAP3K2 / MAP3K3 EVIDENCE
# ============================================================

def contains_target(protein_ids):
    proteins = set(str(protein_ids).split(";"))

    return (
        MAP3K2_ID in proteins
        or MAP3K3_ID in proteins
    )


x = df[
    df["Protein.Ids"]
    .fillna("")
    .apply(contains_target)
].copy()

# ============================================================
# ASSIGN CONDITION
# ============================================================

x["Condition"] = np.where(
    x["Run"].str.contains("_G_B24_"),
    "DMSO",
    np.where(
        x["Run"].str.contains("_H_B24_"),
        "DTAGv1",
        "Unknown"
    )
)

# ============================================================
# CLASSIFY PEPTIDES
# ============================================================

def classify_peptide(protein_ids):

    proteins = set(str(protein_ids).split(";"))

    has_map3k2 = MAP3K2_ID in proteins
    has_map3k3 = MAP3K3_ID in proteins

    if has_map3k2 and has_map3k3:
        return "Shared MAP3K2/MAP3K3"

    if has_map3k2:
        return "MAP3K2-specific"

    if has_map3k3:
        return "MAP3K3-specific"

    return "Other"


x["Peptide_Class"] = (
    x["Protein.Ids"]
    .apply(classify_peptide)
)

# ============================================================
# CONFIDENCE FILTERING
# ============================================================

# 1% precursor-level FDR
# AND
# 1% protein-group-level FDR

x_filtered = x[
    (x["Q.Value"] <= 0.01)
    &
    (x["PG.Q.Value"] <= 0.01)
].copy()

print(
    f"MAP3K2/MAP3K3 rows after confidence filtering: "
    f"{len(x_filtered):,}"
)

# ============================================================
# SUMMARIZE EACH PEPTIDE WITHIN EACH CONDITION
# ============================================================

summary = (
    x_filtered
    .groupby(
        [
            "Peptide_Class",
            "Stripped.Sequence",
            "Condition"
        ],
        as_index=False
    )
    .agg(
        Median_Quantity=(
            "Precursor.Quantity",
            "median"
        ),

        Detection_Count=(
            "Run",
            "nunique"
        ),

        Median_Q_Value=(
            "Q.Value",
            "median"
        ),

        Median_PG_Q_Value=(
            "PG.Q.Value",
            "median"
        )
    )
)

# ============================================================
# MAKE ONE ROW PER PEPTIDE
# ============================================================

quantity = summary.pivot_table(
    index=[
        "Peptide_Class",
        "Stripped.Sequence"
    ],
    columns="Condition",
    values="Median_Quantity"
)

detections = summary.pivot_table(
    index=[
        "Peptide_Class",
        "Stripped.Sequence"
    ],
    columns="Condition",
    values="Detection_Count"
)

result = quantity.copy()

# Rename abundance columns

if "DMSO" in result.columns:
    result = result.rename(
        columns={
            "DMSO": "DMSO_Median"
        }
    )

if "DTAGv1" in result.columns:
    result = result.rename(
        columns={
            "DTAGv1": "DTAGv1_Median"
        }
    )

# Add detection counts

if "DMSO" in detections.columns:
    result["DMSO_Detections"] = (
        detections["DMSO"]
    )

if "DTAGv1" in detections.columns:
    result["DTAGv1_Detections"] = (
        detections["DTAGv1"]
    )

result = result.reset_index()

# Missing detection = zero detections

for col in [
    "DMSO_Detections",
    "DTAGv1_Detections"
]:

    if col not in result.columns:
        result[col] = 0

    result[col] = (
        result[col]
        .fillna(0)
        .astype(int)
    )

# ============================================================
# FOLD CHANGE
# ============================================================

result["DTAG_vs_DMSO"] = (
    result["DTAGv1_Median"]
    /
    result["DMSO_Median"]
)

result["log2FC"] = np.log2(
    result["DTAG_vs_DMSO"]
)

# ============================================================
# READABLE DETECTION RATES
# ============================================================

result["DMSO_Detection_Rate"] = (
    result["DMSO_Detections"]
    .astype(str)
    + "/3"
)

result["DTAGv1_Detection_Rate"] = (
    result["DTAGv1_Detections"]
    .astype(str)
    + "/3"
)

# ============================================================
# PRINT PEPTIDE TABLE
# ============================================================

display_columns = [
    "Peptide_Class",
    "Stripped.Sequence",
    "DMSO_Detection_Rate",
    "DTAGv1_Detection_Rate",
    "DMSO_Median",
    "DTAGv1_Median",
    "log2FC"
]

print()
print("MAP3K2 / MAP3K3 PEPTIDE SUMMARY")
print("=" * 100)

print(
    result[
        display_columns
    ].to_string(index=False)
)

# ============================================================
# CLASS SUMMARY
# ============================================================

class_summary = (
    result
    .groupby("Peptide_Class")
    .agg(
        Peptides=(
            "Stripped.Sequence",
            "nunique"
        ),

        DMSO_Total_Detections=(
            "DMSO_Detections",
            "sum"
        ),

        DTAGv1_Total_Detections=(
            "DTAGv1_Detections",
            "sum"
        )
    )
)

print()
print("CLASS SUMMARY")
print("=" * 100)
print(class_summary)

# ============================================================
# SAVE OUTPUT
# ============================================================

OUTPUT_FILE = (
    "map3k2_map3k3_peptide_summary.csv"
)

result.to_csv(
    OUTPUT_FILE,
    index=False
)

print()
print(f"Saved: {OUTPUT_FILE}")

# if "DMSO" in result.columns:
#     result = result.rename(columns={"DMSO": "DMSO_Median"})

# if "DTAGv1" in result.columns:
#     result = result.rename(columns={"DTAGv1": "DTAGv1_Median"})

# if "DMSO" in detections.columns:
#     result["DMSO_Detections"] = detections["DMSO"]

# if "DTAGv1" in detections.columns:
#     result["DTAGv1_Detections"] = detections["DTAGv1"]

# if "DMSO_Median" in result.columns and "DTAGv1_Median" in result.columns:
#     result["DTAG_vs_DMSO"] = (
#         result["DTAGv1_Median"] / result["DMSO_Median"]
#     )
#     result["log2FC"] = np.log2(result["DTAG_vs_DMSO"])

# result = result.reset_index()

# print("\nMAP3K2 / MAP3K3 PEPTIDE SUMMARY")
# print("=" * 80)
# print(result.to_string(index=False))

# print("\nCLASS COUNTS")
# print("=" * 80)

# class_counts = (
#     result.groupby("Peptide_Class")["Stripped.Sequence"]
#     .nunique()
#     .sort_values()
# )

# print(class_counts)

# # Save the result so we can use it later in VS Code / Excel / plotting
# OUTPUT_FILE = "map3k2_map3k3_peptide_summary.csv"
# result.to_csv(OUTPUT_FILE, index=False)

# print(f"\nSaved: {OUTPUT_FILE}")

# Replace missing detection counts with 0
