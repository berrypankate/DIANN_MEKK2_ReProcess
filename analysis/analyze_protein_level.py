import pandas as pd

INPUT_FILE = "../MEKK2_final.parquet"

MAP3K2_ID = "Q9Y2U5"
MAP3K3_ID = "Q99759"

df = pd.read_parquet(INPUT_FILE)

# Keep rows whose inferred protein group is MAP3K2 or MAP3K3
x = df[
    df["Protein.Group"].isin([MAP3K2_ID, MAP3K3_ID])
].copy()

# Keep confident protein-group assignments
x = x[
    x["PG.Q.Value"] <= 0.01
].copy()

# Assign treatment condition
def get_condition(run):
    if "_G_B24_" in run:
        return "DMSO"
    elif "_H_B24_" in run:
        return "DTAGv1"
    return "Unknown"

x["Condition"] = x["Run"].apply(get_condition)

# Give readable protein names
protein_names = {
    MAP3K2_ID: "MAP3K2",
    MAP3K3_ID: "MAP3K3"
}

x["Protein"] = x["Protein.Group"].map(protein_names)

# DIA-NN repeats protein-level quantities on peptide rows,
# so keep one value per protein per run.
protein_run = (
    x[
        [
            "Run",
            "Condition",
            "Protein",
            "Protein.Group",
            "PG.MaxLFQ",
            "Genes.MaxLFQ",
            "Genes.MaxLFQ.Unique",
            "PG.Q.Value"
        ]
    ]
    .drop_duplicates()
    .sort_values(["Protein", "Condition", "Run"])
)

print()
print("DIA-NN PROTEIN-LEVEL VALUES")
print("=" * 110)

print(
    protein_run.to_string(index=False)
)

print()
print("MEDIAN PROTEIN ABUNDANCE BY CONDITION")
print("=" * 110)

summary = (
    protein_run
    .groupby(["Protein", "Condition"])
    .agg(
        Runs=("Run", "nunique"),
        Median_PG_MaxLFQ=("PG.MaxLFQ", "median"),
        Median_Gene_MaxLFQ=("Genes.MaxLFQ", "median"),
        Median_Gene_MaxLFQ_Unique=("Genes.MaxLFQ.Unique", "median")
    )
)

print(summary)

protein_run.to_csv(
    "map3k2_map3k3_protein_level.csv",
    index=False
)

print()
print("Saved: map3k2_map3k3_protein_level.csv")