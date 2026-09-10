import pandas as pd

FILES = {
    "Baseline": "../MEKK2_final.parquet",
    "PF": "../results/MEKK_PF/final/MEKK_PF_final.parquet",
    "NPI": "../results/MEKK2_NPI/final/MEKK2_NPI_final.parquet",
}

MAP3K2 = "Q9Y2U5"
MAP3K3 = "Q99759"

def condition_from_run(run):
    if "_G_B24_" in run:
        return "DMSO"
    if "_H_B24_" in run:
        return "DTAGv1"
    return "Unknown"

def classify(ids):
    proteins = set(str(ids).split(";"))
    has2 = MAP3K2 in proteins
    has3 = MAP3K3 in proteins

    if has2 and has3:
        return "Shared"
    if has2:
        return "MAP3K2-specific"
    if has3:
        return "MAP3K3-specific"
    return "Other"

all_rows = []

for mode, path in FILES.items():
    df = pd.read_parquet(path)

    x = df[
        df["Protein.Ids"].fillna("").apply(
            lambda s: MAP3K2 in str(s).split(";") or MAP3K3 in str(s).split(";")
        )
    ].copy()

    x = x[
        (x["Q.Value"] <= 0.01) &
        (x["PG.Q.Value"] <= 0.01)
    ].copy()

    x["Condition"] = x["Run"].apply(condition_from_run)
    x["Peptide_Class"] = x["Protein.Ids"].apply(classify)

    counts = (
        x.groupby(["Condition", "Peptide_Class"])["Stripped.Sequence"]
        .nunique()
        .unstack(fill_value=0)
    )

    for condition in ["DMSO", "DTAGv1"]:
        row = {
            "Mode": mode,
            "Condition": condition,
            "MAP3K2-specific": counts.loc[condition].get("MAP3K2-specific", 0)
                if condition in counts.index else 0,
            "MAP3K3-specific": counts.loc[condition].get("MAP3K3-specific", 0)
                if condition in counts.index else 0,
            "Shared": counts.loc[condition].get("Shared", 0)
                if condition in counts.index else 0,
        }

        p2 = df[
            (df["Protein.Group"] == MAP3K2) &
            (df["PG.Q.Value"] <= 0.01) &
            (df["Run"].apply(condition_from_run) == condition)
        ]

        p3 = df[
            (df["Protein.Group"] == MAP3K3) &
            (df["PG.Q.Value"] <= 0.01) &
            (df["Run"].apply(condition_from_run) == condition)
        ]

        row["MAP3K2_Reported"] = not p2.empty
        row["MAP3K3_Reported"] = not p3.empty

        row["MAP3K2_Median_PG_MaxLFQ"] = (
            p2.groupby("Run")["PG.MaxLFQ"].first().median()
            if not p2.empty else None
        )

        row["MAP3K3_Median_PG_MaxLFQ"] = (
            p3.groupby("Run")["PG.MaxLFQ"].first().median()
            if not p3.empty else None
        )

        all_rows.append(row)

result = pd.DataFrame(all_rows)

print()
print("DIA-NN INFERENCE MODE COMPARISON")
print("=" * 110)
print(result.to_string(index=False))

result.to_csv("diann_inference_mode_comparison.csv", index=False)

print()
print("Saved: diann_inference_mode_comparison.csv")
