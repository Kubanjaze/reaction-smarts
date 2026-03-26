import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

FAMILY_COLORS = {"benz": "#4C72B0", "naph": "#DD8452", "ind": "#55A868",
                 "quin": "#C44E52", "pyr": "#8172B2", "bzim": "#937860", "other": "#808080"}

# Reaction SMARTS library — single-reactant form: test if compound contains the reactant pattern
REACTIONS = {
    "amide_coupling":       ("[C:1](=O)[OH]",      "acid → amide"),
    "ester_formation":      ("[C:1](=O)[OH]",       "acid → ester"),
    "reductive_amination":  ("[CH1,CH2:1]=O",       "aldehyde/ketone → amine"),
    "sulfonamide_form":     ("[S:1](=O)(=O)[Cl]",   "sulfonyl Cl → sulfonamide"),
    "buchwald_amine":       ("[c:1][Br,Cl]",        "aryl halide → aryl amine (Pd)"),
    "suzuki_coupling":      ("[c:1][Br,Cl]",        "aryl halide → biaryl (Pd)"),
}

# Full reaction SMARTS for product validation
REACTION_SMARTS = {
    "amide_coupling":       "[C:1](=O)[OH:2].[NH2:3]>>[C:1](=O)[N:3]",
    "ester_formation":      "[C:1](=O)[OH:2].[OH:3][c,C:4]>>[C:1](=O)[O:3][C:4]",
    "reductive_amination":  "[CH1,CH2:1]=O.[NH2:2]>>[C:1][N:2]",
    "sulfonamide_form":     "[S:1](=O)(=O)Cl.[NH2:2]>>[S:1](=O)(=O)[N:2]",
    "buchwald_amine":       "[c:1][Br,Cl].[NH2:2]>>[c:1][N:2]",
    "suzuki_coupling":      "[c:1][Br,Cl].[cH:2]>>[c:1][c:2]",
}

def load_compounds(path):
    df = pd.read_csv(path)
    records, n_bad = [], 0
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        if mol is None:
            n_bad += 1
            continue
        fam = str(row["compound_name"]).split("_")[0]
        records.append({
            "compound_name": str(row["compound_name"]),
            "family": fam if fam in FAMILY_COLORS else "other",
            "mol": mol,
        })
    print(f"  {len(records)} valid ({n_bad} skipped)")
    return pd.DataFrame(records)

def test_reaction_applicability(df):
    """Test each compound against each reaction's reactant SMARTS pattern."""
    patterns = {}
    for name, (sma, _) in REACTIONS.items():
        pat = Chem.MolFromSmarts(sma)
        patterns[name] = pat

    rows = []
    for _, row in df.iterrows():
        mol = row["mol"]
        result = {"compound_name": row["compound_name"], "family": row["family"]}
        for name, pat in patterns.items():
            if pat is not None:
                result[name] = 1 if mol.HasSubstructMatch(pat) else 0
            else:
                result[name] = 0
        result["n_applicable"] = sum(result[rxn] for rxn in REACTIONS)
        rows.append(result)
    return pd.DataFrame(rows)

def plot_heatmap(result_df, output_path):
    df = result_df.sort_values(["family", "compound_name"]).reset_index(drop=True)
    heat = df[list(REACTIONS.keys())].astype(int)

    # Rename columns to shorter labels
    col_labels = [k.replace("_", "\n") for k in REACTIONS.keys()]

    fig = plt.figure(figsize=(10, max(8, len(df) * 0.28)))
    ax_h = fig.add_axes([0.22, 0.05, 0.68, 0.88])
    ax_c = fig.add_axes([0.04, 0.05, 0.06, 0.88])

    sns.heatmap(heat, ax=ax_h, annot=True, fmt="d", cmap="Blues",
                linewidths=0.3, linecolor="white",
                xticklabels=col_labels, yticklabels=df["compound_name"],
                vmin=0, vmax=1, cbar=False)
    ax_h.tick_params(axis="y", labelsize=7)
    ax_h.tick_params(axis="x", labelsize=7)
    ax_h.set_title("Reaction SMARTS Applicability per Compound", fontsize=13, fontweight="bold")

    for i, fam in enumerate(df["family"]):
        ax_c.barh(i, 1, color=FAMILY_COLORS.get(fam, "#808080"), edgecolor="none")
    ax_c.set_xlim(0, 1); ax_c.set_ylim(-0.5, len(df) - 0.5)
    ax_c.invert_yaxis(); ax_c.axis("off")

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nLoading: {args.input}")
    df = load_compounds(args.input)

    print("Testing reaction SMARTS applicability...")
    result_df = test_reaction_applicability(df)

    csv_path = os.path.join(args.output_dir, "reaction_applicability.csv")
    result_df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    plot_heatmap(result_df, os.path.join(args.output_dir, "reaction_heatmap.png"))
    print(f"Saved: {args.output_dir}/reaction_heatmap.png")

    print("\n--- Reaction applicability summary ---")
    for rxn in REACTIONS:
        count = result_df[rxn].sum()
        desc = REACTIONS[rxn][1]
        print(f"  {rxn:<22} {count:>3}/{len(result_df)} compounds  ({desc})")

    print("\n--- Applicable reactions per family ---")
    fam_summary = result_df.groupby("family")["n_applicable"].agg(["mean", "min", "max"]).round(2)
    print(fam_summary.to_string())

    print("\n--- Most versatile compounds (most applicable reactions) ---")
    top = result_df.nlargest(5, "n_applicable")[["compound_name", "family", "n_applicable"]]
    print(top.to_string(index=False))
    print("\nDone.")

if __name__ == "__main__":
    main()
