from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ============================================================
# 1. PROJECT DIRECTORIES
# ============================================================

# Project root directory
ROOT = ROOT = Path.cwd()

# Input and output directories
RESULTS_DIR = ROOT / "results" / "tables"
FIGURES_DIR = ROOT / "results" / "figures"

# Create figures directory if it does not exist
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. PLOT SETTINGS
# ============================================================

sns.set_theme(
    style="whitegrid",
    context="talk",
    palette="Set2"
)

plt.rcParams["figure.dpi"] = 120
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.family"] = "DejaVu Sans"


# ============================================================
# 3. LOAD AND VALIDATE RESULTS
# ============================================================

def load_results():

    file_path = RESULTS_DIR / "parameter_sweep.csv"

    if not file_path.exists():
        raise FileNotFoundError(
            f"\nResults file not found: {file_path}\n"
            "Please run the metabolic optimization script first."
        )

    data = pd.read_csv(file_path)

    required = {
        "tissue",
        "percentile",
        "active_reactions",
        "inactive_reactions",
        "intermediate_reactions",
        "imat_objective",
        "runtime_seconds",
        "status",
    }

    missing = required - set(data.columns)

    if missing:
        raise ValueError(
            f"Missing required columns in CSV: {missing}"
        )

    # Ensure numerical columns are numeric
    numerical_columns = [
        "percentile",
        "active_reactions",
        "inactive_reactions",
        "intermediate_reactions",
        "imat_objective",
        "runtime_seconds",
    ]

    for column in numerical_columns:
        data[column] = pd.to_numeric(
            data[column],
            errors="coerce"
        )

    # Remove rows with missing essential values
    data = data.dropna(
        subset=[
            "tissue",
            "percentile",
            "active_reactions",
            "imat_objective",
        ]
    )

    if data.empty:
        raise ValueError("The results CSV contains no valid data.")

    return data


# ============================================================
# 4. PLOT: ACTIVE REACTIONS
# ============================================================

def plot_active_reactions(data):

    plt.figure(figsize=(10, 7))

    sns.lineplot(
        data=data,
        x="percentile",
        y="active_reactions",
        hue="tissue",
        marker="o",
        linewidth=2.5,
        markersize=10,
    )

    plt.title(
        "Effect of Expression Percentile on Active Reactions",
        fontsize=15,
        fontweight="bold"
    )

    plt.xlabel("Expression percentile")
    plt.ylabel("Number of active reactions")

    plt.tight_layout()

    output = FIGURES_DIR / "active_reactions.png"

    plt.savefig(output, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output}")


# ============================================================
# 5. PLOT: REACTION CATEGORIES
# ============================================================

def plot_reaction_categories(data):

    columns = [
        "active_reactions",
        "inactive_reactions",
        "intermediate_reactions",
    ]

    melted = data.melt(
        id_vars=["tissue", "percentile"],
        value_vars=columns,
        var_name="reaction_category",
        value_name="reaction_count",
    )

    melted["reaction_category"] = (
        melted["reaction_category"]
        .str.replace("_reactions", "", regex=False)
        .str.capitalize()
    )

    melted["experiment"] = (
        melted["tissue"]
        + " - "
        + melted["percentile"].astype(str)
        + "%"
    )

    plt.figure(figsize=(14, 8))

    sns.barplot(
        data=melted,
        x="experiment",
        y="reaction_count",
        hue="reaction_category",
    )

    plt.title(
        "Predicted Reaction Categories",
        fontsize=15,
        fontweight="bold"
    )

    plt.xlabel("Tissue and expression percentile")
    plt.ylabel("Number of reactions")

    plt.xticks(rotation=30, ha="right")

    plt.tight_layout()

    output = FIGURES_DIR / "reaction_categories.png"

    plt.savefig(output, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output}")


# ============================================================
# 6. PLOT: iMAT OBJECTIVE VALUES
# ============================================================

def plot_objective_values(data):

    plt.figure(figsize=(10, 7))

    sns.lineplot(
        data=data,
        x="percentile",
        y="imat_objective",
        hue="tissue",
        marker="o",
        linewidth=2.5,
        markersize=10,
    )

    plt.title(
        "iMAT Objective Values Across Percentiles",
        fontsize=15,
        fontweight="bold"
    )

    plt.xlabel("Expression percentile")
    plt.ylabel("iMAT objective value")

    plt.tight_layout()

    output = FIGURES_DIR / "objective_values.png"

    plt.savefig(output, bbox_inches="tight")
    plt.close()

    print(f"Saved: {output}")


# ============================================================
# 7. PRINT OPTIMIZATION SUMMARY
# ============================================================

def print_summary(data):

    print("\n" + "=" * 65)
    print("METABOLIC MODEL OPTIMIZATION SUMMARY")
    print("=" * 65)

    print("\nComplete optimization results:")

    print(data.to_string(index=False))

    print("\nOptimization status:")

    print(data["status"].value_counts())

    print("\nAverage active reactions by tissue:")

    print(
        data.groupby("tissue")["active_reactions"]
        .mean()
        .round(2)
    )

    print("\nAverage iMAT objective by tissue:")

    print(
        data.groupby("tissue")["imat_objective"]
        .mean()
        .round(2)
    )

    print("\nAverage runtime by percentile:")

    print(
        data.groupby("percentile")["runtime_seconds"]
        .mean()
        .round(2)
    )


# ============================================================
# 8. MAIN FUNCTION
# ============================================================

def main():

    print("\nLoading metabolic optimization results...")

    data = load_results()

    print_summary(data)

    print("\nGenerating figures...")

    plot_active_reactions(data)

    plot_reaction_categories(data)

    plot_objective_values(data)

    print("\n" + "=" * 65)
    print("ANALYSIS COMPLETED SUCCESSFULLY")
    print("=" * 65)

    print("\nFigures saved in:")
    print(FIGURES_DIR)


if __name__ == "__main__":
    main()
