import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from cobra.core.configuration import Configuration

from imatpy.model_utils import read_model
from imatpy.parse_gpr import gene_to_rxn_weights
from imatpy.imat import imat

DATA = Path("data")
RESULTS = Path("results")

(RESULTS / "tables").mkdir(
    parents=True, exist_ok=True
)
(RESULTS / "figures").mkdir(
    parents=True, exist_ok=True
)

Configuration().solver = "glpk"

# Load the human metabolic network.
model = read_model(str(DATA / "Recon3D.xml"))


# Test different expression thresholds.
percentiles = [10, 15, 20]

tissues = ["liver", "skeletal muscle"]

summary = []
predictions = {}

print("\n========== GENE ID INSPECTION ==========")

print("\nFirst 20 HPA gene IDs:")
print(expression.index[:20].tolist())

print("\nFirst 20 Recon3D gene IDs:")
print([gene.id for gene in model.genes][:20])

print("\nFirst 20 Recon3D gene names:")
print([gene.name for gene in model.genes][:20])



# DIAGNOSTIC CHECKS

print("\n========== MODEL DIAGNOSTICS ==========")

print("Number of model genes:", len(model.genes))
print("Number of reactions:", len(model.reactions))
print("Number of metabolites:", len(model.metabolites))

print("\n========== EXPRESSION DIAGNOSTICS ==========")

for tissue in tissues:

    expression = pd.read_csv(
        DATA / f"{tissue}_expression.csv",
        index_col=0
    )["expression"]

    model_genes = {
        gene.id for gene in model.genes
    }

    matched = set(expression.index) & model_genes

    print(f"\nTissue: {tissue}")
    print("Total expression genes:", len(expression))
    print("Matched genes:", len(matched))

    if len(expression) > 0:
        print("Expression statistics:")
        print(expression.describe())

    if len(matched) > 0:
        print(
            "Matched expression statistics:"
        )
        print(
            expression[
                expression.index.isin(matched)
            ].describe()
        )

print("\n========== MODEL FEASIBILITY ==========")

solution = model.optimize()

print("Model status:", solution.status)
print("Objective value:", solution.objective_value)

print("\n========== END DIAGNOSTICS ==========")


for tissue in tissues:

    print(f"\nProcessing {tissue}")

    expression = pd.read_csv(
        DATA / f"{tissue}_expression.csv",
        index_col=0
    )["expression"]

    # Match expression genes to genes in Recon3D.
    expression = pd.read_csv(
        DATA / f"{tissue}_expression.csv",
        index_col=0
    )["expression"]

    expression.index = expression.index.astype(str)

    for percentile in percentiles:

        # Classify genes using expression percentiles.
        low_cutoff = expression.quantile(
            percentile / 100
        )

        high_cutoff = expression.quantile(
            1 - percentile / 100
        )

        gene_weights = pd.Series(
            0, index=expression.index
        )

        gene_weights[
            expression <= low_cutoff
        ] = -1

        gene_weights[
            expression >= high_cutoff
        ] = 1
        
        print("\nTissue:", tissue)
        print("Threshold:", percentile)

        print("Low expression cutoff:", low_cutoff)
        print("High expression cutoff:", high_cutoff)

        print("Low-expression genes:", (gene_weights == -1).sum())
        print("High-expression genes:", (gene_weights == 1).sum())
        print("Intermediate genes:", (gene_weights == 0).sum())


        import time

        print("Starting reaction mapping...", flush=True)
        start = time.time()

        reaction_weights = gene_to_rxn_weights(
            model=model,
            gene_weights=gene_weights
            )

        print(
            "Reaction mapping completed:",
            round(time.time() - start, 2),
            "seconds",
            flush=True
        )

        # Convert gene weights to reaction weights.
        reaction_weights = gene_to_rxn_weights(
            model=model,
            gene_weights=gene_weights
        )
        
        print("Reaction weights:", len(reaction_weights))
        print("Nonzero reaction weights:", (reaction_weights != 0).sum())
        
        
        print("Starting iMAT optimization...", flush=True)

        start = time.time()

        result = imat(
            model=model,
            rxn_weights=reaction_weights,
            epsilon=1,
            threshold=0.01
        )

        print(
            "iMAT completed:",
            round(time.time() - start, 2),
            "seconds",
            flush=True
            )

        # Run the iMAT optimization.
        result = imat(
            model=model,
            rxn_weights=reaction_weights,
            epsilon=1,
            threshold=0.01
        )

        fluxes = result.fluxes

        active = fluxes.abs() >= 1
        inactive = fluxes.abs() < 0.01

        key = (tissue, percentile)

        predictions[key] = set(
            fluxes[active].index
        )

        summary.append({
            "tissue": tissue,
            "percentile": percentile,
            "active_reactions": int(active.sum()),
            "inactive_reactions": int(inactive.sum()),
            "intermediate_reactions": int(
                len(fluxes) - active.sum() - inactive.sum()
            ),
            "imat_objective": result.objective_value
        })

        print(
            f"{percentile}% threshold: "
            f"objective={result.objective_value}"
        )

# Save results.
results = pd.DataFrame(summary)

results.to_csv(
    RESULTS / "tables/parameter_sweep.csv",
    index=False
)

# Plot active reactions across thresholds.
plt.figure(figsize=(8, 5))

for tissue in tissues:

    subset = results[
        results["tissue"] == tissue
    ]

    plt.plot(
        subset["percentile"],
        subset["active_reactions"],
        marker="o",
        label=tissue
    )

plt.xlabel("Expression percentile threshold (%)")
plt.ylabel("Predicted active reactions")
plt.title("iMAT parameter sensitivity")
plt.legend()
plt.tight_layout()

plt.savefig(
    RESULTS / "figures/parameter_sensitivity.png",
    dpi=300
)

# Compare tissue-specific active reactions.
liver = predictions[("liver", 15)]
muscle = predictions[("skeletal muscle", 15)]

comparison = pd.DataFrame([
    {
        "category": "Shared",
        "reactions": len(liver & muscle)
    },
    {
        "category": "Liver-specific",
        "reactions": len(liver - muscle)
    },
    {
        "category": "Muscle-specific",
        "reactions": len(muscle - liver)
    }
])

comparison.to_csv(
    RESULTS / "tables/tissue_comparison.csv",
    index=False
)

print("\nAnalysis completed.")
print(results)
print(comparison)
