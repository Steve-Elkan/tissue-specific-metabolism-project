import pandas as pd
import matplotlib.pyplot as plt
import time

from pathlib import Path
from cobra.core.configuration import Configuration

from imatpy.model_utils import read_model
from imatpy.parse_gpr import gene_to_rxn_weights
from imatpy.imat import imat


# ============================================================
# 1. PATHS AND OUTPUT DIRECTORIES
# ============================================================

DATA = Path("/path/to/directory/data")
RESULTS = Path("results")

(RESULTS / "tables").mkdir(parents=True, exist_ok=True)
(RESULTS / "figures").mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. SOLVER CONFIGURATION
# ============================================================

Configuration().solver = "gurobi"

print("Available solvers:")

from cobra.util.solver import solvers
print(list(solvers.keys()))

if "gurobi" not in solvers:
    raise RuntimeError(
        "Gurobi is not available to COBRApy. "
        "Check your installation."
    )


# ============================================================
# 3. LOAD RECON3D
# ============================================================

print("\n========== LOADING RECON3D ==========")

model = read_model(str(DATA / "Recon3D.xml"))

print("Genes:", len(model.genes))
print("Reactions:", len(model.reactions))
print("Metabolites:", len(model.metabolites))

print("Solver:", model.solver.interface.__name__)

# Configure solver output and timeout.
model.solver.configuration.verbosity = 3
model.solver.configuration.timeout = 1800

# Configure Gurobi directly.
model.solver.problem.Params.OutputFlag = 1
model.solver.problem.Params.Threads = 4
model.solver.problem.Params.Presolve = 2
model.solver.problem.Params.TimeLimit = 1800


# ============================================================
# 4. MODEL FEASIBILITY CHECK
# ============================================================

print("\n========== MODEL FEASIBILITY ==========")

solution = model.optimize()

print("Model status:", solution.status)
print("Objective value:", solution.objective_value)

if solution.status != "optimal":
    raise RuntimeError(
        "Recon3D did not pass the initial feasibility check."
    )


# ============================================================
# 5. EXPERIMENTAL PARAMETERS
# ============================================================

percentiles = [10, 15, 20]

tissues = [
    "liver",
    "skeletal muscle"
]

summary = []
predictions = {}

print("\n========== EXPRESSION ANALYSIS ==========")


# ============================================================
# 6. PROCESS EACH TISSUE
# ============================================================

for tissue in tissues:

    print(f"\nProcessing {tissue}")

    expression = pd.read_csv(
        DATA / f"{tissue}_expression.csv",
        index_col=0
    )["expression"]

    # Ensure gene IDs are strings.
    expression.index = expression.index.astype(str)

    # Ensure expression values are numeric.
    expression = pd.to_numeric(
        expression,
        errors="coerce"
    ).dropna()

    # Inspect gene mapping.
    model_genes = {
        gene.id for gene in model.genes
    }

    matched = set(expression.index) & model_genes

    print("Expression genes:", len(expression))
    print("Matched genes:", len(matched))
    print("Unmatched expression genes:",
          len(expression) - len(matched))

    if tissue == "liver":
        print("\nFirst 20 expression gene IDs:")
        print(expression.index[:20].tolist())

        print("\nFirst 20 Recon3D gene IDs:")
        print([gene.id for gene in model.genes][:20])

    # ========================================================
    # 7. TEST DIFFERENT EXPRESSION THRESHOLDS
    # ========================================================

    for percentile in percentiles:

        print("\n" + "=" * 55)
        print(f"Tissue: {tissue}")
        print(f"Percentile: {percentile}")
        print("=" * 55)

        # Calculate expression cutoffs.
        low_cutoff = expression.quantile(
            percentile / 100
        )

        high_cutoff = expression.quantile(
            1 - percentile / 100
        )

        # Classify genes.
        gene_weights = pd.Series(
            0,
            index=expression.index
        )

        gene_weights[
            expression <= low_cutoff
        ] = -1

        gene_weights[
            expression >= high_cutoff
        ] = 1

        print("Low expression cutoff:", low_cutoff)
        print("High expression cutoff:", high_cutoff)

        print(
            "Low-expression genes:",
            (gene_weights == -1).sum()
        )

        print(
            "High-expression genes:",
            (gene_weights == 1).sum()
        )

        print(
            "Intermediate genes:",
            (gene_weights == 0).sum()
        )

        # Save gene classifications.
        gene_weights.to_csv(
            RESULTS / "tables" /
            f"{tissue}_{percentile}_gene_weights.csv"
        )

        # ====================================================
        # 8. MAP GENE WEIGHTS TO REACTION WEIGHTS
        # ====================================================

        print("\nStarting reaction mapping...", flush=True)

        start = time.time()

        reaction_weights = gene_to_rxn_weights(
            model=model,
            gene_weights=gene_weights
        )

        mapping_time = time.time() - start

        print(
            "Reaction mapping completed:",
            round(mapping_time, 2),
            "seconds"
        )

        print("Total reaction weights:",
              len(reaction_weights))

        print(
            "Positive reaction weights:",
            (reaction_weights > 0).sum()
        )

        print(
            "Negative reaction weights:",
            (reaction_weights < 0).sum()
        )

        print(
            "Neutral reaction weights:",
            (reaction_weights == 0).sum()
        )

        # Save reaction classifications.
        reaction_weights.to_csv(
            RESULTS / "tables" /
            f"{tissue}_{percentile}_reaction_weights.csv"
        )

        # ====================================================
        # 9. RUN iMAT USING AN INDEPENDENT MODEL COPY
        # ====================================================

        print("\nStarting Gurobi iMAT optimization...", flush=True)

        start = time.time()

        # Avoid carrying optimisation modifications
        # from one experiment into another.
        working_model = model.copy()

        # Explicitly configure the copied model.
        working_model.solver = "gurobi"

        working_model.solver.configuration.verbosity = 3
        working_model.solver.configuration.timeout = 1800

        working_model.solver.problem.Params.OutputFlag = 1
        working_model.solver.problem.Params.Threads = 4
        working_model.solver.problem.Params.Presolve = 2
        working_model.solver.problem.Params.TimeLimit = 1800

        try:

            result = imat(
                model=working_model,
                rxn_weights=reaction_weights,
                epsilon=1,
                threshold=0.01
            )

            elapsed = time.time() - start

            print(
                "\niMAT completed in:",
                round(elapsed, 2),
                "seconds"
            )

            print("Optimization status:", result.status)
            print("Objective value:", result.objective_value)

            # Only process optimal solutions.
            if result.status != "optimal":

                print(
                    "WARNING: This condition did not "
                    "reach optimality. Skipping predictions."
                )

                continue

            # =================================================
            # 10. EXTRACT FLUXES
            # =================================================

            fluxes = result.fluxes

            active = fluxes.abs() >= 1
            inactive = fluxes.abs() < 0.01

            key = (tissue, percentile)

            predictions[key] = set(
                fluxes[active].index
            )

            active_count = int(active.sum())
            inactive_count = int(inactive.sum())

            intermediate_count = (
                len(fluxes)
                - active_count
                - inactive_count
            )

            summary.append({
                "tissue": tissue,
                "percentile": percentile,
                "active_reactions": active_count,
                "inactive_reactions": inactive_count,
                "intermediate_reactions": intermediate_count,
                "imat_objective": result.objective_value,
                "runtime_seconds": elapsed,
                "status": result.status
            })

            print("Active reactions:", active_count)
            print("Inactive reactions:", inactive_count)
            print("Intermediate reactions:", intermediate_count)

        except Exception as error:

            print(
                f"ERROR during {tissue}, "
                f"{percentile}% optimization:"
            )

            print(error)

            continue


# ============================================================
# 11. SAVE PARAMETER SWEEP RESULTS
# ============================================================

print("\n========== SAVING RESULTS ==========")

results = pd.DataFrame(summary)

results.to_csv(
    RESULTS / "tables/parameter_sweep.csv",
    index=False
)

print(results)


# ============================================================
# 12. PLOT PARAMETER SENSITIVITY
# ============================================================

if not results.empty:

    plt.figure(figsize=(8, 5))

    for tissue in tissues:

        subset = results[
            results["tissue"] == tissue
        ].sort_values("percentile")

        if not subset.empty:

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

    plt.show()


# ============================================================
# 13. COMPARE LIVER AND SKELETAL MUSCLE
# ============================================================

liver_key = ("liver", 15)
muscle_key = ("skeletal muscle", 15)

if liver_key in predictions and muscle_key in predictions:

    liver = predictions[liver_key]
    muscle = predictions[muscle_key]

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

    print("\nTissue comparison:")
    print(comparison)

else:

    print(
        "\nWARNING: Both 15th percentile predictions "
        "are required for tissue comparison."
    )


# ============================================================
# 14. FINISH
# ============================================================

print("\n========== ANALYSIS COMPLETED ==========")

print("Successful optimizations:", len(summary))
print("Expected optimizations:", len(tissues) * len(percentiles))

print("\nResults saved in:", RESULTS.resolve())
