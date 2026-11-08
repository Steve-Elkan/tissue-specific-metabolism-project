import pandas as pd
from pathlib import Path

DATA_DIR = Path("data")

# Load HPA expression dataset
expression = pd.read_csv(
    DATA_DIR / "rna_tissue_consensus.tsv",
    sep="\t"
)

# Load Recon3D model
from cobra.io import read_sbml_model

model = read_sbml_model(
    str(DATA_DIR / "Recon3D.xml")
)

# Normalize gene symbols
expression["Gene name"] = (
    expression["Gene name"]
    .astype(str)
    .str.strip()
    .str.upper()
)

# Tissue names
tissues = ["liver", "skeletal muscle"]

for tissue in tissues:

    print(f"\nProcessing {tissue}")

    # Select tissue expression data
    subset = expression[
        expression["Tissue"].str.lower() == tissue
    ].copy()

    # Average expression for each gene symbol
    symbol_expression = (
        subset.groupby("Gene name")["nTPM"]
        .mean()
    )

    # Map HPA gene symbols to Recon3D gene IDs
    mapped_expression = {}

    for gene in model.genes:

        symbol = gene.name.strip().upper()

        if symbol in symbol_expression.index:

            mapped_expression[gene.id] = (
                symbol_expression[symbol]
            )

    # Create expression dataframe
    mapped = pd.Series(
        mapped_expression,
        name="expression"
    )

    # Save mapped expression
    output = DATA_DIR / f"{tissue}_expression.csv"

    mapped.to_csv(output, header=True)

    # Mapping statistics
    print("Total Recon3D genes:", len(model.genes))
    print("Matched genes:", len(mapped))

    print(
        "Mapping percentage:",
        round(
            len(mapped) / len(model.genes) * 100,
            2
        ),
        "%"
    )

    print("\nExpression statistics:")
    print(mapped.describe())

print("\nGene mapping completed.")

