import pandas as pd
from pathlib import Path

from cobra.io import read_sbml_model

DATA_DIR = Path("/path/to/directory/data")

# Load HPA expression dataset
expression = pd.read_csv(
    DATA_DIR / "rna_tissue_consensus.tsv",
    sep="\t"
)

# Load Recon3D model
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

# Remove missing or invalid gene symbols
expression = expression[
    expression["Gene name"].notna()
    & (expression["Gene name"] != "")
    & (expression["Gene name"] != "NAN")
].copy()

# Ensure expression values are numeric
expression["nTPM"] = pd.to_numeric(
    expression["nTPM"],
    errors="coerce"
)

expression = expression.dropna(
    subset=["nTPM"]
)

tissues = ["liver", "skeletal muscle"]

# Build gene symbol lookup
symbol_expression_data = {}

for tissue in tissues:

    print(f"\nProcessing {tissue}")

    subset = expression[
        expression["Tissue"].str.lower() == tissue
    ].copy()

    symbol_expression = (
        subset.groupby("Gene name")["nTPM"]
        .mean()
    )

    mapped_expression = {}

    for gene in model.genes:

        symbol = gene.name.strip().upper()

        if symbol in symbol_expression.index:

            mapped_expression[gene.id] = (
                symbol_expression[symbol]
            )

    mapped = pd.Series(
        mapped_expression,
        name="expression",
        dtype=float
    )

    output = DATA_DIR / f"{tissue}_expression.csv"

    mapped.to_csv(output, header=True)

    print("Total Recon3D genes:", len(model.genes))
    print("Matched gene IDs:", len(mapped))

    print(
        "Mapping percentage:",
        round(len(mapped) / len(model.genes) * 100, 2),
        "%"
    )

    print("Unique HPA symbols:", len(symbol_expression))

    print("\nExpression statistics:")
    print(mapped.describe())

    print("\nFirst 10 mapped gene IDs:")
    print(mapped.index[:10].tolist())

print("\nGene mapping completed.")
