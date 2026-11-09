Human Tissue-Specific Metabolic Network Reconstruction
A computational study inspired by Shlomi et al. (2008)

Project Description

This project investigates the reconstruction of tissue-specific human metabolic networks by integrating gene-expression information with a genome-scale metabolic model.

The approach is inspired by:

Shlomi, T., Cabili, M. N., Herrgård, M. J., Palsson, B. Ø., & Ruppin, E. (2008). Network-based prediction of human tissue-specific metabolism. Nature Biotechnology, 26, 1003–1010.

https://doi.org/10.1038/nbt.1487

The project uses the Recon3D human metabolic reconstruction, Human Protein Atlas tissue-expression data, and the iMAT optimization method.

Research Question

How does changing gene-expression percentile thresholds affect the reconstruction of tissue-specific metabolic networks in human liver and skeletal muscle?

Datasets

Human Protein Atlas

RNA tissue consensus dataset:

https://www.proteinatlas.org/download/

Provides gene-expression information across human tissues.

Recon3D

Human genome-scale metabolic reconstruction:

https://bigg.ucsd.edu/static/models/Recon3D.xml

Provides metabolic reactions, metabolites, genes, and gene-protein-reaction associations.

Both datasets are downloaded using scripts/download_data.py.

Methodology

The analysis consists of the following steps:

Download the gene-expression dataset and Recon3D model.
Load and preprocess gene-expression measurements.
Match expression genes to Recon3D genes.
Classify genes into low, intermediate, and high expression groups.
Map expression information to metabolic reactions using GPR relationships.
Reconstruct tissue-specific networks using iMAT.
Solve the resulting mixed-integer optimization problems using Gurobi.
Compare predicted reaction activity between tissues.
Experimental Modification: Parameter Sensitivity Analysis

The original project workflow was extended to investigate the effect of expression percentile thresholds.

Three percentile settings were tested:

10th percentile
15th percentile
20th percentile

Two tissues were analyzed:

Liver
Skeletal muscle
This resulted in six optimization experiments.

Results

All six optimization experiments achieved optimal solver status.


Tissue	Percentile	Active	Inactive	Intermediate
Liver	10	3806	6780	14
Liver	15	4877	5661	62
Liver	20	5609	4902	89
Skeletal muscle	10	2300	8265	35
Skeletal muscle	15	3017	7528	55
Skeletal muscle	20	4367	6094	139
Interpretation

Increasing the expression percentile resulted in an increase in predicted active reactions in both tissues.

In liver, active reactions increased from 3,806 at percentile 10 to 5,609 at percentile 20.

In skeletal muscle, active reactions increased from 2,300 to 4,367.

This indicates that the reconstructed metabolic networks are sensitive to the expression classification thresholds.

Liver consistently showed more predicted active reactions than skeletal muscle under the tested conditions.

The results demonstrate the importance of parameter selection when integrating gene-expression information into genome-scale metabolic reconstructions.

These reaction classifications represent model-based predictions and should not be interpreted as direct measurements of metabolic flux or biological activity.

Tissue Comparison

The initial network comparison identified:

Shared reactions: 2,227
Liver-specific reactions: 2,650
Skeletal-muscle-specific reactions: 790
These results suggest that the two reconstructed tissues contain distinct predicted metabolic capabilities.

The comparison should be interpreted in conjunction with the percentile used for reconstruction and the definitions of active and tissue-specific reactions.

Software Requirements

Python 3.10 or later
COBRApy
iMATpy
Gurobi Optimizer
Gurobi Python interface
Pandas
NumPy
Matplotlib
Seaborn
Gurobi requires a separately installed solver and an appropriate license.

Installation

Clone the repository:


bash
git clone YOUR_GITHUB_REPOSITORY_URL
cd human-tissue-specific-metabolism
Create a virtual environment:


bash
python -m venv venv
source venv/bin/activate
On Windows:


bash
venv\Scripts\activate
Install dependencies:


bash
pip install -r requirements.txt
Install and configure Gurobi and iMATpy according to their respective documentation.

Download Data

Run:


bash
python scripts/download_data.py
This downloads the expression dataset and Recon3D model into the data/ directory.

Run the Analysis

Run the metabolic reconstruction pipeline:


bash
python scripts/run_analysis.py
Generate figures from the optimization results:


bash
python scripts/analyze_results.py
Repository Structure

Limitations

Gene expression does not directly measure metabolic flux.
Missing gene-expression measurements can affect reaction classification.
Predicted reaction activity depends on expression thresholds and model constraints.
Recon3D is a generic human metabolic reconstruction, not an experimentally validated tissue-specific network.
The optimization results depend on the chosen algorithm and solver configuration.

References

Shlomi et al. (2008). Network-based prediction of human tissue-specific metabolism. Nature Biotechnology, 26, 1003–1010.
Zur, H., Ruppin, E., & Shlomi, T. (2010). iMAT: an integrative metabolic analysis tool. Bioinformatics, 26, 3140–3142.
Human Protein Atlas: https://www.proteinatlas.org/
Recon3D: https://bigg.ucsd.edu/
