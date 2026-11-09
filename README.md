# Tissue-Specific-Metabolism-Project
Systems Biology &amp; Network Modeling (CB 5313) Course Project

Network-Based Prediction of Human Tissue-Specific Metabolism

1. Project overview
This project implements the integrative Metabolic Analysis Tool (iMAT), based on the method proposed by Shlomi et al. (2008), to predict tissue-specific metabolic activity in humans.
The analysis focuses on two human tissues:
* Liver
* Skeletal muscle

The project investigates how gene-expression thresholds influence the predicted tissue-specific metabolic networks.

2. Research question
How does the choice of gene-expression thresholds affect the prediction of tissue-specific metabolic activity in human liver and skeletal muscle?

3. Methodology
The project uses iMAT, a constraint-based metabolic modeling method that integrates gene-expression data with a genome-scale metabolic reconstruction.

The workflow consists of:
(1) Downloading human tissue gene-expression data.
(2) Downloading the Recon3D human metabolic model.
(3) Processing and classifying gene-expression values.
(4) Mapping genes to metabolic reactions using gene-protein-reaction associations.
(5) Running iMAT optimization.
(6) Comparing predicted metabolic activity between tissues.
(7) Performing a parameter-sensitivity analysis.

4. Data sources
Gene-expression data
Human Protein Atlas (HPA), consensus RNA tissue dataset.
https://www.proteinatlas.org/humanproteome/tissue/data
Metabolic network
Recon3D human metabolic reconstruction from BiGG Models.
https://bigg.ucsd.edu/models/Recon3D

5. Experimental extension
The original method is extended through a parameter-sensitivity experiment.
Three gene-expression percentile thresholds are investigated:
* 10%
* 15%
* 20%
The predicted metabolic networks are compared across these settings.
The analysis also compares shared and tissue-specific reactions between liver and skeletal muscle.

6. Installation
Install Python dependencies:
bash
pip install -r requirements.txt
The optimization requires a compatible solver. The project uses GLPK through COBRApy.

7. Running the project
Run the following commands from the repository root:
bash
python scripts/download_data.py

python scripts/preprocess.py

python scripts/run_imat.py

8. Expected outputs
The analysis generates:
* results/tables/parameter_sweep.csv
* results/tables/tissue_comparison.csv
* results/figures/parameter_sensitivity.png
These files summarize predicted metabolic activity and parameter sensitivity.

9. Results and interpretation
The results section examines:
* Changes in predicted reaction activity across expression thresholds.
* Differences between liver and skeletal muscle.
* Shared and tissue-specific metabolic reactions.
* The sensitivity and limitations of iMAT predictions.
Interpretations will be based on the actual computational results rather than assumed biological outcomes.

10. References
Shlomi, T., Cabili, M. N., Herrgård, M. J., Palsson, B. Ø., & Ruppin, E. (2008). Network-based prediction of human tissue-specific metabolism. Nature Biotechnology, 26, 1003–1010.
