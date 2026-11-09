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
(1) Download the gene-expression dataset and Recon3D model.

(2) Load and preprocess gene-expression measurements.

(3) Match expression genes to Recon3D genes.

(4) Classify genes into low, intermediate, and high expression groups.

(5) Map expression information to metabolic reactions using GPR relationships.

(6) Reconstruct tissue-specific networks using iMAT.

(7) Solve the resulting mixed-integer optimization problems using Gurobi.

(8) Compare predicted reaction activity between tissues.

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
pip install -r requirements.txt

Gurobi requires a separately installed solver and an appropriate license.

7. Running the project
Run the following commands:

scripts/download_data.py

scripts/data_preprocess.py

scripts/run_analysis.py

scripts/analyze_results.py

8. Expected outputs

The analysis generates:
* results/tables/parameter_sweep.csv
* results/tables/tissue_comparison.csv
* results/figures/parameter_sensitivity.png

These files summarize predicted metabolic activity and parameter sensitivity.

9. Results

All six optimization experiments achieved optimal solver status.

10. Interpretation

Increasing the expression percentile resulted in an increase in predicted active reactions in both tissues.
In liver, active reactions increased from 3,806 at percentile 10 to 5,609 at percentile 20.
In skeletal muscle, active reactions increased from 2,300 to 4,367.
This indicates that the reconstructed metabolic networks are sensitive to the expression classification thresholds.
Liver consistently showed more predicted active reactions than skeletal muscle under the tested conditions.
The results demonstrate the importance of parameter selection when integrating gene-expression information into genome-scale metabolic reconstructions.
These reaction classifications represent model-based predictions and should not be interpreted as direct measurements of metabolic flux or biological activity.

11. Tissue Comparison

The initial network comparison identified:
* Shared reactions: 2,227
* Liver-specific reactions: 2,650
* Skeletal-muscle-specific reactions: 790
These results suggest that the two reconstructed tissues contain distinct predicted metabolic capabilities.
The comparison should be interpreted in conjunction with the percentile used for reconstruction and the definitions of active and tissue-specific reactions.

12. Limitations

* Gene expression does not directly measure metabolic flux.
* Missing gene-expression measurements can affect reaction classification.
* Predicted reaction activity depends on expression thresholds and model constraints.
* Recon3D is a generic human metabolic reconstruction, not an experimentally validated tissue-specific network.
* The optimization results depend on the chosen algorithm and solver configuration.

13. References

(1) Shlomi et al. (2008). Network-based prediction of human tissue-specific metabolism. Nature Biotechnology, 26, 1003–1010.

(2) Zur, H., Ruppin, E., & Shlomi, T. (2010). iMAT: an integrative metabolic analysis tool. Bioinformatics, 26, 3140–3142.

(3) Human Protein Atlas: https://www.proteinatlas.org/

(4) Recon3D: https://bigg.ucsd.edu/
