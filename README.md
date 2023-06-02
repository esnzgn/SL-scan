# SL-scan
The SL-scan pipeline implementation

<p align="center">
 <img src="graphical abstract_v3.jpg">
</p>

## Data Loading and Filtering

The **1_expression_data_preparations.R** script reads the CCLE expression data from the "../MetabolicSLinput/" directory and applies a filter to include only cancers in RPMI with more than 10 observations, as described in the manuscript.

## FastSL, MCS, gMCS SL pair prediction

The **2.1_SL_pair_prediction_gmcs_mcs_fastsl.m** MATLAB script constructs genome-scale metabolic networks (GSMNs) for each cancer using iMAT, Recon2.0v4, and relative expression data. It performs functionality checks then predicts SL pairs for each cancer model using FastSL, MCS, and gMCS algorithms.

## SL-scan cancer dependency map

The **2.2_SL_pair_prediction_slscan.m** MATLAB script constructs GSMN for each cancer cell line using iMAT, Recon2.0v4, and relative expression data. Then, it constuct a gene dependency vector for each cancer model.

## SL-scan Exhuastive search with CRISPR

The R script **3_MN_exhaustive_ttest_mut_achiles.R** predicts SL pairs for the dependency map of SL-scan and CRISPR data using mutation data through exhaustive T-tests. It then compares the results obtained from other SL pair prediction algorithms with CRISPR data FDR (P-values). Similarly, it performs the same comparison for SL-scan outcomes. Additionally, the script implements a series of hypergeometric tests for each algorithm and each cancer type based on their concordant results with CRISPR.

## SL-scan Exhuastive search with drug perturbation

The R script **4_MN_exhustive_ttest_mut_prism.R** predicts SL pairs for the dependency map of SL-scan and PRISM data using mutation data through exhaustive T-tests. It then compares the results obtained from other SL pair prediction algorithms with PRISM data P-values. Similarly, it performs the same comparison for SL-scan outcomes.

## GSEA and SL pairs cancer generalizationa 

The R script **5_Generalization_and_GeneSetAnalysis_SLscan_and_ttest.R** identifies the driver partnering genes in SL-scan outcomes along with their T-statistic. It then performs a series of Gene Set Enrichment Analysis (GSEA) to identify the most frequent vulnerable biological processes in the SL-scan output. Additionally, the script constructs a bipartite network of SL pairs and cancers to analyze cancer generalization.

## functionality checks violon plot 

The R script **6_functionalities_plot.R** generates violin plots for each cancer metabolic model, comparing their functionality checks to 100 random models.

##
Please note this is a general outline of the script's functionality. For detailed usage instructions and additional information, please refer to the manuscript or contact the authors.




