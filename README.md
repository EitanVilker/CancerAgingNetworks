# CancerAgingNetworks

Getting the TCGA SummarizedExperiments for the BRCA assay data requires access to the SCC, so unfortunately the DFE R script cannot be run on its own. However, the CSVs generated from it are included in the results folder.

FBA_With_DFE_Integration.ipynb is the main script, running FBA and combining the flux results with DFE results. The other scripts are mostly for interpreting the results from here, with the exception of Reformatting_DE_Table.ipynb, which gets the related genes to the top 100 genes found from the DFE table for the purpose of assisting with FBA integration.
