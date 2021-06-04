# An analytical pipeline for illustrating and modelling the structure of wild plant virus communities
An analytical pipeline, related data and an R-package wrapping these up, for analysing 
the meta17-project virus community data of wild Plantago lanceolata populations in the Åland Island, Finland.
The main script contains the workflow for descriptive analysis as well as fitting 
Conditional Markov Random Field models.

# Workflow
The 'main' script contains the full workflow from data processing and descriptive illustrations to model fitting and validation. 
The 'spat_plots' is a standalone script, with which one can make spatial plots of Åland and the study sites. All the original, 
unprocessed data files are in the folder 'data' and the results of data processing, the data used for the analysis, are in folder 'mod_data'. 