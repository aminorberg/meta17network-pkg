# An analytical pipeline for illustrating and modelling the structure of wild plant virus
communities

An analytical pipeline, related data and an R-package wrapping these up, for analysing 
the meta17-project virus community data of wild Plantago lanceolata populations in the
Åland Island, Finland.

# Workflow
1) The 'create_pkg' script includes the lines of code needed for creating the package
locally. Note that the package does not install automatically, so if you clone this repository,
you need to install the package as instructed.

2) The 'main' script contains the full workflow from data processing and descriptive
illustrations to model fitting and validation. In the beginning, the user will also define
the required subdirectories.

3) The 'spat_plots' is a standalone script, with which one can make spatial plots of Åland
and the study sites. 

4) All the original, unprocessed data files are in the folder 'data' and the results of
data processing, the data used for the analysis, are in folder 'mod_data'. 
