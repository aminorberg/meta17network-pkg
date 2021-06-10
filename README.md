# meta17network: an analytical pipeline for illustrating and modelling the structure of wild plant virus communities

A pipeline, related data and an R-package wrapping these up for analysing 
the meta17-project virus community data of wild Plantago lanceolata populations in the
Åland Island, Finland. 

## Installing meta17network

```{r}
### Install devtools (if you do not have them already) and load them
# install.packages("devtools")
library(devtools)

### Install the meta17network package and load it
install_github(repo = "aminorberg/meta17network-pkg", 
			   ref = "HEAD", 
			   subdir = "meta17network")
library(meta17network)
```

## Workflow

All the steps below depend on the succesfull installation of the 'meta17network' package.

1) The 'main' script contains the full workflow from data processing to model fitting and validation. 

2) With the standalone 'coinf_plots' one can make descriptive plots of the local virus communities. 

3) The data used for the illustration and analysis are in folder 'mod_data'.
