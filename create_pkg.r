
rm(list = ls(all = TRUE)) 
gc()
#install.packages("devtools")
#devtools::install_github("klutometis/roxygen")

library("devtools")
library("roxygen2")

root_dir <- "/Users/anorberg/Documents/Zurich/UZH/META"
working_dir <- "/Users/anorberg/Documents/Zurich/UZH/META/meta17network-pkg"
pkg_dir <- "/Users/anorberg/Documents/Zurich/UZH/META/meta17network-pkg/meta17network"

### create the package
#setwd(working_dir)
#create("meta17network")

### document the package
setwd(pkg_dir)
document()

### install and load
setwd(working_dir)
install("meta17network")

#library("meta17network")
