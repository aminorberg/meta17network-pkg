
rm(list = ls(all = TRUE)) 
gc()
#install.packages("devtools")
#devtools::install_github("klutometis/roxygen")

library("devtools")
library("roxygen2")

root_dir <- "/Users/annanorberg/switchdrive/UZH/Projects/META"
working_dir <- file.path(root_dir, "meta17network-pkg")
pkg_dir <- file.path(working_dir, "meta17network")

### create the package
#setwd(working_dir)
#create("meta17network")

### document the package
setwd(pkg_dir)
document()

### install and load
setwd(working_dir)
install("meta17network")

library("meta17network")

