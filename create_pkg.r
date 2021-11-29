
rm(list = ls(all = TRUE)) 
gc()

# required packages for installing the package
#install.packages("devtools")
#devtools::install_github("klutometis/roxygen")

library("devtools")
library("roxygen2")

# you only need to define the directory where you habe cloned the meta17network-pkg repository
root_dir <- "/Users/annanorberg/switchdrive/UZH/Projects/META/"
# after this you define the other directories and install the package
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

# check that you can succesfully open up the library
# after this, the scripts of main.r and spat_plots.r should work!
library("meta17network")
