#' Defining directories
#'
#' Defining directories/paths to be used
#' @param working_dir Full path to the working directory, under which everything else will be created.
#' @param data_fold Folder for the modified data (defaults to "mod_data")
#' @param fit_fold Folder where the model fits will be saved  (defaults to "fits")
#' @return List of directories needed for the modelling pipeline 
#' @export

set_dirs <- function(working_dir,
                     fit_fold = "fits")
{
    res <- structure(list(), class = "dirlist")

    res$wd <- working_dir
    res$fits <- file.path(working_dir, fit_fold)
    res$figs <- file.path(working_dir, "figs")
    res$dat <- file.path(working_dir, "data")
    res$aland <- file.path(working_dir, "data", "aland")
    res$mod_dat <- file.path(working_dir, "mod_data")
    
    return(res)
}
