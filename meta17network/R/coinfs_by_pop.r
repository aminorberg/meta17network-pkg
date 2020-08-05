#' (Co)infections by population
#'
#' Calculate coinfections for host populations
#' @param dirs Directories
#' @param Y Virus community data
#' @param POPs Numerical vector indicating the host populations for Y
#' @param pop_lifes List of selected populations' names, e.g. list("1", "2", ...)
#' @return List of results for the analysis 
#' @export

coinfs_by_pop <- function(Y, pop_lifes, POPs, dirs) {

    for (pop in names(pop_lifes)) {
        pop_life <- Y[which(POPs == pop), ]
        colnames(pop_life) <- sub("sp_", "", colnames(pop_life))
        tmp1 <- toString(colnames(pop_life)[pop_life[1,] == 1])
        for (i in 2:nrow(pop_life)) {
            tmp2 <- toString(colnames(pop_life)[pop_life[i,] == 1])
            tmp1 <- rbind(tmp1, tmp2)
        }
        tmp1[which(tmp1 == "")] <- "No infection"
        pop_lifes[[pop]] <- table(tmp1)
    }
    return(pop_lifes)
}
