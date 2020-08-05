#' Beta diversity for populations
#'
#' Calculate beta diversity for virus populations using the package 'betapart' and function 'beta.multi'
#' @param Y Virus community data
#' @param POPs Numerical vector indicating the host populations for Y
#' @return List of results for the analysis 
#' @export

do_beta <- function(Y, POPs, dirs) {

    tmp <- c(NA, NA, NA, NA)
    for (i in unique(POPs)) {
        tmp <- rbind(tmp, 
                     c(betapart::beta.multi(Y[which(POPs == i), ], "sorensen"), 
                       i) )
    }
    tmp <- tmp[-1,]
    colnames(tmp)[4] <- c("Population")
    return(tmp)

}
