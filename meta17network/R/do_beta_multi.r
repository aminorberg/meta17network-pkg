#' Beta diversity for populations
#'
#' Calculate beta diversity for virus populations using the package 'betapart' and function 'beta.multi'
#' @param Y Virus community data
#' @param POPs Numerical vector indicating the host populations for Y
#' @return List of results for the analysis 
#' @export

do_beta_multi <- function(Y, POPs, dirs) {

    beta_full <- betapart::beta.multi(Y, "sorensen") 

    beta_pop <- c(NA, NA, NA, NA)
    aggr_tmp <- NA
    for (i in unique(POPs)) {
        beta_pop <- rbind(beta_pop, 
                          c(betapart::beta.multi(Y[which(POPs == i), ], "sorensen"), 
                          i) )
        aggr_tmp <- rbind(aggr_tmp, colSums(Y[which(POPs == i), ]))
    }
    beta_pop <- beta_pop[-1,]
    colnames(beta_pop)[4] <- c("Population")
    
    beta_aggr <- betapart::beta.multi.abund(aggr_tmp[-1,], "bray")

    res <- list("full_data" = beta_full,
                "pop_aggd_abund" = beta_aggr,
                "by_pop" = beta_pop)    
    return(res)

}
