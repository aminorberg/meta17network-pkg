#' Pairwise beta diversity for populations
#'
#' Calculate beta diversity for virus populations using the package 'betapart' and function 'beta.pair'
#' @param Y Virus community data
#' @param POPs Numerical vector indicating the host populations for Y
#' @return List of results for the analysis 
#' @export

do_beta_pair <- function(Y, POPs, dirs) {

    beta_full <- betapart::beta.pair(Y, "sorensen") 

    beta_pop <- list()
    aggr_tmp <- NA

    for (i1 in 1:length(unique(POPs))) {
        i <-  unique(POPs)[i1]
        beta_pop[[i1]] <- betapart::beta.pair(Y[which(POPs == i), ], "sorensen") 
        aggr_tmp <- rbind(aggr_tmp, colSums(Y[which(POPs == i), ]))
    }
    names(beta_pop) <- unique(POPs)
    
    beta_aggr <- betapart::beta.pair.abund(aggr_tmp[-1,], "bray")

    res <- list("full_data" = beta_full,
                "pop_aggd_abund" = beta_aggr,
                "by_pop" = beta_pop)    
    return(res)

}
