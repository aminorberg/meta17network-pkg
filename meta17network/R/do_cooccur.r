#' Co-occurrence analysis
#'
#' Do co-occurrence analysis using the package and function 'cooccur'
#' @param Y Virus community data
#' @param POPs Numerical vector indicating the host populations for Y
#' @return List of results for the analysis 
#' @export

do_cooccur <- function(Y, POPs, dirs) {

    # by population
    tmp_aggr <- NA
    cooccurs_bypop <- list()
    for (i in 1:length(unique(POPs))) {
        pop <- unique(POPs)[i]
        tmp_pop <- Y[which(POPs == pop),]
        tmp_aggr <- rbind(tmp_aggr, colSums(tmp_pop))
        try(cooccurs_bypop[[i]] <- cooccur::cooccur(t(tmp_pop), spp_names = TRUE))
    }
    #saveRDS(cooccurs_bypop, file = file.path(dirs$fits, "cooccurs_by_pop.rds"))

    #cooccurs_bypop[[1]]$results # no significant pairs
    cooc_res <- cooccurs_bypop[[2]]$results
    cooccurs_bypop_res <- cooc_res[cooc_res$p_gt <= 0.05 | cooc_res$p_lt <= 0.05, ]
    for (i in 3:length(cooccurs_bypop)) {
        cooc_res <- cooccurs_bypop[[i]]$results
        try(cooccurs_bypop_res <- rbind(cooccurs_bypop_res, 
                                        cooc_res[cooc_res$p_gt <= 0.05 | cooc_res$p_lt <= 0.05, ]))
    }
    toPlot_pop <- cbind(cooccurs_bypop_res$sp1_name, 
                        cooccurs_bypop_res$sp2_name,
                        cooccurs_bypop_res$obs_cooccur,
                        (cooccurs_bypop_res$obs_cooccur / cooccurs_bypop_res$exp_cooccur))

    # by population, aggregated
    tmp_aggr <- 1 * (tmp_aggr[-1,] > 0)
    cooc_bypop_aggr <- cooccur::cooccur(t(tmp_aggr), spp_names = TRUE)
    cooc_bypop_aggr_sig <- cooc_bypop_aggr$results
    cooc_bypop_aggr_sig <- cooc_bypop_aggr_sig[cooc_bypop_aggr_sig$p_gt <= 0.05 | cooc_bypop_aggr_sig$p_lt <= 0.05, ]
    toPlot_pop_aggr <- cbind(cooc_bypop_aggr_sig$sp1_name, 
                             cooc_bypop_aggr_sig$sp2_name,
                             cooc_bypop_aggr_sig$obs_cooccur,
                             (cooc_bypop_aggr_sig$obs_cooccur / cooc_bypop_aggr_sig$exp_cooccur))


    # by host plant (full data)
    cooc_meta <- cooccur::cooccur(t(Y), spp_names = TRUE)
    #saveRDS(cooc_meta, file = file.path(dirs$fits, "cooccurs_full_data.rds"))


    # significant pairs for the whole data
    cooc_res <- cooc_meta$results
    cooc_meta_sig <- cooc_res[cooc_res$p_gt <= 0.05 | cooc_res$p_lt <= 0.05, ]
    toPlot_all <- cbind(cooc_meta_sig$sp1_name, 
                        cooc_meta_sig$sp2_name,
                        cooc_meta_sig$obs_cooccur,
                        (cooc_meta_sig$obs_cooccur / cooc_meta_sig$exp_cooccur))

    # impute missing pairs to population level results
    toPlot_pop <- cbind(paste(toPlot_pop[,1], toPlot_pop[,2], sep = "_"), toPlot_pop)
    toPlot_all <- cbind(paste(toPlot_all[,1], toPlot_all[,2], sep = "_"), toPlot_all)
    toPlot_pop_aggr <- cbind(paste(toPlot_pop_aggr[,1], toPlot_pop_aggr[,2], sep = "_"), toPlot_pop_aggr)
    colnames(toPlot_pop) <- c("pair", "sp1", "sp2", "obs", "obs_div_exp")
    colnames(toPlot_all) <- colnames(toPlot_pop)
    colnames(toPlot_pop_aggr) <- colnames(toPlot_pop)
    rownames(toPlot_pop) <- toPlot_pop[, "pair"]
    rownames(toPlot_all) <- toPlot_all[, "pair"]
    rownames(toPlot_pop_aggr) <- toPlot_pop_aggr[, "pair"]

    toPlot_pop3 <- rbind(toPlot_all, toPlot_pop_aggr)
    toPlot_pop3[, c("obs", "obs_div_exp")] <- -1
    toPlot_pop3[rownames(toPlot_pop_aggr), ] <- toPlot_pop_aggr
    
    toPlot_pop2 <- toPlot_pop3
    toPlot_pop2[, c("obs", "obs_div_exp")] <- -1
    toPlot_pop2[rownames(toPlot_pop), ] <- toPlot_pop
    
    cooccur_res_list <- list("populations" = toPlot_pop2,
                             "populations_aggregated" = toPlot_pop3, 
                             "full_data" = toPlot_all)
    #saveRDS(cooccur_res_list, file = file.path(dirs$fits, "cooccur_res_list.rds"))

    res <- list("cooccur_population" = cooccurs_bypop, 
                "cooccur_population_aggr" = cooc_bypop_aggr, 
                "cooccur_full_data" = cooc_meta, 
                "significants" = cooccur_res_list)
    return(res)
}
