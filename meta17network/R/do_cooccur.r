#' Co-occurrence analysis
#'
#' Do co-occurrence analysis using the package and function 'cooccur'
#' @param Y Virus community data
#' @param POPs Numerical vector indicating the host populations for Y
#' @return List of results for the analysis 
#' @export

do_cooccur <- function(Y, POPs, dirs) {

    # by population
    cooccurs_bypop <- list()
    for (i in 1:length(unique(POPs))) {
        pop <- unique(POPs)[i]
        tmp_pop <- Y[which(POPs == pop),]
        try(cooccurs_bypop[[i]] <- cooccur::cooccur(t(tmp_pop), spp_names = TRUE))
    }
    cooccurs_bypop
    #saveRDS(cooccurs_bypop, file = file.path(dirs$fits, "cooccurs_by_pop.rds"))

    #cooccurs_bypop[[1]]$results # no significant pairs
    cooc_res <- cooccurs_bypop[[2]]$results
    cooccurs_bypop_res <- cooc_res[cooc_res$p_gt <= 0.05 | cooc_res$p_lt <= 0.05, ]
    for (i in 3:length(cooccurs_bypop)) {
        cooc_res <- cooccurs_bypop[[i]]$results
        try(cooccurs_bypop_res <- rbind(cooccurs_bypop_res, 
                                        cooc_res[cooc_res$p_gt <= 0.05 | cooc_res$p_lt <= 0.05, ]))
    }

    # by host plant (full data)
    cooc_meta <- cooccur::cooccur(t(Y), spp_names = TRUE)
    #saveRDS(cooc_meta, file = file.path(dirs$fits, "cooccurs_full_data.rds"))

    # significant pairs for host populations
    toPlot_pop <- cbind(cooccurs_bypop_res$sp1_name, 
                        cooccurs_bypop_res$sp2_name,
                        cooccurs_bypop_res$obs_cooccur,
                        (cooccurs_bypop_res$obs_cooccur / cooccurs_bypop_res$exp_cooccur))


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
    colnames(toPlot_pop) <- c("pair", "sp1", "sp2", "obs", "obs_div_exp")
    colnames(toPlot_all) <- colnames(toPlot_pop)
    rownames(toPlot_pop) <- toPlot_pop[, "pair"]
    rownames(toPlot_all) <- toPlot_all[, "pair"]
    toPlot_pop2 <- toPlot_all
    toPlot_pop2[, c("obs", "obs_div_exp")] <- -1
    toPlot_pop2[rownames(toPlot_pop), ] <- toPlot_pop
    cooccur_res_list <- list("populations" = toPlot_pop2, "full_data" = toPlot_all)
    #saveRDS(cooccur_res_list, file = file.path(dirs$fits, "cooccur_res_list.rds"))

    return(list("cooccur_population" = cooccurs_bypop, 
                "cooccur_full_data" = cooc_meta, 
                "significants" = cooccur_res_list))
}
