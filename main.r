rm(list = ls(all = TRUE)) ; gc()
working_dir <- "/Users/anorberg/Documents/Zurich/UZH/META/meta17network-pkg"
setwd(working_dir)

library("meta17network")
dirs <- set_dirs(working_dir = working_dir)
#saveRDS(dirs, file = file.path(working_dir, "dirs.rds"))

dat <- process_data(dirs = dirs,
                    return_data = TRUE,
                    save_data = TRUE,
                    rmNAs = TRUE)

sp_subset_thr <- 10
sp_subset <- colnames(dat$Y)[colSums(dat$Y) >= sp_subset_thr]
dat1 <- dat
dat1$Y <- dat1$Y[, sp_subset]
y <- as.matrix(dat1$Y)
dim(y)
colSums(y)
y_cocs <- t(y) %*% y

# 2 Markov random fields (MRF)

#devtools::install_github("davharris/rosalia")
#devtools::install_github("hmsc-r/HMSC", build_opts = c("--no-resave-data", "--no-manual"))
#devtools::install_github("davharris/mistnet")
#devtools::install_github("gordy2x/ecoCopula", upgrade = "always")
#devtools::install_github("nicholasjclark/MRFcov")

#library(rosalia)
#library(mistnet)
#library(ecoCopula)
#ft1 <- rosalia(y, prior = make_logistic_prior(scale = 2), trace = FALSE)    
#ft1_1 <- rosalia(y, trace = FALSE)    

library(adespatial)
library(MRFcov)

?MRFcov

# 2.1 Unconditional MRF
mrf1 <- MRFcov(as.data.frame(y), family = "binomial")
saveRDS(mrf1, file = file.path(dirs$fits, "mrf1.rds"))

# 2.2 Conditional MRF: non-spatial, spatial and spatial using Moran's eigenvectors

x <- as.data.frame(dat1$X)

x_nums <- x[, c("plant_size", 
                "kytk17", 
                "prop_agri_area",
                "plm2", 
                "shannon", 
                "severe_winter_days", 
                "temp_eff_days_summer16")]
x_nums_scaled <- scale(x_nums)

x_bools <- x[, c("moth", 
                 "miner",
                 "spittle",
                 "holes",
                 "suck_or_bite",
                 "virus_symptoms")]

x_facs <- as.data.frame(x[, "pop"])
colnames(x_facs) <- "pop"
x_facs$pop <- as.factor(x_facs$pop)

x_spat <- x[, c("x", "y", "x_pop", "y_pop")]

mor_nonnull <- dbmem(x_spat[, c("x", "y")], store.listw = FALSE, MEM.autocor = "non-null")
mor_pos <- dbmem(x_spat[, c("x", "y")], store.listw = FALSE, MEM.autocor = "positive")
#plot(mor1[,1:4], x_spat[, c("x", "y")])

mems1 <- cbind(mor_nonnull$MEM1, mor_nonnull$MEM2, mor_nonnull$MEM3, mor_nonnull$MEM4)

data_df <- as.data.frame(cbind(y, x_nums_scaled, x_bools))
data_mems_df <- as.data.frame(cbind(y, x_nums_scaled, x_bools, mems1))
data_only_mems_df <- as.data.frame(cbind(y, mems1))

# 2.2.1 Non-spatial CMRF

crf1 <- MRFcov(data_df, family = "binomial", n_nodes = 16)
saveRDS(crf1, file = file.path(dirs$fits, "crf1.rds"))

# 2.2.2 Spatial CMRF

# plant level coordinates
crf1_spat1 <- MRFcov_spatial(data_df, coords = x_spat[, c("x", "y")], family = "binomial", n_nodes = 16)
saveRDS(crf1_spat1, file = file.path(dirs$fits, "crf1_spat1.rds"))

# plot level coordinates
crf1_spat2 <- MRFcov_spatial(data_df, 
                             coords = x_spat[, c("x_pop", "y_pop")], 
                             family = "binomial", 
                             n_nodes = 16)
saveRDS(crf1_spat2, file = file.path(dirs$fits, "crf1_spat2.rds"))

# 2.2.3 Spatial CMRF with Morans' eigenvectors

# with environment
crf1_mem <- MRFcov(data_mems_df, family = "binomial", n_nodes = 16)
saveRDS(crf1_mem, file = file.path(dirs$fits, "crf1_mem.rds"))

# only MEMs
crf1_only_mem <- MRFcov(data_only_mems_df, family = "binomial", n_nodes = 16)
saveRDS(crf1_only_mem, file = file.path(dirs$fits, "crf1_only_mem.rds"))

# 2.3 Predictions

preds_mrf1 <- predict_MRF(data_df, mrf1)
preds_crf1 <- predict_MRF(data_df, crf1)
preds_crf1_spat1 <- predict_MRF(data_df, crf1_spat1)
preds_crf1_spat2 <- predict_MRF(data_df, crf1_spat2)
preds_crf1_mem <- predict_MRF(data_mems_df, crf1_mem)
preds_crf1_only_mem <- predict_MRF(data_only_mems_df, crf1_only_mem)

# 2.3 AUCs

pROC::auc(response = matrix(y, ncol = 1), 
          predictor = matrix(preds_mrf1$Probability_predictions, ncol = 1))

pROC::auc(response = matrix(y, ncol = 1),   # seems to be the best
          predictor = matrix(preds_crf1$Probability_predictions, ncol = 1))

pROC::auc(response = matrix(y, ncol = 1), 
          predictor = matrix(preds_crf1_spat1$Probability_predictions, ncol = 1))

pROC::auc(response = matrix(y, ncol = 1), 
          predictor = matrix(preds_crf1_spat2$Probability_predictions, ncol = 1))

pROC::auc(response = matrix(y, ncol = 1), 
          predictor = matrix(preds_crf1_mem$Probability_predictions, ncol = 1))

pROC::auc(response = matrix(y, ncol = 1), 
          predictor = matrix(preds_crf1_only_mem$Probability_predictions, ncol = 1))


# 2.4 MRF Cross Validation And Assessment Of Predictive Performance

cv_comparison_rep100 <- cv_MRF_diag_rep(data_mems_df, 
                                        n_nodes = 16,
                                        compare_null = TRUE,
                                        plot = FALSE,
                                        n_cores = 2, 
                                        family = "binomial",
                                        n_fold_runs = 100)
saveRDS(cv_comparison_rep100, 
        file = file.path(dirs$fits, "cv_comparison_mrf1_crf1_mem_rep100.rds"))

plot(cv_comparison_rep100$mean_specificity[which(cv_comparison_rep100$model != "CRF")], 
     type = "l",
     lwd = 3,
     ylim = c(0.98, 1))
lines(cv_comparison_rep100$mean_specificity[which(cv_comparison_rep100$model == "CRF")],
      lwd = 3,
      col = "red")

# 2.5 Final CRF with bootstrapping

bootedCRF_w_mems <- bootstrap_MRF(data_mems_df,
                                  n_nodes = 16,
                                  family = "binomial",
                                  sample_prop = 0.7,
                                  n_bootstraps = 100)
saveRDS(bootedCRF_w_mems, 
        file = file.path(dirs$fits, "bootedCRF_w_mems_n100.rds"))
        
bootedCRF_w_mems <- readRDS(file = file.path(dirs$fits, "bootedCRF_w_mems_n100.rds"))
names(bootedCRF_w_mems)

# 3 Figures
library(circleplot)

# direct associations
MRF_cov_coefs <- plotMRF_hm(MRF_mod = bootedCRF_w_mems,
                            node_names = colnames(y),
                            main = "Predicted interactions (95% CIs)")

#grid::grid.draw(MRF_cov_coefs)
MRF_cov_coefs

pdf(file = file.path(dirs$figs,
                     "associations.pdf"),
    bg = "transparent", 
    width = 18, 
    height = 6)
    par(family = "serif")
    plot(MRF_cov_coefs)
dev.off()


# key coefficients
colSums(y)
bootedCRF_w_mems$mean_key_coefs$Betaflexiviridae
bootedCRF_w_mems$mean_key_coefs$Fimoviridae
bootedCRF_w_mems$mean_key_coefs$Bromoviridae

viruses <- c("Betaflexiviridae", 
             "Fimoviridae", 
             "Bromoviridae", 
             "Closteroviridae", 
             "Caulimoviridae", 
             "Geminiviridae")
for (whichvirus in viruses) {
    tmp <- data.frame(sp1 = bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]]$Variable, 
                      sp2 = rep(whichvirus, nrow(bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]])),
                      coef = bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]]$Rel_importance)
    pdf(file = file.path(dirs$figs,
                         paste0("key_ints_", whichvirus, ".pdf")),
        bg = "transparent", 
        width = 6, 
        height = 6)
        par(family = "serif")
        circleplot(tmp,
                   cluster = FALSE,
                   style = "classic",
                   plot.control = list(point.labels = TRUE,
                                   cex.point = 15,
                                   line.breaks = c(-2,-1,0,0.25,1,100),
                                   line.cols = c("#004080",
                                                 "#62b1ff",
                                                 "#f5d3dc", 
                                                 "#e795aa", 
                                                 "#c60032"),
                                   line.widths = 5))
    dev.off()
}


