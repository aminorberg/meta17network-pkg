
rm(list = ls(all = TRUE)) ; gc()
working_dir <- "/Users/anorberg/Documents/Zurich/UZH/META/meta17network-pkg"
setwd(working_dir)
library("meta17network")
dirs <- set_dirs(working_dir = working_dir)
#saveRDS(dirs, file = file.path(working_dir, "dirs.rds"))

#devtools::install_github("nicholasjclark/MRFcov")

# 1 Data processing and preliminary analysis
dat <- process_data(dirs = dirs,
                    return_data = TRUE,
                    save_data = TRUE,
                    rmNAs = TRUE)

# virus prevalences
png(file.path(dirs$figs, "prevalences.png"),
    height = 4, 
    width = 6, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(8,3,1,1))
    barplot(colSums(dat$Y)[order(colSums(dat$Y))], las = 2, ylim = c(0, 150))
dev.off()

# 1.1 Co-occurrence analysis
res_cooc <- do_cooccur(Y = dat$Y, POPs = dat$X$pop, dirs = dirs)
saveRDS(res_cooc, file = file.path(dirs$fits, "res_cooc.rds"))

# 1.2 Host population level beta diversity
betas <- do_beta(Y = dat$Y, POPs = dat$X$pop, dirs = dirs)
betas <- apply(betas, 2, as.numeric)
write.table(round(betas, 3), 
            file = file.path(dirs$fits, "betas.csv"), 
            quote = FALSE, 
            sep = ",", 
            dec = ".")
saveRDS(betas, file = file.path(dirs$fits, "betas.rds"))
#betas <- readRDS(file = file.path(dirs$fits, "betas.rds"))
# populations 861 vs. 3225 and 946 are the most contrasting ones in terms of Simpson

# coinfections per population
coinfs <- coinfs_by_pop(Y = dat$Y, 
                        POPs = dat$X$pop, 
                        dirs = dirs, 
                        pop_lifes = list("861" = NA, "3225" = NA))
saveRDS(coinfs, file = file.path(dirs$fits, "coinfs.rds"))

# coinfections in the full data
y_coinfs <- dat$Y
colnames(y_coinfs) <- sub("sp_", "", colnames(y_coinfs))
colnames(y_coinfs) <- sub("viridae", "", colnames(y_coinfs))
colnames(y_coinfs) <- sub("viroidae", "", colnames(y_coinfs))
colnames(y_coinfs) <- sub("tidae", "", colnames(y_coinfs))
tmp1 <- toString(colnames(y_coinfs)[y_coinfs[1,] == 1])
for (i in 2:nrow(y_coinfs)) {
    tmp2 <- toString(colnames(y_coinfs)[y_coinfs[i,] == 1])
    tmp1 <- rbind(tmp1, tmp2)
}
tmp1[which(tmp1 == "")] <- "No infection"
all_coinfs <- lapply(tmp1, sort)
all_coinfs <- table(sort(unlist(all_coinfs)))
all_coinfs <- all_coinfs[order(all_coinfs)]
saveRDS(all_coinfs, file = file.path(dirs$fits, "all_coinfs.rds"))

# 1.3 Data processing for (C)MRF mocdelling

# 1.3.1 Individual plant level

## Y
sp_subset_thr <- 10
sp_subset <- colnames(dat$Y)[colSums(dat$Y) >= sp_subset_thr]
dat1 <- dat
dat1$Y <- dat1$Y[, sp_subset]
y <- as.matrix(dat1$Y)
dim(y)
colSums(y)
y_cocs <- t(y) %*% y

## X
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
                 "suck_or_bite")]
x_facs <- as.data.frame(x[, "pop"])
colnames(x_facs) <- "pop"
x_facs$pop <- as.factor(x_facs$pop)
x_spat <- x[, c("x", "y", "x_pop", "y_pop")]

#1.3.2 Population level
# Y
y_pop <- cbind(x_facs$pop, y)

y_pop <- NA
for (pop in unique(x_facs$pop)) {
   tmp1 <- colSums(y[which(x_facs$pop == pop),])
   y_pop <- rbind(y_pop, tmp1)
}
y_pop <- y_pop[-1,]
y_pop_pa <- (y_pop > 0) * 1
colSums(y_pop_pa)/nrow(y_pop_pa)
# co-occurrence analysis
cooc_aggr <- cooccur::cooccur(t(y_pop_pa), spp_names = TRUE)


# X
x_nums_pop <- NA
for (pop in unique(x_facs$pop)) {
   tmp1 <- colMeans(x_nums[which(x_facs$pop == pop),])
   x_nums_pop <- rbind(x_nums_pop, tmp1)
}
x_nums_pop <- x_nums_pop[-1,]
x_nums_pop_scaled <- scale(x_nums_pop)

x_bools_pop <- NA
for (pop in unique(x_facs$pop)) {
   tmp1 <- colSums(x_bools[which(x_facs$pop == pop),])
   x_bools_pop <- rbind(x_bools_pop, tmp1)
}
x_bools_pop <- x_bools_pop[-1,]


# 1.3.3 Spatial eigenvector maps

## sample level
mor_nonnull <- adespatial::dbmem(x_spat[, c("x", "y")], 
                                 store.listw = FALSE,
                                 MEM.autocor = "non-null")
plot(mor_nonnull)

mems1 <- cbind(mor_nonnull$MEM1, 
               mor_nonnull$MEM2, 
               mor_nonnull$MEM3, 
               mor_nonnull$MEM4)
     
## population level
#mor_pop_nonnull <- adespatial::dbmem(unique(x_spat[, c("x_pop", "y_pop")]), 
#                                     store.listw = FALSE, 
#                                     MEM.autocor = "non-null")
#mems1_pop <- cbind(mor_pop_nonnull$MEM1, 
#                   mor_pop_nonnull$MEM2, 
#                   mor_pop_nonnull$MEM3, 
#                   mor_pop_nonnull$MEM4)

## combine Xs
x_and_mems <- as.data.frame(cbind(x_nums_scaled, x_bools, mems1))
#x_and_mems_pop <- as.data.frame(cbind(x_nums_pop_scaled, x_bools_pop, mems1_pop))

# 1.3.4 Check for correlations

# sample level
var_cors <- cor(x_and_mems, use = "na.or.complete")
max(var_cors[upper.tri(var_cors, diag = FALSE)]) # highest correlation 0.56

# population level
#var_pop_cors <- cor(x_and_mems_pop, use = "na.or.complete")
#max(var_pop_cors[upper.tri(var_pop_cors, diag = FALSE)]) # highest correlation 0.63

# 1.3.5 Final data frames

## sample level
data_df <- as.data.frame(cbind(y, cbind(x_nums_scaled, x_bools)))
data_mems_df <- as.data.frame(cbind(y, x_and_mems))
data_only_mems_df <- as.data.frame(cbind(y, mems1))

## population level
#data_pop_df <- as.data.frame(cbind(y_pop, cbind(x_nums_pop_scaled, x_bools_pop)))
#data_pop_mems_df <- as.data.frame(cbind(y_pop, x_and_mems_pop))
#data_pop_only_mems_df <- as.data.frame(cbind(y_pop, mems1_pop))

# 2 (Conditional) Markov Random Fields, (C)MRF
library(MRFcov)
?MRFcov

# 2.1 Unconditional MRF

## sample level
mrf1 <- MRFcov(as.data.frame(y), family = "binomial")
saveRDS(mrf1, file = file.path(dirs$fits, "mrf1.rds"))

## population level
#mrf1_pop <- MRFcov(as.data.frame(y_pop), family = "poisson")
#saveRDS(mrf1_pop, file = file.path(dirs$fits, "mrf1_pop.rds"))

# 2.2 Conditional MRF: non-spatial, spatial and spatial using Moran's eigenvectors

# 2.2.1 Non-spatial CMRF

## sample level
crf1 <- MRFcov(data_df, family = "binomial", n_nodes = 16)
saveRDS(crf1, file = file.path(dirs$fits, "crf1.rds"))

#3 population level
#crf1_pop <- MRFcov(data_pop_df, family = "poisson", n_nodes = 16)
#saveRDS(crf1_pop, file = file.path(dirs$fits, "crf1_pop.rds"))

# 2.2.2 Spatial CMRF

## sample level
### with plant level coordinates
crf1_spat1 <- MRFcov_spatial(data_df, 
                             coords = x_spat[, c("x", "y")], 
                             family = "binomial", 
                             n_nodes = 16)
saveRDS(crf1_spat1, file = file.path(dirs$fits, "crf1_spat1.rds"))

### population level coordinates
#crf1_spat2 <- MRFcov_spatial(data_df, 
#                             coords = x_spat[, c("x_pop", "y_pop")], 
#                             family = "binomial", 
#                             n_nodes = 16)
#saveRDS(crf1_spat2, file = file.path(dirs$fits, "crf1_spat2.rds"))
#
### population level
#crf1_pop_spat1 <- MRFcov_spatial(data_pop_df, 
#                                 coords = unique(x_spat[, c("x_pop", "y_pop")]), 
#                                 family = "poisson", 
#                                 n_nodes = 16)
#saveRDS(crf1_pop_spat1, file = file.path(dirs$fits, "crf1_pop_spat1.rds"))

# 2.2.3 Spatial CRF with Morans' eigenvectors

## sample level
### with environment
crf1_mem <- MRFcov(data_mems_df, family = "binomial", n_nodes = 16)
saveRDS(crf1_mem, file = file.path(dirs$fits, "crf1_mem.rds"))

### only MEMs
crf1_only_mem <- MRFcov(data_only_mems_df, family = "binomial", n_nodes = 16)
saveRDS(crf1_only_mem, file = file.path(dirs$fits, "crf1_only_mem.rds"))

#population level
## with environment
#crf1_pop_mem <- MRFcov(data_pop_mems_df, family = "poisson", n_nodes = 16)
#saveRDS(crf1_pop_mem, file = file.path(dirs$fits, "crf1_pop_mem.rds"))
#
## only MEMs
#crf1_pop_only_mem <- MRFcov(data_pop_only_mems_df, family = "poisson", n_nodes = 16)
#saveRDS(crf1_pop_only_mem, file = file.path(dirs$fits, "crf1_pop_only_mem.rds"))

# 2.3 Predictions
# (load the previously fitted models)

## (sample level)
mrf1 <- readRDS(file = file.path(dirs$fits, "mrf1.rds"))
crf1 <- readRDS(file = file.path(dirs$fits, "crf1.rds"))
crf1_spat1 <- readRDS(file = file.path(dirs$fits, "crf1_spat1.rds"))
crf1_spat2 <- readRDS(file = file.path(dirs$fits, "crf1_spat2.rds"))
crf1_mem <- readRDS(file = file.path(dirs$fits, "crf1_mem.rds"))
crf1_only_mem <- readRDS(file = file.path(dirs$fits, "crf1_only_mem.rds"))

## (population level)
#mrf1_pop <- readRDS(file = file.path(dirs$fits, "mrf1_pop.rds"))
#crf1_pop <- readRDS(file = file.path(dirs$fits, "crf1_pop.rds"))
#crf1_pop_spat1 <- readRDS(file = file.path(dirs$fits, "crf1_pop_spat1.rds"))
#crf1_pop_mem <- readRDS(file = file.path(dirs$fits, "crf1_pop_mem.rds"))
#crf1_pop_only_mem <- readRDS(file = file.path(dirs$fits, "crf1_pop_only_mem.rds"))

## sample level
preds_mrf1 <- predict_MRF(data_df, mrf1)
preds_crf1 <- predict_MRF(data_df, crf1)
preds_crf1_spat1 <- predict_MRF(data_df, crf1_spat1)
preds_crf1_spat2 <- predict_MRF(data_df, crf1_spat2)
preds_crf1_mem <- predict_MRF(data_mems_df, crf1_mem)
preds_crf1_only_mem <- predict_MRF(data_only_mems_df, crf1_only_mem)

## population level
#preds_mrf1_pop <- predict_MRF(data_pop_df, mrf1_pop)
#preds_crf1_pop <- predict_MRF(data_pop_df, crf1_pop)
#preds_crf1_pop_spat1 <- predict_MRF(data_pop_df, crf1_pop_spat1)
#preds_crf1_pop_mem <- predict_MRF(data_pop_mems_df, crf1_pop_mem)
#preds_crf1_pop_only_mem <- predict_MRF(data_pop_only_mems_df, crf1_pop_only_mem)


# 2.4 Model fit evaluation

## AUCs for sample level
pROC::auc(response = matrix(y, ncol = 1), 
          predictor = matrix(preds_mrf1$Probability_predictions, ncol = 1))
pROC::auc(response = matrix(y, ncol = 1),
          predictor = matrix(preds_crf1$Probability_predictions, ncol = 1))
pROC::auc(response = matrix(y, ncol = 1), 
          predictor = matrix(preds_crf1_spat1$Probability_predictions, ncol = 1))
pROC::auc(response = matrix(y, ncol = 1), 
          predictor = matrix(preds_crf1_spat2$Probability_predictions, ncol = 1))
pROC::auc(response = matrix(y, ncol = 1),       # seems to be the best (marginally)
          predictor = matrix(preds_crf1_mem$Probability_predictions, ncol = 1))
pROC::auc(response = matrix(y, ncol = 1),  
          predictor = matrix(preds_crf1_only_mem$Probability_predictions, ncol = 1))

tmp2 <- NA
for (i in 1:ncol(y)) {
    tmp1 <- pROC::auc(response = y[,i],
                      predictor = preds_crf1_mem$Probability_predictions[,i])    
    tmp2 <- c(tmp2, tmp1)
}
aucs <- matrix(tmp2[-1], ncol = 1)
rownames(aucs) <- colnames(y)
mean(aucs[,1])
# all above 0.75, except Closteroviridae 0.66

# Tjur R2 : for each of the two categories of the dependent variable, 
# calculate the mean of the predicted probabilities of an event,
# then, take the difference between those two means.

# crf
probs <- preds_crf1_mem$Probability_predictions
# mrf
probs <- preds_mrf1$Probability_predictions

mean(probs[y == 1]) - mean(probs[y == 0])   # for the whole data
tmp <- NA
for (i in 1:ncol(y)){                       # for all species separately
    print(c(colnames(y)[i], 
            round(mean(probs[,i][y[,i] == 1]) - mean(probs[,i][y[,i] == 0]), 2)))
            tmp <- c(tmp, mean(probs[,i][y[,i] == 1]) - mean(probs[,i][y[,i] == 0]))
            
}
tjuR2s <- matrix(tmp[-1], ncol = 1)
rownames(tjuR2s) <- colnames(y)
mean(tjuR2s[,1])

## correlations for population level
#cor(matrix(y_pop, ncol = 1), 
#    matrix(preds_mrf1_pop, ncol = 1))
#cor(matrix(y_pop, ncol = 1), 
#    matrix(preds_crf1_pop, ncol = 1))
#cor(matrix(y_pop, ncol = 1), 
#    matrix(preds_crf1_pop_spat1, ncol = 1))
#cor(matrix(y_pop, ncol = 1), 
#    matrix(preds_crf1_pop_mem, ncol = 1))       # seems to be the best
#cor(matrix(y_pop, ncol = 1), 
#    matrix(preds_crf1_pop_only_mem, ncol = 1))

### the best model at both levels is CRF with MEMs
##########################################################################################

# 2.5 Final CRF with bootstrapping

# 2.5.1 Sample level
bootedCRF_w_mems <- bootstrap_MRF(data_mems_df,
                                  n_nodes = 16,
                                  family = "binomial",
                                  sample_prop = 0.7,
                                  n_bootstraps = 1000)
saveRDS(bootedCRF_w_mems, 
        file = file.path(dirs$fits, "bootedCRF_w_mems_n1000.rds"))

bootedCRF_w_mems <- bootstrap_MRF(data_mems_df,
                                  n_nodes = 16,
                                  family = "binomial",
                                  sample_prop = 1,
                                  n_bootstraps = 100)
saveRDS(bootedCRF_w_mems, 
        file = file.path(dirs$fits, "bootedCRF_w_mems_sampleprop100_n100.rds"))

# corresponding MRF
bootedMRF <- bootstrap_MRF(as.data.frame(y),
                           n_nodes = 16,
                           family = "binomial",
                           sample_prop = 1,
                           n_bootstraps = 100)
saveRDS(bootedMRF, 
        file = file.path(dirs$fits, "bootedMRF_sampleprop100_n100.rds"))

# 2.5.1 Population level
#bootedCRF_w_mems_pop <- bootstrap_MRF(data_pop_mems_df,
#                                  n_nodes = 16,
#                                  family = "poisson",
#                                  sample_prop = 0.7,
#                                  n_bootstraps = 1000)
#saveRDS(bootedCRF_w_mems_pop, 
#        file = file.path(dirs$fits, "bootedCRF_w_mems_pop_n1000.rds"))


# 2.6 MRF Cross Validation And Assessment Of Predictive Performance

# 2.6.1 Sample level
cv_comparison_rep1000 <- cv_MRF_diag_rep(data_mems_df, 
                                         n_nodes = 16,
                                         compare_null = TRUE,
                                         plot = FALSE,
                                         n_cores = 2, 
                                         family = "binomial",
                                         n_fold_runs = 1000)
                                        
saveRDS(cv_comparison_rep1000, 
        file = file.path(dirs$fits, "cv_comparison_mrf1_crf1_mem_rep1000.rds"))
cv_comparison_rep1000 <- readRDS(file = file.path(dirs$fits, "cv_comparison_mrf1_crf1_mem_rep1000.rds"))

mod_fits <- cv_comparison_rep1000
names(mod_fits)

# sensitivity (aka the true positive rate) is
# the proportion of presences that are correctly classified as such
# unconditional model
quantile(mod_fits$mean_sensitivity[which(mod_fits$model != "CRF")],
         probs = c(0.01, 0.5, 0.99))
# conditional model
quantile(mod_fits$mean_sensitivity[which(mod_fits$model == "CRF")],
         probs = c(0.01, 0.5, 0.99))

# specificity (aka the true negative rate) is
# the proportion of absences that are correctly classified as such
# unconditional model
quantile(mod_fits$mean_specificity[which(mod_fits$model != "CRF")],
         probs = c(0.01, 0.5, 0.99))
# conditional model
quantile(mod_fits$mean_specificity[which(mod_fits$model == "CRF")],
         probs = c(0.01, 0.5, 0.99))

# proportions of true infections that were correctly predicted
# unconditional model
quantile(mod_fits$mean_tot_pred[which(mod_fits$model != "CRF")],
         probs = c(0.01, 0.5, 0.99))
# conditional model
quantile(mod_fits$mean_tot_pred[which(mod_fits$model == "CRF")],
         probs = c(0.01, 0.5, 0.99))

# model's ability to correctly classify positive infections
# unconditional model
quantile(mod_fits$mean_pos_pred[which(mod_fits$model != "CRF")],
         probs = c(0.01, 0.5, 0.99),
         na.rm = TRUE)
# conditional model
quantile(mod_fits$mean_pos_pred[which(mod_fits$model == "CRF")],
         probs = c(0.01, 0.5, 0.99),
         na.rm = TRUE)

# 2.6.2 Population level
#cv_comparison_pop_rep1000 <- cv_MRF_diag_rep(data_mems_df, 
#                                             n_nodes = 16,
#                                             compare_null = TRUE,
#                                             plot = FALSE,
#                                             n_cores = 2, 
#                                             family = "poisson",
#                                             n_fold_runs = 1000)
#                                        
#saveRDS(cv_comparison_pop_rep1000, 
#        file = file.path(dirs$fits, "cv_comparison_mrf1_crf1_mem_pop_rep1000.rds"))

# 3 Results & Figures
library(circleplot)

# 3.1 Beta diversity

###* ALL COINFECTIONS *###

all_coinfs <- readRDS(file = file.path(dirs$fits, "all_coinfs.rds"))

png(file.path(dirs$figs,
              "coinfections.png"), 
    height = 25, 
    width = 20, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(mar = c(70, 3, 1, 1), family = "serif")
    barplot(all_coinfs, las = 2)
dev.off()

###* end of ALL COINFECTIONS *###


###* FIGURE 1 IN MAIN TEXT *###

# maps
make_maps(dirs = dirs, POPs = dat$X$pop)
# (co)infection profile piecharts
plot_pops(dat = dat, dirs = dirs)

# beta diversities
coinfs <- readRDS(file = file.path(dirs$fits, "coinfs.rds"))
all_combs_frqs <- table(unlist(lapply(coinfs, names)))
all_combs_frqs <- all_combs_frqs[order(all_combs_frqs, decreasing = TRUE)]
all_combs_frqs_df <- as.data.frame(all_combs_frqs)
viruses_split <- strsplit(as.character(all_combs_frqs_df[,1]), ",")
viruses_nos <- unlist(lapply(viruses_split, length))

library(wesanderson)
wesandcols_sel <- c(wes_palette("GrandBudapest1"),
                    wes_palette("GrandBudapest2"),
                    wes_palette("Zissou1"),
                    wes_palette("Rushmore"))

combs_cols <- wesandcols_sel[1:nrow(all_combs_frqs_df)]
all_coinf_combs_cols <- cbind(as.character(all_combs_frqs_df[,1]), combs_cols)
rownames(all_coinf_combs_cols) <- all_coinf_combs_cols[,1]
all_coinf_combs_cols <- all_coinf_combs_cols[order(viruses_nos), ]
all_coinf_combs_cols[which(rownames(all_coinf_combs_cols) == "No infection"), 2] <- "white"

#pop <- "3225"
#par(family = "serif", mar = c(0,0,0,0), mfrow = c(1,2))
for (pop in names(coinfs)) {
    png(file.path(dirs$figs,
                 "coinfections_by_pop", 
                 paste0("coinfection_profile_pop_", pop,".png")), 
        height = 3, 
        width = 3, 
        bg = "transparent",
        units = "in", 
        res = 300)

    patch_cols <- all_coinf_combs_cols[sort(match(names(coinfs[[pop]]), rownames(all_coinf_combs_cols))), 2]
    par(mar = rep(0, 4))
    pie(coinfs[[pop]][names(patch_cols)], labels = NA, col = patch_cols)

    dev.off()
}
all_coinf_combs_cols_leg <- all_coinf_combs_cols[-which(all_coinf_combs_cols[,1] == "No infection"),]
png(file = file.path(dirs$figs, 
                     "coinfections_by_pop", 
                     "coinfection_profile_legend.png"),
    bg = "white", 
    width = 12, 
    height = 6,
    units = "in", 
    res = 300)
    par(family = "serif", mar = rep(0, times = 4), bty = "n")
    plot(x = 1:10, 1:10, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    legend("topleft", 
           legend = all_coinf_combs_cols_leg[,1], 
            fill = all_coinf_combs_cols_leg[,2], 
            bty = "n")
dev.off()
###* end of figure 1 in main text *###

# 3.2 Co-occurrences
res_cooc <- readRDS(file = file.path(dirs$fits, "res_cooc.rds"))
names(res_cooc)
names(res_cooc$significants)

###* FIGURE 2 IN MAIN TEXT *###

# percentage of significant pairs our of analysed pairs (85%)
res_cooc$cooccur_full_data$co_occurrences / (res_cooc$cooccur_full_data$pot_pairs - res_cooc$cooccur_full_data$omitted)
# percentage of significant pairs our of all possible pairs (31%)
res_cooc$cooccur_full_data$co_occurrences / res_cooc$cooccur_full_data$pot_pairs

plotops <- c("populations", "all", "aggr_by_pop")
for (whichplot in plotops) {
    if (whichplot == "populations") {
        toPlot <- as.data.frame(res_cooc$significants$populations[, -1])
    }
    if (whichplot == "all") {
        toPlot <- as.data.frame(res_cooc$significants$full_data[, -1])
    }
    if (whichplot == "aggr_by_pop") {
        toPlot <- as.data.frame(res_cooc$significants$populations_aggregated[, -1])
    }

    colnames(toPlot) <- c("sp1", "sp2", "cooccs", "obs_div_exp")
    toPlot[, "obs_div_exp"] <- as.numeric(toPlot[, "obs_div_exp"])
    toPlot[, "cooccs"] <- as.numeric(toPlot[, "cooccs"])
    toPlot[, "sp1"] <- sub("sp_", "", toPlot[, "sp1"])
    toPlot[, "sp2"] <- sub("sp_", "", toPlot[, "sp2"])

    toPlot_obs <- toPlot[, c("sp1", "sp2", "cooccs")]

    pdf(file = file.path(dirs$figs,
                         paste0("cooccurrences_", whichplot, ".pdf")),
        bg = "transparent", 
        width = 6, 
        height = 6)
        par(family = "serif")
        circleplot(toPlot_obs, 
                   cluster = FALSE,
                   style = "classic",
                   plot.control = list(point.labels = TRUE,
                                       cex.point = 15,
                                       line.breaks = c(-2,0,4,14,100),
                                       line.cols = c("#ffffff",
                                                    "#f5d3dc", 
                                                     "#e795aa", 
                                                     "#c60032"),
                                       line.widths = 5))
    dev.off()
}

pdf(file = file.path(dirs$figs, 
                     "raw_coocc_legend.pdf"),
    bg = "transparent", 
    width = 5, 
    height = 3)
    par(family = "serif")
    plot(x = 1:10, 1:10, type = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("topleft", 
           legend = c("1-4 co-occurrences", 
                      "5-14 co-occurrences", 
                      ">15 co-occurrences"), 
           fill = c("#f5d3dc", 
                    "#e795aa", 
                    "#c60032"), 
           bty = "n")
dev.off()
###* end of figure 2 in main text *###


# 3.3 Markov Random Fields

# (read the final fitted models)
bootedCRF_w_mems <- readRDS(file = file.path(dirs$fits, "bootedCRF_w_mems_sampleprop100_n100.rds"))
bootedMRF <- readRDS(file = file.path(dirs$fits, "bootedMRF_sampleprop100_n100.rds"))

# 3.1 Sample level

# associations btwn viruses
MRF_cov_coefs <- plotMRF_hm(MRF_mod = bootedMRF,
                            node_names = colnames(y),
                            main = "Predicted associations (95% CIs)")
CRF_cov_coefs <- plotMRF_hm(MRF_mod = bootedCRF_w_mems,
                            node_names = colnames(y),
                            main = "Predicted associations (95% CIs)")

sig_assocs <- list("mrf" = NA, "crf" = NA)
models <- c("mrf", "crf")
for (whichmodel in models) {

    if (whichmodel == "mrf") {
        assocs <- MRF_cov_coefs$data
    }
    if (whichmodel == "crf") {
        assocs <- CRF_cov_coefs$data
    }
    ass_means <- assocs[which(assocs$Factor == "Mean"), c("value", "Var1", "Var2")]
    ass_means_xtab <- xtabs(value~Var1+Var2, data = ass_means)
    ass_means_mat <- matrix(ass_means_xtab, nrow(ass_means_xtab), ncol(ass_means_xtab))
    colnames(ass_means_mat) <- colnames(ass_means_xtab)
    rownames(ass_means_mat) <- rownames(ass_means_xtab)
    signs <- sign(assocs[which(assocs$Factor == "Lower (5%)"), "value"]) + sign(assocs[which(assocs$Factor == "Upper (95%)"), "value"])

    ass_mat <- cbind(assocs[which(assocs$Factor == "Lower (5%)"), c("value", "Var1", "Var2")],
                     assocs[which(assocs$Factor == "Upper (95%)"), c("value", "Var1", "Var2")],
                     assocs[which(assocs$Factor == "Mean"), c("value")])
    all(ass_mat[, 2] == ass_mat[, 5]) & all(ass_mat[, 3] == ass_mat[, 6])
    ass_mat <- ass_mat[, c(2,3,1,4,7)]
    colnames(ass_mat) <- c("sp1", "sp2", "5%", "95%", "Mean")
    ass_tmp <- assocs[assocs$Factor == "Mean", 1:3]
    ass_tmp$Var1 <- as.character(ass_tmp$Var1)
    ass_tmp$Var2 <- as.character(ass_tmp$Var2)
    toplot_pos <- ass_tmp[which(signs == 2),]
    toplot_pos <- toplot_pos[which(toplot_pos$value > 0),]
    toplot_neg <- ass_tmp[which(signs == -2),]
    toplot_neg <- toplot_neg[which(toplot_neg$value < 0),]
    toplot <- rbind(toplot_pos, toplot_neg)
    sig_assocs[[whichmodel]] <- toplot
}

#grid::grid.draw(CRF_cov_coefs)
#pdf(file = file.path(dirs$figs,
#                     "associations.pdf"),
#    bg = "transparent", 
#    width = 18, 
#    height = 6)
#    par(family = "serif")
#    plot(MRF_cov_coefs)
#dev.off()

###* FIGURE 3 IN MAIN TEXT *###
library(circleplot)
for (whichmodel in models) {
    pdf(file = file.path(dirs$figs,
                         paste0("associations_circle_", whichmodel, ".pdf")),
        bg = "transparent", 
        width = 10, 
        height = 10)
        par(family = "serif")
        circleplot(sig_assocs[[whichmodel]],
                   cluster = FALSE,
                   style = "classic",
                   plot.control = list(point.labels = TRUE,
                                   cex.point = 15,
                                   line.breaks = c(-10,-1,-0.5,0,0.5,1,10),
                                   line.cols = c("#004080",
                                                 "#62b1ff",
                                                 "#afd7ff",
                                                 "#f5d3dc", 
                                                 "#e795aa", 
                                                 "#c60032"),
                                   line.widths = 5))
    dev.off()
}
#par(mfrow = c(1,2))

#net <- graph.adjacency(ass_means_mat, weighted = TRUE, mode = "undirected")
#library(igraph)
#pdf(file = file.path(dirs$figs,
#                     paste0("associations_network_", whichmodel, ".pdf")),
#    bg = "transparent", 
#    width = 10, 
#    height = 10)
#    par(family = "serif")
#    igraph::plot.igraph(net, 
#                        layout = igraph::layout.circle,
#                        vertex.label.dist = 0.1,
#                        vertex.label.cex = 1.5,                    
#                        vertex.label.color = "transparent",                    
#                        vertex.size = 0,
#                        #vertex.color = "white",
#                        edge.width = abs(igraph::E(net)$weight) ^2,
#                        edge.color = ifelse(igraph::E(net)$weight < 0.1, "red2", "red3"))
#dev.off()

###* end of figure 3 in main text *###

#pdf(file = file.path(dirs$figs,
#                     "associations_circle.pdf"),
#    bg = "transparent", 
#    width = 18, 
#    height = 6)
#    par(family = "serif")
#    plot(MRF_cov_coefs)
#dev.off()

# direct covariate effects
cov_inds <- 18:33
colnames(bootedCRF_w_mems$direct_coef_means)
bootedCRF_w_mems$direct_coef_means[,cov_inds]
colSums(bootedCRF_w_mems$direct_coef_means != 0)

signs <- sign(bootedCRF_w_mems$direct_coef_upper90) + sign(bootedCRF_w_mems$direct_coef_lower90)
sigs <- which(signs == -2 | signs == 2, arr.ind = TRUE)

cbind(rownames(bootedCRF_w_mems$direct_coef_means)[sigs[,1]], 
      colnames(bootedCRF_w_mems$direct_coef_means)[sigs[,2]], 
      bootedCRF_w_mems$direct_coef_means[sigs], 
      signs[sigs])

min(bootedCRF_w_mems$direct_coef_means[,cov_inds])
max(bootedCRF_w_mems$direct_coef_means[,cov_inds])
bootedCRF_w_mems$direct_coef_means[,cov_inds]

round(apply(bootedCRF_w_mems$direct_coef_means[,cov_inds], 2, min), 3)
round(apply(bootedCRF_w_mems$direct_coef_means[,cov_inds], 2, max), 3)
dim(bootedCRF_w_mems$direct_coef_means)

graysc <- paste0("gray", round(seq(10, 90, length.out = 16)))
for (j in 1:2) {
    pdf(file = file.path(dirs$figs,
                         paste0("direct_main_effects", j, ".pdf")),
        bg = "transparent", 
        width = 10, 
        height = 6)

    yLm <- switch(j, 
                  "1" = 1.5,
                  "2" = 0.3)

    par(mar = c(10, 2, 1, 1), family = "serif")
    plot(0, 0, 
         xlim = c(0, length(cov_inds) + 1), ylim = c(-0.3, yLm), 
         xlab = "",
         ylab = "",
         xaxt = "n",
         type = "n")
    abline(h = 0)
    axis(1, labels = colnames(bootedCRF_w_mems$direct_coef_means)[cov_inds], at = 1:length(cov_inds), las = 2)
    for (i in 1:16) {
        i2 <- cov_inds[i]
        points(y = bootedCRF_w_mems$direct_coef_means[, i2], 
               x = jitter(rep(i, length(cov_inds))),
               col = "black",
               bg = graysc, 
               pch = 21)
    }
    dev.off()
}

# key coefficients
library(circleplot)
for (whichvirus in colnames(y)) {
    tmp <- data.frame(sp1 = bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]]$Variable, 
                      sp2 = rep(whichvirus, 
                            nrow(bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]])),
                      coef = bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]]$Rel_importance)
    pdf(file = file.path(dirs$figs,
                         "key_ints",
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
                                   line.breaks = c(-20,-1,0,0.1,0.5,20),
                                   line.cols = c("#004080",
                                                 "#62b1ff",
                                                 "#f5d3dc", 
                                                 "#e795aa", 
                                                 "#c60032"),
                                   line.widths = 5))
    dev.off()
}

# How many times each virus appears as a key effect for other viruses?
for (v1 in colnames(y)) {
    tmp1 <- NA
    for (v2 in colnames(y)) {
        tmp2 <- length(grep(v1, bootedCRF_w_mems$mean_key_coefs[[v2]]$Variable))
        tmp1 <- c(tmp1, tmp2)
    }
    bromo <- tmp1[-1]
    print(c(v1, sum(bromo > 0)))
}
# Bromoviridae is a key coefficient for 9 viruses, 
# the rest, max 6, mostly less, all at least 1


##########################################################################################

key_sigs <- list()
for (i in 1:ncol(y)) {
    whichvirus <- colnames(y)[i]
    low_key_signs <- sign(bootedCRF_w_mems$direct_coef_upper90[whichvirus,bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]]$Variable])
    hig_key_signs <- sign(bootedCRF_w_mems$direct_coef_lower90[whichvirus,bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]]$Variable])
    sigs <- (low_key_signs + hig_key_signs) == 2 | (low_key_signs + hig_key_signs) == -2
    #bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]]$Variable[sigs]
    
    tmp <- data.frame(var1 = bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]]$Variable[sigs], 
                      var22 = rep(whichvirus, 
                            nrow(bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]][sigs,])),
                      coef = bootedCRF_w_mems$mean_key_coefs[whichvirus][[1]]$Rel_importance[sigs])
    key_sigs[[i]] <- tmp
}


pdf(file = file.path(dirs$figs,
                     "mrf_associations.pdf"),
    bg = "transparent", 
    width = 6, 
    height = 6)
    par(family = "serif")
    plotMRF_hm(mrf1)
dev.off()

# predictions
preds <- predict_MRF(data = data_mems_df,
                     MRF_mod = bootedCRF_w_mems)

netw <- predict_MRFnetworks(data = data_mems_df, 
                            MRF_mod = bootedCRF_w_mems, 
                            #metric = "degree",
                            cutoff = 0.1)
length(netw)
plot(netw[[4]])

# 3.2 Population level
names(bootedCRF_w_mems_pop)
bootedCRF_w_mems_pop$direct_coef_means[,18:33]

# associations btwn viruses
MRF_pop_cov_coefs <- plotMRF_hm(MRF_mod = bootedCRF_w_mems_pop,
                                node_names = colnames(y),
                                main = "Predicted associations (95% CIs)")

library(igraph)
net <- igraph::graph.adjacency(crf1_pop_mem$graph, weighted = TRUE, mode = "undirected")

pdf(file = file.path(dirs$figs,
                     "associations_network.pdf"),
    bg = "transparent", 
    width = 10, 
    height = 10)
    par(family = "serif")
    igraph::plot.igraph(net, 
                        layout = igraph::layout.circle,
                        vertex.label.dist = 0.1,
                        vertex.label.cex = 1.5,                    
                        vertex.label.color = "black",                    
                        vertex.size = 0,
                        vertex.color = "white",
                        edge.width = abs(igraph::E(net)$weight) * 3,
                        edge.color = ifelse(igraph::E(net)$weight < 0, "blue", "red3"))
dev.off()

key_pop_sigs <- list()
for (i in 1:ncol(y_pop)) {
    whichvirus <- colnames(y_pop)[i]
    low_key_signs <- sign(bootedCRF_w_mems_pop$direct_coef_upper90[whichvirus,bootedCRF_w_mems_pop$mean_key_coefs[whichvirus][[1]]$Variable])
    hig_key_signs <- sign(bootedCRF_w_mems_pop$direct_coef_lower90[whichvirus,bootedCRF_w_mems_pop$mean_key_coefs[whichvirus][[1]]$Variable])
    sigs <- (low_key_signs + hig_key_signs) == 2 | (low_key_signs + hig_key_signs) == -2
    #bootedCRF_w_mems_pop$mean_key_coefs[whichvirus][[1]]$Variable[sigs]

    tmp <- data.frame(var1 = bootedCRF_w_mems_pop$mean_key_coefs[whichvirus][[1]]$Variable[sigs], 
                      var22 = rep(whichvirus, 
                            nrow(bootedCRF_w_mems_pop$mean_key_coefs[whichvirus][[1]][sigs,])),
                      coef = bootedCRF_w_mems_pop$mean_key_coefs[whichvirus][[1]]$Rel_importance[sigs])
    key_pop_sigs[[i]] <- tmp
}
