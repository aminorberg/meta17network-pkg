
rm(list = ls(all = TRUE)) ; gc()

# define the working directory again
working_dir <- "/Users/annanorb/Documents/UZH/Projects/META/meta17network-pkg"

setwd(working_dir)
library("meta17network")

dirs <- set_dirs(working_dir = working_dir)
#saveRDS(dirs, file = file.path(working_dir, "dirs.rds"))

#devtools::install_github("nicholasjclark/MRFcov")

# 0 DATA PROCESSING ######################################################################
dat <- process_data(dirs = dirs,
                    return_data = TRUE,
                    save_data = TRUE,
                    rmNAs = TRUE)
warnings()
# these warnings are related to the missing covariate data for four plants,
# which are consequently removed from the final data set for analysis (400 -> 396 plants)
names(dat)
lapply(dat, dim)

dim(dat$Y)
dim(dat$Y)[1] * dim(dat$Y)[2]
head(dat$Y)
sum(dat$Y) / (dim(dat$Y)[1] * dim(dat$Y)[2])
max(colSums(dat$Y)/nrow(dat$Y)) # max prevalence
min(colSums(dat$Y)/nrow(dat$Y)) # min prevalence
sum(rowSums(dat$Y) == 0) / nrow(dat$Y) # % empty host plants
sum(rowSums(dat$Y) == 1) / nrow(dat$Y) # % single infections
sum(rowSums(dat$Y) > 1) / nrow(dat$Y) # % multi-infections
(sum(rowSums(dat$Y) >= 5) / nrow(dat$Y)) / (sum(rowSums(dat$Y) >= 1) / nrow(dat$Y)) # % multi-infections 5 or more viruses

### 1 DATA EXPLORATION AND DESCRIPTIVE ANALYSIS ##########################################

# all the descriptive illustrations of the data structure for figures 2-3 
# are done in the stand-alone script 'coinf_plots.r'

# plot virus prevalences
# par(family = "serif", mar = c(8, 3, 1, 1))
# barplot(colSums(dat$Y)[order(colSums(dat$Y))], las = 2, ylim = c(0, 280))

# 1.1 Nestedness analysis
	# local communities are regarded as nested 
	# if they (all) could be subsets of the same community
	# (see e.g. Rynkiewicz et al 2019)
# C.scores : nestedchecker gives the number of checkerboard units, 
	# or 2x2 submatrices where both species occur once but on different sites 
	# (Stone & Roberts 1990)
library(vegan)

# host plant level
nest_plants <- nestednodf(dat$Y[rowSums(dat$Y) != 0, ])
plot(nest_plants, names = TRUE, col = "grey50", bty = "n")
# C score
oecosimu(dat$Y, nestedchecker, "r0", statistic = "C.score")

# aggregated over populations
aggr_tmp <- NA
for (i in unique(dat$X$pop)) {
    aggr_tmp <- rbind(aggr_tmp, colSums(dat$Y[which(dat$X$pop == i), ]))
}
Y_aggr_pops <- aggr_tmp[-1,]
Y_aggr_pops[which(Y_aggr_pops > 0)] <- 1
rownames(Y_aggr_pops) <- unique(dat$X$pop)
nest_pops <- nestednodf(Y_aggr_pops)
plot(nest_pops, names = TRUE, col = "grey50")
saveRDS(nest_pops, file = file.path(dirs$fits, "nest_pops.rds"))
oecosimu(Y_aggr_pops, nestedchecker, "r0", statistic = "C.score")

# C score
oecosimu(Y_aggr_pops, nestedchecker, "r0", statistic = "C.score")

# (within populations)
nests <- list()
pops <- unique(dat$X[, "pop"])
for (i in 1:length(pops)) {
    nests[[i]] <- nestednodf(dat$Y[which(dat$X[, "pop"] == pops[i]),])
}
names(nests) <- pops
#e.g.
plot(nests[["861"]], names = TRUE, col = "grey50") # most taxon-rich
plot(nests[["946"]], names = TRUE, col = "grey50") # least taxon-rich


#&&& FIGURE 2B &&&&&
png(file.path(dirs$figs, "nestedness_bypop.png"),
    height = 6, 
    width = 6, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(1,3,8,1))
    plot(nest_pops, names = TRUE, col = "grey50", bty = "l")
dev.off()
#&&&&&&&&&&&&&&&&&&&&

# 1.2 Co-occurrence analysis

# how many times do viruses occur and co-occur?
incidences <- t(dat$Y) %*% as.matrix(dat$Y)
length(which(incidences == 0)) / 2
# 18 not co-occurring pairs

t_y <- data.frame(t(dat$Y))
y_ls <- as.list(t_y)
# number of unique virus combinations 238
length(unique(y_ls))

# Coexistence and species-area curves
set.seed(7)
coexs <- coexistence_curves(Y = dat$Y, 
							xy = dat$X[, c("x", "y")], 
							sample_ids = dat$X[, "sample_ID"], 
							dirs = dirs, 
							nsim = 100)
names(coexs)
# all species were encountered when 275 plants were sampled
which(round(rowMeans(coexs[[3]]), 1) == ncol(dat$Y), arr.ind = TRUE)[1]
# all coexistences were encountered when 319 plants were sampled
which(round(rowMeans(coexs[[1]]), 1) == max(rowMeans(coexs[[1]])), arr.ind = TRUE)[1]

# where are we at when 100 plants have been sampled?
rowMeans(coexs[[1]])[100]
rowMeans(coexs[[3]])[100]
((21^2) - 21) / 2
((25^2) - 25) / 2

#&&& FIGURE 3B &&&&&
# coexistence mean curve
png(file.path(dirs$figs, "coex_mean_curve.png"),
    height = 4, 
    width = 7, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(3,1,1,3))
    plot(y = coexs[[1]][,1], x = 1:nrow(coexs[[1]]), 
         type = "n",
         lwd = 1.5, 
         ylab = "", 
         xlab = "",
         yaxt = "n",
         ylim = c(0, ((ncol(dat$Y)^2)-ncol(dat$Y))/2))
    lines(y = rowMeans(coexs[[1]]),
          x = 1:nrow(coexs[[1]]),
          lwd = 2)
    lines(y = ((rowMeans(coexs[[3]])^2) - rowMeans(coexs[[3]])) / 2,
          x = 1:nrow(coexs[[3]]),
          lwd = 2,
          col = "grey",
          lty = 2)
	axis(4, at = c(10, 50, 100, 282, 300), las = 2)
dev.off()

# species-area accumulation curve
png(file.path(dirs$figs, "sp_area_mean_curve.png"),
    height = 4, 
    width = 7, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(3,3,1,1))
    plot(y = rowMeans(coexs[[3]]), x = 1:nrow(coexs[[3]]), 
         type = "n",
         lty = 2, 
         lwd = 1.5, 
         ylab = "", 
         xlab = "",
         yaxt = "n",
         ylim = c(0, ncol(dat$Y)))
    lines(y = rowMeans(coexs[[3]]),
          x = 1:nrow(coexs[[3]]),
          lwd = 2,
          lty = 3)
	axis(2, at = c(5, 10, 15, ncol(dat$Y)), las = 2)
dev.off()
#&&&&&&&&&&&&&&&&&&&&

# all simulated coexistence curves
#curve_type <- 1
for (curve_type in 1:length(coexs)) {
	png(file.path(dirs$figs, paste0("coex", curve_type, "_all_curves.png")),
		height = 4, 
		width = 7, 
		bg = "transparent",
		units = "in", 
		res = 300)
		par(family = "serif", mar = c(3,3,1,1))
		plot(y = coexs[[curve_type]][,1], x = 1:nrow(coexs[[curve_type]]), 
			 type = "l",
			 lwd = 1.5,
			 lty = 3, 
			 ylab = "", 
			 xlab = "")
		for (i in 2:ncol(coexs[[curve_type]])) {
			lines(y = coexs[[curve_type]][,i], 
			x = 1:nrow(coexs[[curve_type]]),
			lwd = 1.5)
		}
	dev.off()
}

# 2 DATA PROCESSING FOR MRF and CRFs, Moran's eigenvectors ###############################

### Y
y_all <- as.matrix(dat$Y)
y_all_cocs <- t(y_all) %*% y_all
colSums(y_all[rowSums(y_all) == 1,]) # single infections

# subset 
sp_subset_thr <- round(nrow(dat$Y) * 0.025) # 10
sp_subset <- colnames(dat$Y)[colSums(dat$Y) >= sp_subset_thr]
dat1 <- dat
dat1$Y <- dat1$Y[, sp_subset]
y <- as.matrix(dat1$Y)
y_cocs <- t(y) %*% y
dim(y)

### X
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

# numbers of unique values
lapply(apply(x_nums_scaled, 2, unique), length)

### Spatial eigenvector maps
library(adespatial)
library(spdep)

mor_plants <- dbmem(x_spat[, c("x", "y")], 
					store.listw = TRUE,
					MEM.autocor = "positive")

# test the significance of the eigenvectors
#mor_test <- moran.randtest(mor_plants, listw = attr(mor_plants, "listw"))
mor_test <- apply(mor_plants, 2, moran.mc, listw = attr(mor_plants, "listw"), nsim = 999)
# we select three significant vectors representing different spatial scales
mor_test_sig <- mor_test[which(lapply(mor_test, "[[", "p.value") <= 0.05)]
sel_vecs <- round(seq(1, length(mor_test_sig), length.out = 4))
mems_sel <- mor_plants[, sel_vecs]
apply(mems_sel, 2, function(x) { length(unique(x)) })

# &&& FIGURE S3 &&&&&
png(file.path(dirs$figs,
              "MEMs_squares_map.png"), 
    height = 10, 
    width = 10, 
    bg = "transparent",
    units = "in", 
    res = 300)
	par(mfrow = c(2, 2), family = "serif")
	for (i in 1:ncol(mems_sel)) {
		amnt <- 1500
		set.seed(7)
		pointscols <- mems_sel[,i]
		pointscols[which(mems_sel[,i] >= 0)] <- "white"
		pointscols[which(mems_sel[,i] < 0)] <- "black"
		pointscols2 <- mems_sel[,i]
		pointscols2[which(mems_sel[,i] >= 0)] <- "black"
		pointscols2[which(mems_sel[,i] < 0)] <- "white"
		plot(x = jitter(x_spat$x, amount = amnt), 
			 y = jitter(x_spat$y, amount = amnt), 
			 col = pointscols2, 
			 bg = pointscols, 
			 cex = abs(mems_sel[,i]) + 0.25, 
			 main = paste("MEM", sel_vecs[i]),
			 xlab = "",
			 ylab = "",
			 pch = 22)
	}
dev.off()

### combine Xs
x_and_mems <- as.data.frame(cbind(x_nums_scaled, x_bools, mems_sel))

# check for correlations
var_cors <- cor(x_and_mems, use = "na.or.complete")
max(var_cors[upper.tri(var_cors, diag = FALSE)]) # highest correlation 0.51
which(var_cors == max(var_cors[upper.tri(var_cors, diag = FALSE)]), arr.ind = TRUE)

# final data frames
data_df <- as.data.frame(cbind(y, cbind(x_nums_scaled, x_bools)))
data_host_df <- data_habi_df <- data_df[, c(1:ncol(y))]
data_host_df <- cbind(data_host_df, 
					  data_df[, c("plant_size", 
					  			  "moth", 
					  			  "miner", 
					  			  "spittle", 
					  			  "holes", 
					  			  "suck_or_bite")])
data_habi_df <- cbind(data_habi_df,
					  data_df[, c("plm2",
					  			  "kytk17",
					  			  "prop_agri_area",
					  			  "shannon",
					  			  "severe_winter_days",
					  			  "temp_eff_days_summer16")]) 
data_mems_df <- as.data.frame(cbind(y, x_and_mems))
data_only_mems_df <- as.data.frame(cbind(y, mems_sel))


# 3 (CONDITIONAL) MARKOV RANDOM FIELDS, (C)MRF ###########################################
library(MRFcov)

# 3.1 Unconditional MRF
mrf1 <- MRFcov(as.data.frame(y), family = "binomial")
saveRDS(mrf1, file = file.path(dirs$fits, "mrf1.rds"))

# 3.2 Conditional MRF: non-spatial, spatial, spatial using Moran's eigenvectors,
# the last one with and without other environmental covariates

# 3.2.1 Non-spatial CRF

### with all environmental variables
crf1 <- MRFcov(data_df, family = "binomial", n_nodes = ncol(y))
saveRDS(crf1, file = file.path(dirs$fits, "crf1.rds"))

### with only host variables
crf1_host <- MRFcov(data_host_df, family = "binomial", n_nodes = ncol(y))
saveRDS(crf1_host, file = file.path(dirs$fits, "crf1_host.rds"))

### with only habitat variables
crf1_habi <- MRFcov(data_habi_df, family = "binomial", n_nodes = ncol(y))
saveRDS(crf1_habi, file = file.path(dirs$fits, "crf1_habi.rds"))

# 3.2.2 Spatial CRF with Morans' eigenvectors

### only MEMs
crf1_only_mem <- MRFcov(data_only_mems_df, family = "binomial", n_nodes = ncol(y))
saveRDS(crf1_only_mem, file = file.path(dirs$fits, "crf1_only_mem.rds"))

### with environment
crf1_mem <- MRFcov(data_mems_df, family = "binomial", n_nodes = ncol(y))
saveRDS(crf1_mem, file = file.path(dirs$fits, "crf1_mem.rds"))

# 3.3 Predictions
# (load the previously fitted models)
mrf1 <- readRDS(file = file.path(dirs$fits, "mrf1.rds"))
crf1 <- readRDS(file = file.path(dirs$fits, "crf1.rds"))
crf1_host <- readRDS(file = file.path(dirs$fits, "crf1_host.rds"))
crf1_habi <- readRDS(file = file.path(dirs$fits, "crf1_habi.rds"))
crf1_only_mem <- readRDS(file = file.path(dirs$fits, "crf1_only_mem.rds"))
crf1_mem <- readRDS(file = file.path(dirs$fits, "crf1_mem.rds"))

# make predictions
preds_mrf1 <- predict_MRF(data_df, mrf1)
preds_crf1 <- predict_MRF(data_df, crf1)
preds_crf1_host <- predict_MRF(data_host_df, crf1_host)
preds_crf1_habi <- predict_MRF(data_habi_df, crf1_habi)
preds_crf1_only_mem <- predict_MRF(data_only_mems_df, crf1_only_mem)
preds_crf1_mem <- predict_MRF(data_mems_df, crf1_mem)

# 3.4 Model fit evaluation

# 3.4.1 AUCs for sample level
aucs <- list("mrf" = NA, 
             "crf_env" = NA, 
             "crf_host" = NA, 
             "crf_habi" = NA, 
             "crf_only_mem" = NA, 
             "crf_env_mem" = NA)
aucs$mrf <- pROC::auc(response = matrix(y, ncol = 1), 
                      predictor = matrix(preds_mrf1$Probability_predictions, ncol = 1))
aucs$crf_env <- pROC::auc(response = matrix(y, ncol = 1),
                          predictor = matrix(preds_crf1$Probability_predictions, ncol = 1))
aucs$crf_env_mem <- pROC::auc(response = matrix(y, ncol = 1),       # seems to be the best (marginally)
                              predictor = matrix(preds_crf1_mem$Probability_predictions, ncol = 1))
aucs$crf_only_mem <- pROC::auc(response = matrix(y, ncol = 1),  
                               predictor = matrix(preds_crf1_only_mem$Probability_predictions, ncol = 1))
aucs$crf_host <- pROC::auc(response = matrix(y, ncol = 1),  
                           predictor = matrix(preds_crf1_host$Probability_predictions, ncol = 1))
aucs$crf_habi <- pROC::auc(response = matrix(y, ncol = 1),  
                           predictor = matrix(preds_crf1_habi$Probability_predictions, ncol = 1))
aucs
# CRF is better than just MRF, but there is no big difference
# whether one includes spatial predictors and/or landscape and/or host variables

# 3.4.2 Cross-validation and bootstrapping
nFolds <- 10
nCores <- 4

# environment and MEMs
crf_env_mems_cv <- lapply(seq_len(100), 
                          function(x) { cv_MRF_diag(data = data_mems_df, 
                                                    n_nodes = ncol(y), 
                                                    n_folds = nFolds, 
                                                    n_cores = nCores, 
                                                    family = "binomial", 
                                                    compare_null = FALSE, 
                                                    plot = FALSE, 
                                                    cached_model = crf1_mem, 
                                                    cached_predictions = list(predictions = preds_crf1_mem), 
                                                    sample_seed = 10)
                                      })
crf_env_and_mems_cv <- do.call(rbind, crf_env_mems_cv)

# only MEMs
crf_only_mems_cv <- lapply(seq_len(100), 
                           function(x) { cv_MRF_diag(data = data_only_mems_df, 
                                                     n_nodes = ncol(y), 
                                                     n_folds = nFolds, 
                                                     n_cores = nCores, 
                                                     family = "binomial", 
                                                     compare_null = FALSE, 
                                                     plot = FALSE, 
                                                     cached_model = crf1_only_mem, 
                                                     cached_predictions = list(predictions = preds_crf1_only_mem), 
                                                     sample_seed = 10)
                                       })
crf_only_mems_cv <- do.call(rbind, crf_only_mems_cv)

# only env
crf_only_env_cv <- lapply(seq_len(100), 
                         function(x) { cv_MRF_diag(data = data_df, 
                                                   n_nodes = ncol(y), 
                                                   n_folds = nFolds, 
                                                   n_cores = nCores, 
                                                   family = "binomial", 
                                                   compare_null = FALSE, 
                                                   plot = FALSE, 
                                                   cached_model = crf1, 
                                                   cached_predictions = list(predictions = preds_crf1), 
                                                   sample_seed = 10)
                                     })
crf_only_env_cv <- do.call(rbind, crf_only_env_cv)

# only host
crf_only_host_cv <- lapply(seq_len(100), 
                           function(x) { cv_MRF_diag(data = data_host_df, 
                                                     n_nodes = ncol(y), 
                                                     n_folds = nFolds, 
                                                     n_cores = nCores, 
                                                     family = "binomial", 
                                                     compare_null = FALSE, 
                                                     plot = FALSE, 
                                                     cached_model = crf1_host, 
                                                     cached_predictions = list(predictions = preds_crf1_host), 
                                                     sample_seed = 10)
                                        })
crf_only_host_cv <- do.call(rbind, crf_only_host_cv)

# only habi
crf_only_habi_cv <- lapply(seq_len(100), 
                           function(x) { cv_MRF_diag(data = data_habi_df, 
                                                     n_nodes = ncol(y), 
                                                     n_folds = nFolds, 
                                                     n_cores = nCores, 
                                                     family = "binomial", 
                                                     compare_null = FALSE, 
                                                     plot = FALSE, 
                                                     cached_model = crf1_habi, 
                                                     cached_predictions = list(predictions = preds_crf1_habi), 
                                                     sample_seed = 10)
                                        })
crf_only_habi_cv <- do.call(rbind, crf_only_habi_cv)

# MRF
mrf_cv <- lapply(seq_len(100), 
                 function(x) { cv_MRF_diag(data = data_df, 
										   n_nodes = ncol(y), 
										   n_folds = nFolds, 
										   n_cores = nCores, 
                                           family = "binomial", 
                                           compare_null = FALSE, 
                                           plot = FALSE, 
                                           cached_model = crf1_only_mem, 
                                           cached_predictions = list(predictions = preds_mrf1), 
                                           sample_seed = 10)
                           })
mrf_cv <- do.call(rbind, mrf_cv)

cv_res <- list("env_and_mems" = crf_env_and_mems_cv, 
               "only_mems" = crf_only_mems_cv,
               "only_env" = crf_only_env_cv,
               "only_host" = crf_only_host_cv,
               "only_habi" = crf_only_habi_cv,
               "mrf" = mrf_cv)
saveRDS(cv_res, file = file.path(dirs$fits, "cv_res.rds"))

cv_res <- readRDS(file = file.path(dirs$fits, "cv_res.rds"))

lapply(cv_res, apply, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
#pos_pred = the proportion of true positives out of all predicted positives
#sensitivity = prop. of true positives out of all real positives
    # i.e. probability of a positive prediction, given that the virus is present
#neg_pred = prop. of true negatives out of all predicted negatives
#specificity = prop. of true negatives out of all real negatives 
    # i.e. probability of a negative prediction, given that the virus is absent
#tot_pred = prop. of all true predictions out of the whole data

# MRFs predicts more false positives than false negatives
# higher sensitivity than CRFs, but lower specificity

# CRFs have very little difference among each other,
# in general a high proportion of true positives out of all predicted positives
# lowish sensitivity, so a fair amount of false negatives
# high specificity, so very few false positives

# 3.5 Fitting final CRFs with bootstrapping (for obtaining parameter uncertainties)
#library(MRFcov)

sampleProp <- 1
nBoot <- 100
nodes <- ncol(y)

# MRF
bootedMRF <- bootstrap_MRF(as.data.frame(y),
                           n_nodes = nodes,
                           family = "binomial",
                           sample_prop = sampleProp,
                           n_bootstraps = nBoot)
fileName <- paste0("bootedMRF", "_sampleprop", sampleProp, "_boots", nBoot, ".rds")
saveRDS(bootedMRF, 
        file = file.path(dirs$fits, fileName))

# CRF with env and mems
bootedCRF_w_env_mems <- bootstrap_MRF(data_mems_df,
                                      n_nodes = nodes,
                                      family = "binomial",
                                      sample_prop = sampleProp,
                                      n_bootstraps = nBoot)
fileName <- paste0("bootedCRF_w_env_mems", "_sampleprop", sampleProp, "_boots", nBoot, ".rds")
saveRDS(bootedCRF_w_env_mems, 
        file = file.path(dirs$fits, fileName))

# CRF with mems only
bootedCRF_w_mems_only <- bootstrap_MRF(data_only_mems_df,
                                       n_nodes = nodes,
                                       family = "binomial",
                                       sample_prop = sampleProp,
                                       n_bootstraps = nBoot)
fileName <- paste0("bootedCRF_w_mems_only", "_sampleprop", sampleProp, "_boots", nBoot, ".rds")
saveRDS(bootedCRF_w_mems_only, 
        file = file.path(dirs$fits, fileName))

# CRF with env
bootedCRF_w_env <- bootstrap_MRF(data_df,
                                 n_nodes = nodes,
                                 family = "binomial",
                                 sample_prop = sampleProp,
                                 n_bootstraps = nBoot)
fileName <- paste0("bootedCRF_w_env", "_sampleprop", sampleProp, "_boots", nBoot, ".rds")
saveRDS(bootedCRF_w_env, 
        file = file.path(dirs$fits, fileName))

# CRF with habitat
bootedCRF_w_habi <- bootstrap_MRF(data_habi_df,
                                  n_nodes = nodes,
                                  family = "binomial",
                                  sample_prop = sampleProp,
                                  n_bootstraps = nBoot)
fileName <- paste0("bootedCRF_w_habi", "_sampleprop", sampleProp, "_boots", nBoot, ".rds")
saveRDS(bootedCRF_w_habi, 
        file = file.path(dirs$fits, fileName))

# CRF with host vars
bootedCRF_w_host <- bootstrap_MRF(data_host_df,
                                  n_nodes = nodes,
                                  family = "binomial",
                                  sample_prop = sampleProp,
                                  n_bootstraps = nBoot)
fileName <- paste0("bootedCRF_w_host", "_sampleprop", sampleProp, "_boots", nBoot, ".rds")
saveRDS(bootedCRF_w_host, 
        file = file.path(dirs$fits, fileName))

fileBody <- paste0("_sampleprop", sampleProp, "_boots", nBoot, ".rds")

# bootedMRF <- readRDS(file = file.path(dirs$fits, paste0("bootedMRF", fileBody)))
# bootedCRF_w_env_mems <- readRDS(file = file.path(dirs$fits, 
# 												 paste0("bootedCRF_w_env_mems", fileBody)))
# bootedCRF_w_mems_only <- readRDS(file = file.path(dirs$fits, 
# 												  paste0("bootedCRF_w_mems_only", 
# 												  		 fileBody)))
# bootedCRF_w_env <- readRDS(file = file.path(dirs$fits, 
# 											paste0("bootedCRF_w_env", fileBody)))
# bootedCRF_w_habi <- readRDS(file = file.path(dirs$fits, 
# 											 paste0("bootedCRF_w_habi", fileBody)))
# bootedCRF_w_host <- readRDS(file = file.path(dirs$fits, 
# 											 paste0("bootedCRF_w_host", fileBody)))

booted_models <- list("MRF" = bootedMRF, 
                      "CRF_env_mems" = bootedCRF_w_env_mems, 
                      "CRF_env" = bootedCRF_w_env, 
                      "CRF_mems" = bootedCRF_w_mems_only,
                      "CRF_habi" = bootedCRF_w_habi, 
                      "CRF_host" = bootedCRF_w_host)
saveRDS(booted_models, file = file.path(dirs$fits, paste0("booted_models", fileBody)))

booted_models <- readRDS(file = file.path(dirs$fits, paste0("booted_models", fileBody)))

tmp <- plotMRF_hm(MRF_mod = booted_models$CRF_env_mems, 
                 #node_names, 
                 #main, 
                 plot_observed_vals = FALSE, 
                 data = data_mems_df)

# association significances (does the 90% interval overlap with zero?)
significances <- list()
for (i in 1:length(booted_models)) {
    significances[[i]] <- (sign(booted_models[[i]]$direct_coef_upper90) + sign(booted_models[[i]]$direct_coef_lower90))[, 2:(nodes + 1)]
}
names(significances) <- names(booted_models)

# mean associations
associations <- list()
for (i in 1:length(booted_models)) {
    associations[[i]] <- booted_models[[i]]$direct_coef_means[, 2:(nodes + 1)]
}
names(associations) <- names(booted_models)
# round to one decimal to get rid of very small values
associations <- lapply(associations, round, digits = 3)

# leave significant associations, others -> 0
associations_sig <- associations
for (i in 1:length(associations_sig)) {
    associations_sig[[i]][which(abs(significances[[i]]) != 2)] <- 0
}

# check are there conflicting associations 
for (i in 1:length(associations_sig)) {
    asss <- associations_sig[[i]]
    print(all(asss[which((asss != 0) + (t(asss) != 0) == 2)] == t(asss)[which((asss != 0) + (t(asss) != 0) == 2)]))
}

# check if the association matrix is symmetric
ass_sig_symm <- associations_sig
for (i in 1:length(ass_sig_symm)) {
    for (k in 1:nrow(ass_sig_symm[[i]])) {
        for (j in 1:ncol(ass_sig_symm[[i]])) {
            if (ass_sig_symm[[i]][k,j] != 0 | ass_sig_symm[[i]][j,k] != 0) {    
                if (ass_sig_symm[[i]][k,j] != ass_sig_symm[[i]][j,k]) {
                    print(c(i,
                    		ass_sig_symm[[i]][k,j],
                            ass_sig_symm[[i]][j,k]))
#                     ass_sig_symm[[i]][k,j] <- ass_sig_symm[[i]][k,j] + ass_sig_symm[[i]][j,k]
#                     ass_sig_symm[[i]][j,k] <- ass_sig_symm[[i]][k,j]
                }
            }

        }
    }
}

# directions of the significant associations
assoc_signs <- list()
for (i in 1:length(ass_sig_symm)) {
    assoc_signs[[i]] <- sign(ass_sig_symm[[i]])
}
names(assoc_signs) <- names(ass_sig_symm)
# only MRF results in direct negative virus-virus associations

# no changes in the directions of significant associations
all(assoc_signs[[6]][which(assoc_signs[[1]] > 0, arr.ind = TRUE)] >= 0)
all(assoc_signs[[6]][which(assoc_signs[[1]] < 0, arr.ind = TRUE)] <= 0)

ass_sig_symm_cleared <- ass_sig_symm
for (i in 1:length(ass_sig_symm_cleared)) {
    ass_sig_symm_cleared[[i]][upper.tri(ass_sig_symm_cleared[[i]])] <- 0
}

# how many associations in total
nonzero_assocs <- lapply(lapply(ass_sig_symm_cleared, sign), abs)
lapply(nonzero_assocs, sum)

assoc_sum <- 0
for (i in 1:length(nonzero_assocs)) {
    assoc_sum <- nonzero_assocs[[i]] + assoc_sum
}
which(assoc_sum == 6, arr.ind = TRUE)

perm_assoc <- ass_sig_symm$CRF_env_mems
perm_assoc[which(assoc_sum != 6, arr.ind = TRUE)] <- 0
perm_assoc_cleared <- perm_assoc
perm_assoc_cleared[upper.tri(perm_assoc_cleared)] <- 0
sum(perm_assoc_cleared > 0)

# comparison to phylogeny
phyl <- ape::keep.tip(dat$phylogeny, colnames(nonzero_assocs[[1]]))
#plot(igraph::as.igraph(phyl))
phylmat <- cophenetic(phyl)
phyl_comp <- list()
for (i in 1:length(ass_sig_symm)) {
    phyl_comp[[i]] <- vegan::mantel(xdis = dist(ass_sig_symm[[i]]), 
                                    ydis = dist(phylmat))
    print(phyl_comp[[i]]$signif)
}
names(phyl_comp) <- names(ass_sig_symm)
vegan::mantel(xdis = dist(perm_assoc), ydis = dist(phylmat))
# no connection between phylogeny and associations

# create adjacency matrices for plotting the networks
ass_for_adj <- ass_sig_symm_cleared
for (i in 1:length(ass_for_adj)) {
	rownames(ass_for_adj[[i]]) <- colnames(ass_for_adj[[i]]) <- sub(pattern = "viridae", 
																	replacement = "", 
																	colnames(ass_for_adj[[i]]))
	rownames(ass_for_adj[[i]]) <- colnames(ass_for_adj[[i]]) <- sub(pattern = "viroidae", 
																	replacement = "", 
																	colnames(ass_for_adj[[i]]))
}
rownames(perm_assoc) <- colnames(perm_assoc) <- sub(pattern = "viridae", 
													replacement = "", 
													colnames(perm_assoc))
rownames(perm_assoc) <- colnames(perm_assoc) <- sub(pattern = "viroidae", 
													replacement = "", 
													colnames(perm_assoc))
								  
# create adjacency matrices and delete edges that represent very weak interactions
adj_mats <- list()
for (i in 1:length(ass_for_adj)) {
    adj_mats[[i]] <- igraph::graph.adjacency(ass_for_adj[[i]], 
                                             weighted = TRUE, 
                                             mode = "undirected")
    adj_mats[[i]] <- igraph::delete.edges(adj_mats[[i]], 
                                          which(abs(igraph::E(adj_mats[[i]])$weight) <= 0.05))
}
adj_mats[[length(adj_mats) + 1]] <- igraph::graph.adjacency(perm_assoc, 
                                                            weighted = TRUE, 
                                                            mode = "undirected")
adj_mats[[length(adj_mats)]] <- igraph::delete.edges(adj_mats[[length(adj_mats)]], 
                                                     which(abs(igraph::E(adj_mats[[length(adj_mats)]])$weight) <= 0.05))
names(adj_mats) <- c(names(ass_for_adj), "permanent") 

cols <- c(grDevices::adjustcolor("blue2", alpha.f = 0.95), 
          grDevices::adjustcolor("red", alpha.f = 0.95))

for (i in 1:length(adj_mats)) {
    igraph::E(adj_mats[[i]])$color <- ifelse(igraph::E(adj_mats[[i]])$weight < 0, 
                                             cols[1], 
                                             cols[2])
    igraph::E(adj_mats[[i]])$width <- abs(igraph::E(adj_mats[[i]])$weight) * 5
    adj_mats[[i]] <- igraph::delete.vertices(adj_mats[[i]], 
                                             igraph::degree(adj_mats[[i]]) == 0)
}

vert_lab_deg <- list()
for (i in 1:length(adj_mats)) {
	vert_lab_deg[[i]] <- rep(180, times = length(igraph::vertex_attr(adj_mats[[i]])$name))
}
names(vert_lab_deg) <- names(adj_mats)

# plotting settings
vert_lab_deg[["MRF"]][which(igraph::vertex_attr(adj_mats[["MRF"]])$name == "Alphaflexi")] <- pi/2

vert_lab_deg[["CRF_env_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_env_mems"]])$name == "Pospi")] <- pi/2
vert_lab_deg[["CRF_env_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_env_mems"]])$name == "Seco")] <- pi/4
vert_lab_deg[["CRF_env_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_env_mems"]])$name == "Gemini")] <- pi/2
vert_lab_deg[["CRF_env_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_env_mems"]])$name == "Bromo")] <- pi/2
vert_lab_deg[["CRF_env_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_env_mems"]])$name == "Metaxy")] <- pi/2

vert_lab_deg[["permanent"]][which(igraph::vertex_attr(adj_mats[["permanent"]])$name == "Bromo")] <- pi/2
vert_lab_deg[["permanent"]][which(igraph::vertex_attr(adj_mats[["permanent"]])$name == "Alphaflexi")] <- pi/2
vert_lab_deg[["permanent"]][which(igraph::vertex_attr(adj_mats[["permanent"]])$name == "Metaxy")] <- pi/2
vert_lab_deg[["permanent"]][which(igraph::vertex_attr(adj_mats[["permanent"]])$name == "Avsun")] <- pi/2
vert_lab_deg[["permanent"]][which(igraph::vertex_attr(adj_mats[["permanent"]])$name == "Rhabdo")] <- pi/2
# vert_lab_deg[["CRF_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_mems"]])$name == "Bromo")] <- 90

names(adj_mats)
#&&& FIGURE 4A-B, F in main text &&&&&
# create subdirectory for igraphs
dir.create(file.path(dirs$figs, "igraphs"))

for (i in 1:length(adj_mats)) {
    png(file.path(dirs$figs, 
    			  "igraphs",
                  paste0("igraph_", names(adj_mats)[i], ".png")), 
        height = 8, 
        width = 8, 
        bg = "transparent",
        units = "in", 
        res = 500)
        par(mar = c(1, 4, 1, 10), family = "serif")
        set.seed(77)
        plot(adj_mats[[i]], 
             edge.curved = 0.4,
             vertex.size = 6,
             vertex.color = "grey80",
             vertex.label.color = "black",
			 vertex.label.cex = 1.75,
             vertex.label.font = 2,
             vertex.frame.color = grDevices::adjustcolor("black", alpha.f = 0.8),
			 vertex.label.dist = 1.1,
			 vertex.label.degree = vert_lab_deg[[i]],
             #layout = igraph::layout.circle,
             label = TRUE)
    dev.off()
}
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# Compare associations, where can they accounted to
library(NetComp)

# returns a matrix containing edges present in matrix1 that are not present in matrix2
host_effect <- netDiff(matrix1 = ass_sig_symm_cleared$CRF_habi, 
                       matrix2 = ass_sig_symm_cleared$CRF_env_mems, 
                       cutoff = 0.05)
habi_effect <- netDiff(matrix1 = ass_sig_symm_cleared$CRF_host, 
                       matrix2 = ass_sig_symm_cleared$CRF_env_mems, 
                       cutoff = 0.05)
env_effect <- netDiff(matrix1 = ass_sig_symm_cleared$CRF_mems, 
                      matrix2 = ass_sig_symm_cleared$CRF_env_mems, 
                      cutoff = 0.05)
mems_effect <- netDiff(matrix1 = ass_sig_symm_cleared$CRF_env, 
                       matrix2 = ass_sig_symm_cleared$CRF_env_mems, 
                       cutoff = 0.05)
full_crf_effect <- netDiff(matrix1 = ass_sig_symm_cleared$MRF, 
                           matrix2 = ass_sig_symm_cleared$CRF_env_mems, 
                           cutoff = 0.05)
adj_mats_diffs <- list("mems" = mems_effect, 
					   "host" = host_effect, 
					   "habi" = habi_effect, 
					   "env" = env_effect)

for (i in 1:length(adj_mats_diffs)) {	
	colnames(adj_mats_diffs[[i]]) <- sub(pattern = "viridae", 
										 replacement = "", 
										 colnames(adj_mats_diffs[[i]]))
	colnames(adj_mats_diffs[[i]]) <- sub(pattern = "viroidae", 
										 replacement = "", 
										 colnames(adj_mats_diffs[[i]]))
	rownames(adj_mats_diffs[[i]]) <- sub(pattern = "viridae", 
										 replacement = "", 
										 rownames(adj_mats_diffs[[i]]))
	rownames(adj_mats_diffs[[i]]) <- sub(pattern = "viroidae", 
										 replacement = "", 
										 rownames(adj_mats_diffs[[i]]))
	adj_mats_diffs[[i]] <- igraph::graph.adjacency(adj_mats_diffs[[i]], 
												   weighted = TRUE, 
												   mode = "undirected")
}

# create adjacency matrices for plotting the networks
for (i in 1:length(adj_mats_diffs)) {
    # delete edges that represent weak interactions
    adj_mats_diffs[[i]] <- igraph::delete.edges(adj_mats_diffs[[i]], 
                                                which(abs(igraph::E(adj_mats_diffs[[i]])$weight) <= 0.05))
}
cols <- c(grDevices::adjustcolor("blue2", alpha.f = 0.95), 
          grDevices::adjustcolor("red", alpha.f = 0.95))

for (i in 1:length(adj_mats_diffs)) {
    igraph::E(adj_mats_diffs[[i]])$color <- ifelse(igraph::E(adj_mats_diffs[[i]])$weight < 0, 
                                                   cols[1], 
                                                   cols[2])
    igraph::E(adj_mats_diffs[[i]])$width <- abs(igraph::E(adj_mats_diffs[[i]])$weight) * 5
    adj_mats_diffs[[i]] <- igraph::delete.vertices(adj_mats_diffs[[i]], 
                                                   igraph::degree(adj_mats_diffs[[i]]) == 0)
}

vert_lab_deg <- list()
for (i in 1:length(adj_mats_diffs)) {
	vert_lab_deg[[i]] <- rep(180, 
							 times = length(igraph::vertex_attr(adj_mats_diffs[[i]])$name))
}
names(vert_lab_deg) <- names(adj_mats_diffs)

# plotting settings
#vert_lab_deg[["env"]][which(igraph::vertex_attr(adj_mats_diffs[["env"]])$name == "Solemo")] <- pi/2

names(adj_mats_diffs)

#&&& FIGURE 4C-E in main text &&&&&
for (i in 1:length(adj_mats_diffs)) {
    png(file.path(dirs$figs,
    			  "igraphs",
                  paste0("igraph_", names(adj_mats_diffs)[i], "_effect.png")), 
        height = 8, 
        width = 8, 
        bg = "transparent",
        units = "in", 
        res = 500)
        par(mar = c(1, 4, 1, 10), family = "serif")
        set.seed(77)
        plot(adj_mats_diffs[[i]], 
             edge.curved = 0.4,
             edge.lty = 2,
             vertex.size = 6,
             vertex.color = "grey80",
             vertex.label.color = "black",
			 vertex.label.cex = 1.75,
             vertex.label.font = 2,
             vertex.frame.color = grDevices::adjustcolor("black", alpha.f = 0.8), 
			 vertex.label.dist = 1.1,
			 vertex.label.degree = vert_lab_deg[[i]],
             #layout = igraph::layout.circle,
             label = TRUE)
    dev.off()
}
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# Direct covariate effects
mod <- booted_models$CRF_env_mems

names(mod)
colnames(mod$direct_coef_means)
cov_inds <- c(1, (nodes + 2):(nodes + 17))
colnames(mod$direct_coef_means)[cov_inds]
mod$direct_coef_means[, cov_inds]

signs <- sign(mod$direct_coef_upper90) + sign(mod$direct_coef_lower90)
sigs <- which(signs == -2 | signs == 2, arr.ind = TRUE)

cbind(rownames(mod$direct_coef_means)[sigs[,1]], 
      colnames(mod$direct_coef_means)[sigs[,2]], 
      mod$direct_coef_means[sigs], 
      signs[sigs])

sig_direct_env_effects <- cbind(rownames(mod$direct_coef_means)[sigs[,1]], 
                                colnames(mod$direct_coef_means)[sigs[,2]], 
                                round(mod$direct_coef_means[sigs], 2))

write.table(sig_direct_env_effects,
            file = file.path(dirs$fits,"direct_sig_env_effs.csv"),
            quote = FALSE,
            sep = ";",
            dec = ".")

##########################################################################################
