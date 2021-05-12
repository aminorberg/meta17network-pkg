
rm(list = ls(all = TRUE)) ; gc()
working_dir <- "/Users/annanorberg/switchdrive/UZH/Projects/META/meta17network-pkg"
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

max(colSums(dat$Y)/nrow(dat$Y)) # max prevalence 33%
min(colSums(dat$Y)/nrow(dat$Y)) # min prevalence 0.5%
sum(rowSums(dat$Y) == 0) / nrow(dat$Y) # 29% empty host plants, 71% infected plants
sum(rowSums(dat$Y) == 1) / nrow(dat$Y) # 32% single infections
sum(rowSums(dat$Y) > 1) / nrow(dat$Y) # 40% multi-infections

### 1 DATA EXPLORATION AND DESCRIPTIVE ANALYSIS ##########################################

# plot virus prevalences
par(family = "serif", mar = c(8, 3, 1, 1))
barplot(colSums(dat$Y)[order(colSums(dat$Y))], las = 2, ylim = c(0, 150))

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
plot(nests[["861"]], names = TRUE, col = "grey50")
plot(nests[["3222"]], names = TRUE, col = "grey50")


#&&& FIGURE 2B &&&&&
png(file.path(dirs$figs, "nestedness_bypop.png"),
    height = 6, 
    width = 6, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(1,3,8,1))
    plot(nest_pops, names = TRUE, col = "grey50", bty = "n")
dev.off()
#&&&&&&&&&&&&&&&&&&&&

# 1.2 Co-occurrence illustration and analysis

# how many times do viruses occur and co-occur?
incidences <- t(dat$Y) %*% as.matrix(dat$Y)
incidences["Avsunviroidae","Alphasatellitidae"] # the only pair not co-occurring together
# hence the maximum amount of co-occurring pairs is 299

t_y <- data.frame(t(dat$Y))
y_ls <- as.list(t_y)
# number of unique virus combinations
length(unique(y_ls))

#&&& FIGURE 2C-D #&&&&&
library(corrplot)
virus_colrs <- load_colour_palette()[[1]][1:25]
# based on species richness and population nestedness, 
# the two most contrasting populations are 861 and 3222
wpop <- "861"
#wpop <- "3222"
toPLot <- dat$Y[which(dat$X$pop == wpop),]
png(file.path(dirs$figs, paste0("viruses_by_plant_pop", wpop, ".png")),
    height = 4, 
    width = 6, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(1,2,1,1))
    barplot(t(toPLot[order(rowSums(toPLot), decreasing = TRUE),]),
            col = virus_colrs,
            ylim = c(0, 25),
            xaxt = "n",
            las = 2)
dev.off()

png(file.path(dirs$figs, "viruses_by_plant_pop_legend.png"),
    height = 8, 
    width = 4, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(1,2,1,1))
    plot(rep(1, 10), y = 1:10, type = "n", xaxt = "n", yaxt = "n", bty = "n")
    legend(x = 1, y = 10, legend = colnames(dat$Y), fill = virus_colrs)
dev.off()
#&&&&&&&&&&&&&&&&&&&&

# (virus (co)incidences at the contrasting populations)
wpop <- "861"
#wpop <- "3222"
incids <- t(dat$Y[which(dat$X$pop == wpop),]) %*% as.matrix(dat$Y[which(dat$X$pop == wpop),])
png(file.path(dirs$figs, paste0("co_occs_pop", wpop, ".png")),
    height = 10, 
    width = 8, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(0,5,5,5))
    corrplot(incids, 
			 is.corr = FALSE,
			 method = "number",
			 type = "lower",
			 order = "AOE",
			 tl.col = "black",
			 cl.pos = "n",
			 addgrid.col = "black",
			 col = colorRampPalette(c("white","black"))(100))
dev.off()

# (coinfections per population)
coinfs <- coinfs_by_pop(Y = dat$Y, 
                        POPs = dat$X$pop, 
                        dirs = dirs, 
                        pop_lifes = list("861" = NA, "3225" = NA))
saveRDS(coinfs, file = file.path(dirs$fits, "coinfs.rds"))

# (coinfections in the full data)
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

# Coexistence and species-area curves
coexs <- coexistence_curves(Y = dat$Y, 
							xy = dat$X[, c("x", "y")], 
							sample_ids = dat$X[, "sampleID"], 
							dirs = dirs, 
							nsim = 100)

#&&& FIGURE 3B &&&&&
# coexistence mean curves
png(file.path(dirs$figs, "coex_mean_curves.png"),
    height = 4, 
    width = 7, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(3,3,1,1))
    plot(y = coexs[[1]][,1], x = 1:nrow(coexs[[1]]), 
         type = "n",
         lwd = 1.5, 
         ylab = "", 
         xlab = "")
    lines(y = rowMeans(coexs[[1]]),
          x = 1:nrow(coexs[[1]]),
          lwd = 2)
    lines(y = rowMeans(coexs[[2]]),
          x = 1:nrow(coexs[[2]]),
          lty = 2,
          lwd = 2)
dev.off()

# species-area accumulation curve
nsp <- rowMeans(sp_area)[100]
((round(ncol(dat$Y))^2) - round(ncol(dat$Y)))/2
rowMeans(insidSumS)[100]
rowMeans(insidSumS)[387]
quartz()
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
         xlab = "")
    lines(y = rowMeans(coexs[[3]]),
          x = 1:nrow(coexs[[3]]),
          lwd = 2,
          lty = 2)
dev.off()
#&&&&&&&&&&&&&&&&&&&&

# all simulated coexistence curves
curve_type <- 1
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


# 2 DATA PROCESSING FOR MRF and CRFs, Moran's eigenvectors ###############################

### Y
y_all <- as.matrix(dat$Y)
y_all_cocs <- t(y_all) %*% y_all
#y_all_cocs[upper.tri(y_all_cocs)] <- NA
#diag(y_all_cocs) <- NA
#which(y_all_cocs == 0, arr.ind = TRUE)
# only one pair never coexists: Avsunviroidae and Alphasatellitidae
colSums(y_all[rowSums(y_all) == 1,]) # single infections

# subset 
sp_subset_thr <- round(nrow(dat$Y) * 0.025) # 10
sp_subset <- colnames(dat$Y)[colSums(dat$Y) >= sp_subset_thr]
dat1 <- dat
dat1$Y <- dat1$Y[, sp_subset]
y <- as.matrix(dat1$Y)
y_cocs <- t(y) %*% y
#y_cocs[upper.tri(y_cocs)] <- NA
#diag(y_cocs) <- NA

#&&& FIGURE 3A &&&&&
library(corrplot)
png(file.path(dirs$figs, "co_occs_all.png"),
    height = 11, 
    width = 9, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(0,5,5,5))
    corrplot(y_all_cocs,
             is.corr = FALSE,
             method = "shade",
             type = "lower",
             order = "FPC",
             tl.col = "black",
             cl.pos = "n",
             addgrid.col = "black",
			 col = colorRampPalette(c(rep("#FFFFFF", 4), 
									  "#D1E5F0",
									  "#92C5DE", 
									  "#4393C3", 
									  "#2166AC"))(200),
             bg = "transparent")
    corrplot(y_all_cocs, 
             is.corr = FALSE,
			 add = TRUE,
             method = "number",
             type = "lower",
             order = "FPC",
             tl.col = "black",
             cl.pos = "n",
             addgrid.col = "black",
             bg = "transparent",
             col = "black")
dev.off()
#&&&&&&&&&&&&&&&&&&&&

# (coinfections in the subsetted data)
y_coinfs <- y
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
y_coinfs <- lapply(tmp1, sort)
y_coinfs <- table(sort(unlist(y_coinfs)))
y_coinfs <- all_coinfs[order(y_coinfs)]
#saveRDS(y_coinfs, file = file.path(dirs$fits, "y_coinfs.rds"))
par(mar = c(70, 3, 1, 1), family = "serif")
barplot(y_coinfs, las = 2)


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
library(adegraphics)
library(spdep)

mor_plants <- dbmem(x_spat[, c("x", "y")], 
					store.listw = TRUE,
					MEM.autocor = "positive")

# test the significance of the eigenvectors
#mor_test <- moran.randtest(mor_plants, listw = attr(mor_plants, "listw"))
mor_test <- apply(mor_plants, 2, moran.mc, listw = attr(mor_plants, "listw"), nsim = 999)
# we select the vectors with > 0.7 statistic (coarse scale) 
# and the last significant one (fine scale)
vecs <- c(1:3, 9)
mems_sel <- mor_plants[, vecs]
lapply(apply(mems_sel, 2, unique), length)

# colrs <- load_colour_palette()
# par(mfrow = c(1, length(vecs)))
# fac <- 500
# for (vec in vecs) {
# 	plot(jitter(x_spat[, "x"], factor = fac),jitter(x_spat[, "y"], factor = fac), 
# 		 pch = 21,
# 		 #cex = abs(mor_nonnull[, vec]) * 1.5,
# 		 cex = (mor_plants[, vec] + abs(min(mor_plants[, vec])) + 0.5) * 1,
# 		 bg = colrs[[1]][as.numeric(x_facs$pop)])
# 		 #bg = round(mor_plants[, vec] + abs(min(mor_plants[, vec]) * 10)))
# }

#&&& FIGURE S2 &&&&&
png(file.path(dirs$figs,
              paste0("MEMs_", paste(vecs, collapse = "_"), ".png")), 
    height = 10, 
    width = 10, 
    bg = "transparent",
    units = "in", 
    res = 300)
	par(mfrow = c(2, 2))
	for (vec in vecs) {
		plot(mems_sel, 
			 type = "l",
			 ylab = "Eigenvalues",
			 xlab = "Host plant")
	}
dev.off()
#&&&&&&&&&&&&&&&&&&&&

### combine Xs
x_and_mems <- as.data.frame(cbind(x_nums_scaled, x_bools, mems_sel))

# check for correlations
var_cors <- cor(x_and_mems, use = "na.or.complete")
max(var_cors[upper.tri(var_cors, diag = FALSE)]) # highest correlation 0.51
which(var_cors == max(var_cors[upper.tri(var_cors, diag = FALSE)]), arr.ind = TRUE)

# final data frames
data_df <- as.data.frame(cbind(y, cbind(x_nums_scaled, x_bools)))
data_host_df <- data_habi_df <- data_df[, c(1:16)]
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
crf1 <- MRFcov(data_df, family = "binomial", n_nodes = 16)
saveRDS(crf1, file = file.path(dirs$fits, "crf1.rds"))

### with only host variables
crf1_host <- MRFcov(data_host_df, family = "binomial", n_nodes = 16)
saveRDS(crf1_host, file = file.path(dirs$fits, "crf1_host.rds"))

### with only habitat variables
crf1_habi <- MRFcov(data_habi_df, family = "binomial", n_nodes = 16)
saveRDS(crf1_habi, file = file.path(dirs$fits, "crf1_habi.rds"))

# 3.2.2 Spatial CRF with Morans' eigenvectors

### only MEMs
crf1_only_mem <- MRFcov(data_only_mems_df, family = "binomial", n_nodes = 16)
saveRDS(crf1_only_mem, file = file.path(dirs$fits, "crf1_only_mem.rds"))

### with environment
crf1_mem <- MRFcov(data_mems_df, family = "binomial", n_nodes = 16)
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

# 3.4.2 calculating Tjur R2s
    # Tjur R2 : for each of the two categories of the dependent variable, 
    # calculate the mean of the predicted probabilities of an event,
    # then, take the difference between those two means.
probs <- list("env_and_mems" = preds_crf1_mem$Probability_predictions,
              "only_mems" = preds_crf1_only_mem$Probability_predictions,
              "only_env" = preds_crf1$Probability_predictions,
              "only_host" = preds_crf1_host$Probability_predictions,
              "only_habi" = preds_crf1_habi$Probability_predictions,
              "mrf" = preds_mrf1$Probability_predictions)
tjuR2s <- list()
for (p in 1:length(probs)) {
    tmp <- NA
    for (i in 1:ncol(y)){                       # for all species separately
        tmp <- c(tmp, mean(probs[[p]][,i][y[,i] == 1]) - mean(probs[[p]][,i][y[,i] == 0]))            
    }
    tjuR2s[[p]] <- matrix(tmp[-1], ncol = 1)
    rownames(tjuR2s[[p]]) <- colnames(y)
}
names(tjuR2s) <- names(probs)
lapply(tjuR2s, mean)

plot(x = 1:length(tjuR2s), y = seq(0, 0.5, length.out = length(tjuR2s)), type = "n", xaxt = "n", xlab = "", ylab = "Tjur R^2")
for (i in 1:length(tjuR2s)) {
    points(x = rep(i, times = length(tjuR2s[[i]])), y = tjuR2s[[i]], pch = 21, bg = "grey")
    points(x = i, y = mean(tjuR2s[[i]]), pch = 21, bg = "red3", cex = 1.5)
}
#mean(probs[[1]][y == 1]) - mean(probs[[1]][y == 0])   # for the whole data

# 3.4.3 Cross-validation

# environment and MEMs
crf_env_mems_cv <- lapply(seq_len(100), 
                          function(x) { cv_MRF_diag(data = data_mems_df, 
                                                    n_nodes = 16, 
                                                    n_folds = 4, 
                                                    n_cores = 1, 
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
                                                     n_nodes = 16, 
                                                     n_folds = 4, 
                                                     n_cores = 1, 
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
                                                   n_nodes = 16, 
                                                   n_folds = 4, 
                                                   n_cores = 1, 
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
                                                     n_nodes = 16, 
                                                     n_folds = 4, 
                                                     n_cores = 1, 
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
                                                     n_nodes = 16, 
                                                     n_folds = 4, 
                                                     n_cores = 1, 
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
                                           n_nodes = 16, 
                                           n_folds = 4, 
                                           n_cores = 1, 
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
str(cv_res)
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
library(MRFcov)

sampleProp <- 1
nBoot <- 100
nodes <- 16

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
bootedMRF <- readRDS(file = file.path(dirs$fits, paste0("bootedMRF", fileBody)))
bootedCRF_w_env_mems <- readRDS(file = file.path(dirs$fits, 
												 paste0("bootedCRF_w_env_mems", fileBody)))
bootedCRF_w_mems_only <- readRDS(file = file.path(dirs$fits, 
												  paste0("bootedCRF_w_mems_only", 
												  		 fileBody)))
bootedCRF_w_env <- readRDS(file = file.path(dirs$fits, 
											paste0("bootedCRF_w_env", fileBody)))
bootedCRF_w_habi <- readRDS(file = file.path(dirs$fits, 
											 paste0("bootedCRF_w_habi", fileBody)))
bootedCRF_w_host <- readRDS(file = file.path(dirs$fits, 
											 paste0("bootedCRF_w_host", fileBody)))

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
    significances[[i]] <- (sign(booted_models[[i]]$direct_coef_upper90) + sign(booted_models[[i]]$direct_coef_lower90))[,1:nodes]
}
names(significances) <- names(booted_models)

# mean associations
associations <- list()
for (i in 1:length(booted_models)) {
    associations[[i]] <- booted_models[[i]]$direct_coef_means[1:nodes, 2:(nodes + 1)]
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

# make the association matrix symmetric
ass_sig_symm <- associations_sig
for (i in 1:length(ass_sig_symm)) {
    for (k in 1:nrow(ass_sig_symm[[i]])) {
        for (j in 1:ncol(ass_sig_symm[[i]])) {
            if (ass_sig_symm[[i]][k,j] != 0 | ass_sig_symm[[i]][j,k] != 0) {    
                if (ass_sig_symm[[i]][k,j] != ass_sig_symm[[i]][j,k]) {
                    print(c(i,
                    		ass_sig_symm[[i]][k,j],
                            ass_sig_symm[[i]][j,k]))
                    ass_sig_symm[[i]][k,j] <- ass_sig_symm[[i]][k,j] + ass_sig_symm[[i]][j,k]
                    ass_sig_symm[[i]][j,k] <- ass_sig_symm[[i]][k,j]
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
	colnames(ass_for_adj[[i]]) <- sub(pattern = "viridae", 
									  replacement = "", 
									  colnames(ass_for_adj[[i]]))
	colnames(ass_for_adj[[i]]) <- sub(pattern = "viroidae", 
									  replacement = "", 
									  colnames(ass_for_adj[[i]]))
	rownames(ass_for_adj[[i]]) <- sub(pattern = "viridae", 
									  replacement = "", 
									  rownames(ass_for_adj[[i]]))
	rownames(ass_for_adj[[i]]) <- sub(pattern = "viroidae", 
									  replacement = "", 
									  rownames(ass_for_adj[[i]]))
}

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
    igraph::E(adj_mats[[i]])$width <- abs(igraph::E(adj_mats[[i]])$weight) * 4
    adj_mats[[i]] <- igraph::delete.vertices(adj_mats[[i]], 
                                             igraph::degree(adj_mats[[i]]) == 0)
}

vert_lab_deg <- list()
for (i in 1:length(adj_mats)) {
	vert_lab_deg[[i]] <- rep(180, times = length(igraph::vertex_attr(adj_mats[[i]])$name))
}
names(vert_lab_deg) <- names(adj_mats)

vert_lab_deg[["MRF"]][which(igraph::vertex_attr(adj_mats[["MRF"]])$name == "Betaflexi")] <- pi/2
vert_lab_deg[["MRF"]][which(igraph::vertex_attr(adj_mats[["MRF"]])$name == "Bromo")] <- pi/2
vert_lab_deg[["MRF"]][which(igraph::vertex_attr(adj_mats[["MRF"]])$name == "Alphaflexi")] <- pi/2
vert_lab_deg[["MRF"]][which(igraph::vertex_attr(adj_mats[["MRF"]])$name == "Poty")] <- pi/2
vert_lab_deg[["MRF"]][which(igraph::vertex_attr(adj_mats[["MRF"]])$name == "Partiti")] <- pi/2
vert_lab_deg[["CRF_host"]][which(igraph::vertex_attr(adj_mats[["CRF_host"]])$name == "Poty")] <- pi/2
vert_lab_deg[["CRF_habi"]][which(igraph::vertex_attr(adj_mats[["CRF_habi"]])$name == "Gemini")] <- 90
vert_lab_deg[["CRF_habi"]][which(igraph::vertex_attr(adj_mats[["CRF_habi"]])$name == "Tombus")] <- pi/2
vert_lab_deg[["CRF_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_mems"]])$name == "Bromo")] <- 90
vert_lab_deg[["CRF_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_mems"]])$name == "Luteo")] <- pi/2
vert_lab_deg[["CRF_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_mems"]])$name == "Gemini")] <- pi/2
vert_lab_deg[["CRF_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_mems"]])$name == "Partiti")] <- pi/2
vert_lab_deg[["CRF_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_mems"]])$name == "Tombus")] <- pi/2
vert_lab_deg[["CRF_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_mems"]])$name == "Alphaflexi")] <- pi/2
vert_lab_deg[["CRF_env"]][which(igraph::vertex_attr(adj_mats[["CRF_env"]])$name == "Tombus")] <- pi/2
vert_lab_deg[["CRF_env_mems"]][which(igraph::vertex_attr(adj_mats[["CRF_env_mems"]])$name == "Tombus")] <- pi/2

#&&& FIGURE 4A-B, F in main text &&&&&

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
        set.seed(73)
        #set.seed(10)
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
names(ass_sig_symm_cleared)
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
    igraph::E(adj_mats_diffs[[i]])$width <- abs(igraph::E(adj_mats_diffs[[i]])$weight) * 4
    adj_mats_diffs[[i]] <- igraph::delete.vertices(adj_mats_diffs[[i]], 
                                                   igraph::degree(adj_mats_diffs[[i]]) == 0)
}

vert_lab_deg <- list()
for (i in 1:length(adj_mats_diffs)) {
	vert_lab_deg[[i]] <- rep(180, 
							 times = length(igraph::vertex_attr(adj_mats_diffs[[i]])$name))
}
names(vert_lab_deg) <- names(adj_mats_diffs)

vert_lab_deg[["habi"]][which(igraph::vertex_attr(adj_mats_diffs[["habi"]])$name == "Alphaflexi")] <- 0
vert_lab_deg[["habi"]][which(igraph::vertex_attr(adj_mats_diffs[["habi"]])$name == "Betaflexi")] <- pi
vert_lab_deg[["habi"]][which(igraph::vertex_attr(adj_mats_diffs[["habi"]])$name == "Partiti")] <- pi/2
vert_lab_deg[["habi"]][which(igraph::vertex_attr(adj_mats_diffs[["habi"]])$name == "Tombus")] <- pi/2
vert_lab_deg[["habi"]][which(igraph::vertex_attr(adj_mats_diffs[["habi"]])$name == "Luteo")] <- pi/2

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
        set.seed(73)
        #set.seed(10)
        plot(adj_mats_diffs[[i]], 
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
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# Direct covariate effects
mod <- booted_models$CRF_env_mems

names(mod)
cov_inds <- c(1, 18:31)
colnames(mod$direct_coef_means)
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
head(sig_direct_env_effects)
write.table(sig_direct_env_effects,
            file = file.path(dirs$fits,"direct_sig_env_effs.csv"),
            quote = FALSE,
            sep = ";",
            dec = ".")            


##########################################################################################
# 110521
