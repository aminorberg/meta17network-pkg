
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

max(colSums(dat$Y)/nrow(dat$Y)) # max prevalence 33%
min(colSums(dat$Y)/nrow(dat$Y)) # min prevalence 0.5%

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

# order plants according to their spatial proximity (first plant selected randomly)
distsXY <- as.matrix(dist(dat$X[, c("x", "y")], diag = FALSE))
rownames(distsXY) <- dat$X$sampleID
colnames(distsXY) <- rownames(distsXY)
diag(distsXY) <- NA

Y_dist_ord_ls <- list()
startPlant <- list()
for (n in 1:100) {
    startPlant[[n]] <- sample(1:nrow(dat$Y), 1)
    Y_dist_ord <- dat$Y[startPlant[[n]],]
    distsXY_tmp <- distsXY
    for (i in 1:(nrow(distsXY)-1)) {
        focalPlant <- which(distsXY_tmp[startPlant[[n]], ] == min(distsXY_tmp[startPlant[[n]], ], 
                                                              na.rm = TRUE))
        #print(focalPlant)
        if (length(focalPlant) > 1) {
            focalPlant <- focalPlant[1]
        }
        distsXY_tmp[startPlant[[n]],focalPlant] <- NA
        Y_dist_ord <- rbind(Y_dist_ord, dat$Y[focalPlant,])        
    }
Y_dist_ord_ls[[n]] <- Y_dist_ord
}
Y_dist_ord_ls_mat <- lapply(Y_dist_ord_ls, as.matrix)

insidSumS <- NA
insidSumS2 <- NA
insids <- list()
insids2 <- list()
sp_area <- list()
for (m in 1:length(Y_dist_ord_ls)) {        

    insids[[1]] <- as.matrix(Y_dist_ord_ls_mat[[m]][1,]) %*% t(Y_dist_ord_ls_mat[[m]][1,])
    tmp <- insids[[1]]
    tmp[upper.tri(tmp, diag = TRUE)] <- 0
    insidsSum <- sum(tmp > 0)

    tmp_comm_mat <- 1 * (apply(Y_dist_ord_ls_mat[[m]], 2, cumsum) > 0)
    insids2[[1]] <- as.matrix(tmp_comm_mat[1,]) %*% t(tmp_comm_mat[1,])
    tmp2 <- insids2[[1]]
    tmp2[upper.tri(tmp2, diag = TRUE)] <- 0
    insidsSum2 <- sum(tmp2 > 0)
    
    
    for (n in 2:nrow(Y_dist_ord_ls_mat[[1]])) {
        insids[[n]] <- t(Y_dist_ord_ls_mat[[m]][1:n,]) %*% Y_dist_ord_ls_mat[[m]][1:n,]        
        tmp <- insids[[n]]
        tmp[upper.tri(tmp, diag = TRUE)] <- 0
        tmp <- sum(tmp > 0)
        insidsSum <- c(insidsSum, tmp)
        
        insids2[[n]] <- t(tmp_comm_mat[1:n,]) %*% tmp_comm_mat[1:n,]
        tmp2 <- insids2[[n]]
        tmp2[upper.tri(tmp2, diag = TRUE)] <- 0
        tmp2 <- sum(tmp2 > 0)
        insidsSum2 <- c(insidsSum2, tmp2)
    }
sp_area[[m]] <- apply(Y_dist_ord_ls_mat[[m]], 2, cumsum)
insidSumS <- cbind(insidSumS, insidsSum)
insidSumS2 <- cbind(insidSumS2, insidsSum2)
}
insidSumS <- insidSumS[,-1]
insidSumS2 <- insidSumS2[,-1]

plot(rowMeans(insidSumS2))
lines(rowMeans(insidSumS2))

tmp1 <- t(dat$Y) %*% as.matrix(dat$Y)
tmp1["Avsunviroidae","Alphasatellitidae"] # the only viruses that never co-occur together
# hence the maximum amount of co-occurring pairs is 299

png(file.path(dirs$figs, "coex2_mean_curve.png"),
    height = 4, 
    width = 7, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(3,3,1,1))
    plot(y = insidSumS2[,1], x = 1:nrow(insidSumS2), 
         type = "l",
         lwd = 1.5,
         lty = 3, 
         ylab = "", 
         xlab = "")
#    for (i in 2:ncol(insidSumS)) {
#        lines(y = insidSumS[,i], 
#        x = 1:nrow(insidSumS),
#        lwd = 1.5)
#    }
dev.off()

png(file.path(dirs$figs, "coex_mean_curve_new.png"),
    height = 4, 
    width = 7, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(3,3,1,1))
    plot(y = insidSumS[,1], x = 1:nrow(insidSumS), 
         type = "l",
         lwd = 1.5,
         lty = 1, 
         ylab = "", 
         xlab = "")
#    for (i in 2:ncol(insidSumS)) {
#        lines(y = insidSumS[,i], 
#        x = 1:nrow(insidSumS),
#        lwd = 1.5)
#    }
dev.off()

png(file.path(dirs$figs, "coex_mean_curve.png"),
    height = 4, 
    width = 7, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(3,3,1,1))
    plot(y = insidSumS[,1], x = 1:nrow(insidSumS), 
         type = "n",
         #ylim = c(0, 300), xlim = c(1, 400), 
         #xaxt = "n",
         #yaxt = "n", 
         lwd = 1.5, 
         ylab = "", 
         xlab = "")
#    for (i in 2:ncol(insidSumS)) {
#        lines(y = insidSumS[,i], 
#        x = 1:nrow(insidSumS),
#        lwd = 1.5)
#    }
    lines(y = rowMeans(insidSumS),
          x = 1:nrow(insidSumS),
          #col = "red3",
          lwd = 2)
dev.off()

tmp <- lapply(lapply(sp_area, function(x){(x > 0) * 1}), rowSums)
sp_area <- simplify2array(tmp)

png(file.path(dirs$figs, "sp_area_mean_curve.png"),
    height = 4, 
    width = 7, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(3,3,1,1))
    plot(y = rowMeans(sp_area), x = 1:nrow(sp_area), 
         type = "n",
         #ylim = c(0, 300), xlim = c(1, 400), 
         #xaxt = "n",
         #yaxt = "n",
         lty = 2, 
         lwd = 1.5, 
         ylab = "", 
         xlab = "")
#    for (i in 2:ncol(insidSumS)) {
#        lines(y = rowSums((sp_area[[i]] > 0)), 
#        x = 1:length(rowSums((sp_area[[i]] > 0))),
#        lwd = 1.5,
#        lty = 2)
#    }
    lines(y = rowMeans(sp_area),
          x = 1:nrow(sp_area),
          lwd = 2,
          #col = "red3",
          lty = 2)
dev.off()


# 1.1 Nestedness analysis
# (see e.g. Rynkiewicz et al 2019)

library(vegan)
nest_plants <- nestednodf(dat$Y)
oecosimu(dat$Y, nestedchecker, "r0", statistic = "C.score")
# Function nestedchecker gives the number of checkerboard units, 
# or 2x2 submatrices where both species occur once but on different sites 
# (Stone & Roberts 1990)
 
# nestedness within populations
nests <- list()
pops <- unique(dat$X[, "pop"])
for (i in 1:length(pops)) {
    nests[[i]] <- nestednodf(dat$Y[which(dat$X[, "pop"] == pops[i]),])
}
plot(nests[[20]], names = TRUE, col = "grey50")

# nestedness over populations
aggr_tmp <- NA
for (i in unique(dat$X$pop)) {
    aggr_tmp <- rbind(aggr_tmp, colSums(dat$Y[which(dat$X$pop == i), ]))
}
Y_aggr_pops <- aggr_tmp[-1,]
Y_aggr_pops[which(Y_aggr_pops > 0)] <- 1
rownames(Y_aggr_pops) <- unique(dat$X$pop)

nest_pops <- nestednodf(Y_aggr_pops)
oecosimu(Y_aggr_pops, nestedchecker, "r0", statistic = "C.score")

### FIGURE 2B #####
png(file.path(dirs$figs, "nestedness_bypop.png"),
    height = 6, 
    width = 6, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(1,3,8,1))
    plot(nest_pops, names = TRUE, col = "grey50", bty = "n")
dev.off()
###################

# 1.1 Co-occurrence illustration and analysis
library(EcoSimR)

c_score(t(dat$Y))
species_combo(t(dat$Y[which(dat$X[, "pop"] == pops[10]),]))
species_combo(t(dat$Y))
checker(t(dat$Y))

# which of the populations are the most species poor and species rich?
which(rowSums(Y_aggr_pops) == max(rowSums(Y_aggr_pops)))
which(rowSums(Y_aggr_pops) == min(rowSums(Y_aggr_pops)))

### FIGURE 2C-D #####
library(corrplot)
virus_colrs <- load_colour_palette()
virus_colrs <- virus_colrs[[1]][1:25]
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
#####################

# virus (co)incidences at the contrasting populations
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

# 1.3 Data processing for incidence-matrix, ordinations and (C)MRF mocdelling

## Y
y_all <- as.matrix(dat$Y)
y_all_cocs <- t(y_all) %*% y_all
y_all_cocs[upper.tri(y_all_cocs)] <- NA
diag(y_all_cocs) <- NA
which(y_all_cocs == 0, arr.ind = TRUE)
# only one pair never coexists: Avsunviroidae and Alphasatellitidae
colSums(y_all[rowSums(y_all) == 1,]) # single infections

sp_subset_thr <- round(nrow(dat$Y) * 0.025) # 10
sp_subset <- colnames(dat$Y)[colSums(dat$Y) >= sp_subset_thr]
dat1 <- dat
dat1$Y <- dat1$Y[, sp_subset]
y <- as.matrix(dat1$Y)
y_cocs <- t(y) %*% y
#y_cocs[upper.tri(y_cocs)] <- NA
#diag(y_cocs) <- NA

# incidences in the subsetted data
### FIGURE 3A #####
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
             method = "number",
             type = "lower",
             order = "FPC",
             tl.col = "black",
             cl.pos = "n",
             addgrid.col = "black",
             col = colorRampPalette(c("white","black"))(100))
dev.off()
###################

# coinfections in the subsetted data
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
saveRDS(y_coinfs, file = file.path(dirs$fits, "y_coinfs.rds"))

png(file.path(dirs$figs,
              "y_coinfections.png"), 
    height = 25, 
    width = 20, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(mar = c(70, 3, 1, 1), family = "serif")
    barplot(y_coinfs, las = 2)
dev.off()


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

# 1.3.3 Spatial eigenvector maps
mor_nonnull <- adespatial::dbmem(x_spat[, c("x", "y")], 
                                 store.listw = FALSE,
                                 MEM.autocor = "non-null")
mems1 <- cbind(mor_nonnull$MEM1, 
               mor_nonnull$MEM2, 
               mor_nonnull$MEM3, 
               mor_nonnull$MEM4)
apply(mems1,2,max)

ade4::s.value(x_spat[, c("x", "y")], mor_nonnull[, 1])



png(file.path(dirs$figs,
              "MEMs.png"), 
    height = 10, 
    width = 10, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(mar = c(1, 1, 1, 1), family = "serif")
    plot(mor_nonnull[, 1:4], 
         SpORcoords = x_spat[, c("x", "y")], 
         "pbackground.col" = "transparent")
dev.off()
     
## combine Xs
x_and_mems <- as.data.frame(cbind(x_nums_scaled, x_bools, mems1))

# 1.3.4 Check for correlations
var_cors <- cor(x_and_mems, use = "na.or.complete")
max(var_cors[upper.tri(var_cors, diag = FALSE)]) # highest correlation 0.56


# 1.3.5 Final data frames
data_df <- as.data.frame(cbind(y, cbind(x_nums_scaled, x_bools)))
data_host_df <- data_df[, c(1:17,20,24:28)]
data_habi_df <- data_df[, c(1:16,18:19,21:23)]
data_mems_df <- as.data.frame(cbind(y, x_and_mems))
data_only_mems_df <- as.data.frame(cbind(y, mems1))

# 2 NMDS
library(vegan)
set.seed(7)
y_no0rows <- y[-which(rowSums(y) == 0),]

nmds_k3_jac <- metaMDS(y_no0rows, 
                       distance = "jaccard",
                       k = 3, 
                       try = 100) 
stressplot(nmds_k3_jac, main = "3D")

nmds_k2_jac <- metaMDS(y_no0rows, 
                       distance = "jaccard",
                       k = 2, 
                       try = 100) 
stressplot(nmds_k2_jac, main = "2D")

# k = 3 is better
saveRDS(nmds_k3_jac, file = file.path(dirs$fits, "nmds_k3_jac.rds"))
nmds_k3_jac <- readRDS(file = file.path(dirs$fits, "nmds_k3_jac.rds"))

png(file.path(dirs$figs, "nmds_k3_jac.png"),
    height = 8, 
    width = 8, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(5,5,1,1))
    plot(nmds_k3_jac, 
         #ylim = c(-1.75, 0.25),
         type = "n")
    points(nmds_k3_jac, dis = "sites", col = popcols, pch = 16, cex = 2)
    text(nmds_k3_jac, dis = "sp", cex = 1.5)
dev.off()

# calculate axis loadings for each virus
vec_sp <- envfit(nmds_k3_jac$points, 
                 y_no0rows,
                 perm = 999)
# check r values
vec_sp <- data.frame(r = vec_sp$vectors$r, viruses = names(vec_sp$vectors$r))
vec_sp[which(vec_sp$r >= 0.1),]


# 3 (Conditional) Markov Random Fields, (C)MRF
library(MRFcov)
?MRFcov

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
lapply(tjuR2s, mean)

# plot(x=tjuR2s[[6]],tjuR2s[[1]])

plot(x = 1:length(tjuR2s), y = seq(0, 0.5, length.out = length(tjuR2s)), type = "n", xaxt = "n", xlab = "", ylab = "Tjur R^2")
for (i in 1:length(tjuR2s)) {
    points(x = rep(i, times = length(tjuR2s[[i]])), y = tjuR2s[[i]], pch = 21, bg = "grey")
    points(x = i, y = mean(tjuR2s[[i]]), pch = 21, bg = "red3", cex = 1.5)
}
#mean(probs[[1]][y == 1]) - mean(probs[[1]][y == 0])   # for the whole data

# 3.4.3 Cross-validation
?cv_MRF_diag

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
saveRDS(cv_res, file = file.path(dirs$fits, "cv_res.rds"))
cv_res <- readRDS(file = file.path(dirs$fits, "cv_res.rds"))
lapply(cv_res, apply, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
#pos_pred = the proportion of true positives out of all predicted positives
    #sum(true_pos, na.rm = TRUE) / (sum(true_pos, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
#sensitivity = prop. of true positives out of all real positives
    # i.e. probability of a positive prediction, given that the virus is present
    # sum(true_pos, na.rm = TRUE)/(sum(true_pos, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
#neg_pred = prop. of true negatives out of all predicted negatives
    #sum(true_neg, na.rm = TRUE) / (sum(true_neg, na.rm = TRUE) + sum(false_neg, na.rm = TRUE))
#specificity = prop. of true negatives out of all real negatives 
    # i.e. probability of a negative prediction, given that the virus is absent
    #sum(true_neg, na.rm = TRUE)/(sum(true_neg, na.rm = TRUE) + sum(false_pos, na.rm = TRUE))
#tot_pred = prop. of all true predictions out of the whole data
    #(sum(true_pos, na.rm = TRUE) + sum(true_neg, na.rm = TRUE))/(length(as.matrix(test_data)))

# MRFs predicts more false positives than false negatives
# higher sensitivity than CRFs, but lower specificity

# CRFs have very little difference among each other,
# in general a high proportion of true positives out of all predicted positives
# lowish sensitivity, so a fair amount of false negatives
# high specificity, so very few false positives

MRFcov:::cv_MRF_diag

# 3.5 Fitting final CRFs with bootstrapping (for obtaining parameter uncertainties)
library(MRFcov)
MRFcov:::bootstrap_MRF
?bootstrap_MRF

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

booted_models <- list("MRF" = bootedMRF, 
                      "CRF_env_mems" = bootedCRF_w_env_mems, 
                      "CRF_env" = bootedCRF_w_env, 
                      "CRF_mems" = bootedCRF_w_mems_only,
                      "CRF_habi" = bootedCRF_w_habi, 
                      "CRF_host" = bootedCRF_w_host)
saveRDS(booted_models, 
        file = file.path(dirs$fits, 
                         paste0("booted_models",
                         "_sampleprop", 
                         sampleProp, 
                         "_boots", 
                         nBoot,
                         ".rds")))
booted_models <- readRDS(file = file.path(dirs$fits, 
                                          paste0("booted_models",
                                          "_sampleprop", 
                                          sampleProp, 
                                          "_boots", 
                                          nBoot,
                                          ".rds")))
MRFcov:::plotMRF_hm
tmp <- plotMRF_hm(MRF_mod = booted_models[[2]], 
                 #node_names, 
                 #main, 
                 plot_observed_vals = FALSE, 
                 data = data_mems_df)

head(booted_models$CRF_env_mems$direct_coef_means)
head(tmp$data)
?bootstrap_MRF

# association significances (does the 90% interval overlap with zero)
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

# no changes in the directions of significant associations
assoc_signs[[6]][which(assoc_signs[[1]] > 0, arr.ind = TRUE)]
assoc_signs[[6]][which(assoc_signs[[1]] < 0, arr.ind = TRUE)]

# make the association matrix symmetric
ass_sig_symm <- associations_sig
for (i in 1:length(ass_sig_symm)) {
    for (k in 1:nrow(ass_sig_symm[[i]])) {
        for (j in 1:ncol(ass_sig_symm[[i]])) {
            if (ass_sig_symm[[i]][k,j] != 0 | ass_sig_symm[[i]][j,k] != 0) {    
                if (ass_sig_symm[[i]][k,j] != ass_sig_symm[[i]][j,k]) {
                    print(c(ass_sig_symm[[i]][k,j],
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

# comparison to phylogeny
phyl <- ape::keep.tip(dat$phylogeny, colnames(nonzero_assocs[[1]]))
#plot(igraph::as.igraph(phyl))
phylmat <- cophenetic(phyl)
phyl_comp <- list()
for (i in 1:length(ass_sig_symm)) {
    phyl_comp[[i]] <- vegan::mantel(xdis = dist(ass_sig_symm[[i]]), 
                                    ydis = dist(phylmat))
}
names(phyl_comp) <- names(ass_sig_symm)
# no connection between phylogeny and associations

perm_assoc <- ass_sig_symm[[2]]
perm_assoc[which(assoc_sum != 6, arr.ind = TRUE)] <- 0
vegan::mantel(xdis = dist(perm_assoc), ydis = dist(phylmat))
perm_assoc_cleared <- perm_assoc
perm_assoc_cleared[upper.tri(perm_assoc_cleared)] <- 0

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


# create adjacency matrices for plotting the networks
adj_mats <- list()
for (i in 1:length(ass_sig_symm_cleared)) {
    adj_mats[[i]] <- igraph::graph.adjacency(ass_sig_symm_cleared[[i]], 
                                             weighted = TRUE, 
                                             mode = "undirected")
    # delete edges that represent weak interactions
    adj_mats[[i]] <- igraph::delete.edges(adj_mats[[i]], 
                                          which(abs(igraph::E(adj_mats[[i]])$weight) <= 0.05))
}

adj_mats[[length(adj_mats) + 1]] <- igraph::graph.adjacency(perm_assoc, 
                                                            weighted = TRUE, 
                                                            mode = "undirected")
# delete edges that represent weak interactions
adj_mats[[length(adj_mats)]] <- igraph::delete.edges(adj_mats[[length(adj_mats)]], 
                                                     which(abs(igraph::E(adj_mats[[length(adj_mats)]])$weight) <= 0.05))

names(adj_mats) <- c(names(ass_sig_symm_cleared), "permanent") 
plot(adj_mats[[7]])

cols <- c(grDevices::adjustcolor("blue2", alpha.f = 0.95), 
          grDevices::adjustcolor("red", alpha.f = 0.95))

for (i in 1:length(adj_mats)) {
    igraph::E(adj_mats[[i]])$color <- ifelse(igraph::E(adj_mats[[i]])$weight < 0, 
                                             cols[1], 
                                             cols[2])
    igraph::E(adj_mats[[i]])$width <- abs(igraph::E(adj_mats[[i]])$weight) * 1.5
    adj_mats[[i]] <- igraph::delete.vertices(adj_mats[[i]], 
                                             igraph::degree(adj_mats[[i]]) == 0)
}

### FIGURE 4 in main text ######

for (i in 1:length(adj_mats)) {
    png(file.path(dirs$figs,
                  paste0("igraph_", names(adj_mats)[i], ".png")), 
        height = 10, 
        width = 10, 
        bg = "transparent",
        units = "in", 
        res = 300)
        par(mar = c(1, 1, 1, 10), family = "serif")
        set.seed(73)
        #set.seed(10)
        plot(adj_mats[[i]], 
             edge.curved = 0.4,
             vertex.size = 6,
             vertex.color = "grey80",
             vertex.label.color = "black",
             vertex.label.font = 2,
             vertex.frame.color = grDevices::adjustcolor("black", alpha.f = 0.8), 
             label = TRUE)
    dev.off()
}

#################################

# Compare associations, where can they accounted to
library(NetComp)
?netDiff
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

adj_mats_diffs <- list("mems" = NA, "host" = NA, "habi" = NA, "env" = NA)
adj_mats_diffs$mems <- igraph::graph.adjacency(mems_effect, 
                                               weighted = TRUE, 
                                               mode = "undirected")
adj_mats_diffs$host <- igraph::graph.adjacency(host_effect, 
                                               weighted = TRUE, 
                                               mode = "undirected")
adj_mats_diffs$habi <- igraph::graph.adjacency(habi_effect, 
                                               weighted = TRUE, 
                                               mode = "undirected")
adj_mats_diffs$env <- igraph::graph.adjacency(env_effect, 
                                               weighted = TRUE, 
                                               mode = "undirected")

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
    igraph::E(adj_mats_diffs[[i]])$width <- abs(igraph::E(adj_mats_diffs[[i]])$weight) * 1.5
    adj_mats_diffs[[i]] <- igraph::delete.vertices(adj_mats_diffs[[i]], 
                                                   igraph::degree(adj_mats_diffs[[i]]) == 0)
}

names(adj_mats_diffs)
for (i in 1:length(adj_mats_diffs)) {
    png(file.path(dirs$figs,
                  paste0("igraph_", names(adj_mats_diffs)[i], "_effect.png")), 
        height = 10, 
        width = 10, 
        bg = "transparent",
        units = "in", 
        res = 300)
        par(mar = c(1, 1, 1, 10), family = "serif")
        set.seed(73)
        #set.seed(10)
        plot(adj_mats_diffs[[i]], 
             edge.curved = 0.4,
             vertex.size = 6,
             vertex.color = "grey80",
             vertex.label.color = "black",
             vertex.label.font = 2,
             vertex.frame.color = grDevices::adjustcolor("black", alpha.f = 0.8), 
             label = TRUE)
    dev.off()
}


# direct covariate effects
mod <- bootedCRF_w_env_mems

names(mod)
cov_inds <- c(1, 18:33)
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
            

# key coefficients
library(circleplot)
for (whichvirus in colnames(y)) {
    tmp <- data.frame(sp1 = mod$mean_key_coefs[whichvirus][[1]]$Variable, 
                      sp2 = rep(whichvirus, 
                            nrow(mod$mean_key_coefs[whichvirus][[1]])),
                      coef = mod$mean_key_coefs[whichvirus][[1]]$Rel_importance)
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
        tmp2 <- length(grep(v1, mod$mean_key_coefs[[v2]]$Variable))
        tmp1 <- c(tmp1, tmp2)
    }
    infs <- tmp1[-1]
    print(c(v1, sum(infs > 0)))
}
# Bromoviridae is a key coefficient for 9 viruses, 
# the rest, max 6 (Virga), all at least 1



# 210221 Tähän asti katsottu koodit
# alla olevat plottauskoodit vanhoja ja sekasin...

















# 4 Results & Figures
library(circleplot)

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


###* FIGURE 2 IN MAIN TEXT *###

# maps
make_maps(dirs = dirs, POPs = dat$X$pop)
# (co)infection profile piecharts
plot_pops(dat = dat, dirs = dirs)
##### 150221
##### include pop barplots to plot_pops

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
###* end of figure 2 in main text *###

# 3.2 Co-occurrences
res_cooc <- readRDS(file = file.path(dirs$fits, "res_cooc.rds"))
names(res_cooc)
names(res_cooc$significants)

###* FIGURE 2 IN MAIN TEXT *###

# percentage of significant pairs our of analysed pairs (85%)
res_cooc$cooccur_full_data$co_occurrences / (res_cooc$cooccur_full_data$pot_pairs - res_cooc$cooccur_full_data$omitted)
# percentage of significant pairs our of all possible pairs (31%)
res_cooc$cooccur_full_data$co_occurrences / res_cooc$cooccur_full_data$pot_pairs

plotops <- c("all", "aggr_by_pop")
for (whichplot in plotops) {
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
bootedCRF_w_mems <- readRDS(file = file.path(dirs$fits, 
                                             "bootedCRF_w_mems_sampleprop100_n100.rds"))
bootedMRF <- readRDS(file = file.path(dirs$fits, 
                                      "bootedMRF_sampleprop100_n100.rds"))

names(bootedCRF_w_mems)
colnames(bootedCRF_w_mems$direct_coef_means)
head(bootedCRF_w_mems$direct_coef_means)
names(bootedCRF_w_mems$indirect_coef_mean)
bootedCRF_w_mems$indirect_coef_mean[[14]]

# IGRAPH 091120
nodes <- 16

# CRF
interactions <- bootedCRF_w_mems$direct_coef_means[1:nodes, 2:(nodes + 1)]
interactions_crf <- interactions
#bootedCRF_w_mems$indirect_coef_mean[[1]]

# MRF
interactions <- bootedMRF$direct_coef_means[, 2:(nodes + 1)]
interactions_mrf <- interactions

sum(interactions_crf[upper.tri(interactions_crf, diag = FALSE)] != 0)
sum(interactions_mrf[upper.tri(interactions_mrf, diag = FALSE)] != 0)

all_ints <- matrix(NA, nrow = nrow(interactions_mrf), ncol = ncol(interactions_mrf))
for (i in 1:nrow(interactions_mrf)) {
    for (j in 1:ncol(interactions_mrf)) {
    
        all_ints[i,j] <- abs(sign(interactions_mrf[i,j])) + abs(sign(interactions_crf[i,j]))
    
    }
}
MRFcov:::MRFcov
interactions_mrf[(all_ints[upper.tri(all_ints, diag = FALSE)] != 1)] <- 0

colSums(interactions_mrf != 0)


# 121120 mieti miten plottaat ympäristövaikutukset !

# comparison to phylogeny
phyl <- ape::keep.tip(dat$phylogeny, colnames(interactions))
#plot(igraph::as.igraph(phyl))
phylmat <- cophenetic(phyl)
vegan::mantel(xdis = dist(interactions), ydis = dist(phylmat))
# no connection between phylogeny and interactions


interactions <- interactions_mrf
# Create an adjacency matrix
adj_mat <- igraph::graph.adjacency(interactions, weighted = TRUE, mode = "undirected")

#^^^

# Delete edges that represent weak interactions
adj_mat <- igraph::delete.edges(adj_mat, which(abs(igraph::E(adj_mat)$weight) <= 0.06))

cols <- c(grDevices::adjustcolor("blue2", alpha.f = 0.95), 
          grDevices::adjustcolor("red", alpha.f = 0.95))
igraph::E(adj_mat)$color <- ifelse(igraph::E(adj_mat)$weight < 0, cols[1], cols[2])



# Add colours to show Groups
#cols <- c("#000000", "#009E73", "#e79f00", 
#    "#9ad0f3", "#0072B2", "#D55E00", "#CC79A7")
#vertex_groups <- summary_table_Phylo$Group[which(summary_table_Phylo$OTU %in% 
#    igraph::V(adj.mat)$name)]
#igraph::V(adj.mat)$color <- cols[as.factor(vertex_groups)]
#
#igraph::V(adj.mat)$label <- ""  #this removes labels

# Adjust weights for easier visualisation of interaction strengths
igraph::E(adj_mat)$width <- abs(igraph::E(adj_mat)$weight) * 1.5

# Remove vertices with no edges (no interactions)
adj_mat <- igraph::delete.vertices(adj_mat, igraph::degree(adj_mat) == 0)

quartz()

png(file.path(dirs$figs,
              "igraph_mrf.png"), 
    height = 10, 
    width = 10, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(mar = c(1, 1, 1, 10), family = "serif")
    set.seed(73)
    #set.seed(10)
    plot(adj_mat, 
         edge.curved = 0.4,
         vertex.size = 6,
         vertex.color = "grey80",
         vertex.label.color = "black",
         vertex.label.font = 2,
         vertex.frame.color = grDevices::adjustcolor("black", alpha.f = 0.8), 
         label = TRUE)
dev.off()

#legend("bottomright", bty = "n", legend = levels(as.factor(vertex_groups)), 
#    col = cols, border = NA, pch = 16, cex = 0.9, 
#    inset = c(-0.45, 0))
    
    
    


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

names(bootedCRF_w_mems)
bootedCRF_w_mems$indirect_coef_mean[1]

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
names(bootedCRF_w_mems)
cov_inds <- 18:33
colnames(bootedCRF_w_mems$direct_coef_means)
bootedCRF_w_mems$direct_coef_means[,cov_inds]

signs <- sign(bootedCRF_w_mems$direct_coef_upper90) + sign(bootedCRF_w_mems$direct_coef_lower90)
sigs <- which(signs == -2 | signs == 2, arr.ind = TRUE)

cbind(rownames(bootedCRF_w_mems$direct_coef_means)[sigs[,1]], 
      colnames(bootedCRF_w_mems$direct_coef_means)[sigs[,2]], 
      bootedCRF_w_mems$direct_coef_means[sigs], 
      signs[sigs])

sig_direct_env_effects <- cbind(rownames(bootedCRF_w_mems$direct_coef_means)[sigs[,1]], 
                                colnames(bootedCRF_w_mems$direct_coef_means)[sigs[,2]], 
                                round(bootedCRF_w_mems$direct_coef_means[sigs], 4))

write.table(sig_direct_env_effects,
            file = file.path(dirs$fits,"direct_sig_env_effs.csv"),
            quote = FALSE,
            sep = ";",
            dec = ".")

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
    infs <- tmp1[-1]
    print(c(v1, sum(infs > 0)))
}
# Bromoviridae is a key coefficient for 9 viruses, 
# the rest, max 5, all at least 1


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
