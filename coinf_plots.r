rm(list = ls(all = TRUE)) ; gc()

# define the working directory again
working_dir <- "/Users/annanorb/Documents/UZH/Projects/META/meta17network-pkg"

setwd(working_dir)
library("meta17network")

dirs <- set_dirs(working_dir = working_dir)
#saveRDS(dirs, file = file.path(working_dir, "dirs.rds"))

# DATA PROCESSING
dat <- process_data(dirs = dirs,
                    return_data = TRUE,
                    save_data = TRUE,
                    rmNAs = TRUE)

# co-occurrence profiles
cooccs <- t(dat$Y) %*% as.matrix(dat$Y)
occs <- dat$Y[, order(colSums(dat$Y), decreasing = TRUE)]
colnames(occs) <- sub("sp_", "", colnames(occs))

coinfs <- list()
for (i in 1:length(unique(dat$X$pop))) {

    patch <- unique(dat$X$pop)[i]
    occs_patch <- occs[which(dat$X$pop == patch), ]

    tmp1 <- toString(names(occs_patch[1, ])[occs_patch[1, ] > 0])
    for (j in 2:nrow(occs_patch)) {
        tmp2 <- names(occs_patch[j, ])[occs_patch[j, ] > 0]
        tmp2 <- toString(tmp2)
        tmp1 <- rbind(tmp1, tmp2)
    }
    coinfs[[i]] <- table(tmp1)
}
names(coinfs) <- unique(dat$X$pop)

all_combs_frqs <- table(unlist(lapply(coinfs, names)))
all_combs_frqs <- all_combs_frqs[order(all_combs_frqs, decreasing = TRUE)]

rares <- as.data.frame(all_combs_frqs[which(all_combs_frqs < 2)])
commons <- as.data.frame(all_combs_frqs[-which(all_combs_frqs < 2)])
rares_ord <- as.data.frame(all_combs_frqs[which(all_combs_frqs < 2)][order(apply(rares, 2, nchar)[,1], decreasing = FALSE)])
all_combs_frqs_ord <- rbind(commons, rares_ord)

viruses_split <- strsplit(as.character(all_combs_frqs_ord[,1]), ",")
viruses_nos <- unlist(lapply(viruses_split, length))

combs_cols <- rep(NA, nrow(all_combs_frqs_ord))
combs_cols[which(viruses_nos == 0)] <- "white"
combs_cols[which(viruses_nos == 1)] <- "grey75"
combs_cols[which(viruses_nos >= 2)] <- "grey50"
combs_cols[which(viruses_nos >= 5)] <- "grey25"
combs_cols[which(viruses_nos >= 10)] <- "grey5"

#combs_cols <- c("white" , wesandcols_sel)
#combs_cols <- c(combs_cols, 
#                rep("gray80", 16),
#                rep("gray60", 21),
#                rep("gray40", 13),
#                rep("gray20", 28))

all_coinf_combs_cols <- cbind(as.character(all_combs_frqs_ord[,1]), combs_cols)
rownames(all_coinf_combs_cols) <- all_coinf_combs_cols[,1]
all_coinf_combs_cols <- all_coinf_combs_cols[order(viruses_nos), ]

all_coinf_combs_cols[which(rownames(all_coinf_combs_cols) == "")] <- "No infection"
rownames(all_coinf_combs_cols)[which(rownames(all_coinf_combs_cols) == "")] <- "No infection"

### FIGURE 2A pies #####
# create subdirectory for coinfections by population
dir.create(file.path(dirs$figs, "coinfections_by_pop"))

for (i in 1:length(coinfs)) {
    cinfs <- coinfs[[i]]
    names(cinfs)[which(names(cinfs) == "")] <- "No infection"
    patch_cols <- all_coinf_combs_cols[sort(match(names(cinfs), rownames(all_coinf_combs_cols))), 2]
    patch <- names(coinfs)[i]
    par(mar = rep(0, times = 4))
    png(file.path(dirs$figs,
                 "coinfections_by_pop", 
                 paste0("coinfections_pop_", patch,".png")), 
        height = 3, 
        width = 3, 
        bg = "transparent",
        units = "in", 
        res = 300)
        pie(cinfs[names(patch_cols)],
            col = patch_cols, 
            labels = NA)
    dev.off()
}

########################
png(file = file.path(dirs$figs, 
                     "coinfections_by_pop", 
                     "coinfections_legend.png"),
    bg = "white", 
    width = 3, 
    height = 3,
    units = "in", 
    res = 300)
    par(family = "serif", mar = rep(0, times = 4), bty = "n")
    plot(x = 1:10, 1:10, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    legend("topleft", 
           legend = c("No infections", 
                      "Single infection", 
                      "2-4 infections", 
                      "5-9 infections", 
                      ">10 infections"), 
			pch = 25,
            col = rep("black", 5),
            pt.bg = unique(combs_cols), 
            bty = "n")
dev.off()


### FIGURE 2C-D #####
nest_pops <- readRDS(file = file.path(dirs$fits, "nest_pops.rds"))
Y_nest_ord <- dat$Y[, colnames(nest_pops$comm)]
virus_colrs <- load_colour_palette()
virus_colrs <- virus_colrs[[1]][1:25]
# based on species richness and population nestedness, 
# the two most contrasting populations are 861 and 946
wpop <- "861"
#wpop <- "946"
toPLot <- Y_nest_ord[which(dat$X$pop == wpop), ]
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
    legend(x = 1, y = 10, legend = colnames(Y_nest_ord), fill = virus_colrs)
dev.off()
#####################

#&&& FIGURE 3A &&&&&
y_all <- as.matrix(dat$Y)
y_all_cocs <- t(y_all) %*% y_all

library(corrplot)

cor_mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

p_mat <- cor_mtest(y_all_cocs, method = "pearson")
p_mat1 <- p_mat
p_mat1[which(p_mat > 0.05)] <- 0.01
p_mat1[which(p_mat <= 0.05)] <- 0.5
diag(p_mat1) <- 0.01
p_mat_diag <- matrix(0, nrow = nrow(p_mat), ncol = ncol(p_mat))
colnames(p_mat_diag) <- rownames(p_mat_diag) <- colnames(y_all_cocs)
diag(p_mat_diag) <- 1
             
png(file.path(dirs$figs, "co_occs_all_rev1.png"),
    height = 11, 
    width = 9, 
    bg = "transparent",
    units = "in", 
    res = 300)
    par(family = "serif", mar = c(0,5,5,5))
	# grid and significant correlations
    corrplot(y_all_cocs,
             is.corr = FALSE,
             method = "shade",
             type = "lower",
             order = "FPC",
             tl.col = "black",
             cl.pos = "n",
             addgrid.col = "black",
# 			 col = colorRampPalette(c("#FFFFFF", 
# 									  "#4393C3", 
# 									  "#2166AC"))(100),
			 col = "white",
             p.mat = p_mat1, 
			 sig.level = 0.05,
			 insig = "pch",
			 pch = 15,
			 pch.cex = 3.5,
			 pch.col = "pink2",
             bg = "transparent")
	# add diagonal colours
    corrplot(y_all_cocs,
             is.corr = FALSE,
			 add = TRUE,
             method = "shade",
             type = "lower",
             order = "FPC",
             tl.col = "black",
             cl.pos = "n",
             addgrid.col = "black",
			 col = "transparent",
             p.mat = p_mat_diag, 
			 sig.level = 0.05,
			 insig = "pch",
			 pch = 15,
			 pch.cex = 3.5,
			 pch.col = "grey75",
             bg = "transparent")
	# add numbers
    corrplot(y_all_cocs, 
			 number.digits = 0,
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

# single infections (inside brackets next to species names)
colSums(dat$Y[rowSums(dat$Y) == 1,])
colSums(dat$Y[rowSums(dat$Y) == 1,]) / colSums(dat$Y)

#&&&&&&&&&&&&&&&&&&&


# (coinfections per population)
coinfs <- coinfs_by_pop(Y = dat$Y, 
                        POPs = dat$X$pop, 
                        dirs = dirs, 
                        pop_lifes = list("861" = NA, "946" = NA))
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
length(all_coinfs)
# saveRDS(all_coinfs, file = file.path(dirs$fits, "all_coinfs.rds"))
