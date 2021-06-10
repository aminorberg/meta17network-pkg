rm(list = ls(all = TRUE)) ; gc()

# define the working directory again
working_dir <- "/Users/annanorberg/switchdrive/UZH/Projects/META/meta17network-pkg"

setwd(working_dir)
library("meta17network")

dirs <- set_dirs(working_dir = working_dir)
#install.packages(c("spMaps", "smoothr"))

dat <- readRDS(file.path(dirs$mod_dat, "dat.rds"))

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
            fill = unique(combs_cols), 
            bty = "n")
dev.off()


### FIGURE 2C-D #####
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
