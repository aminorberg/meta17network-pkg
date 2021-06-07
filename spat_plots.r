rm(list = ls(all = TRUE)) ; gc()

# define the working directory again
working_dir <- "/Users/annanorberg/switchdrive/UZH/Projects/META/meta17network-pkg"

setwd(working_dir)
library("meta17network")

dirs <- set_dirs(working_dir = working_dir)
#saveRDS(dirs, file = file.path(working_dir, "dirs.rds"))
#install.packages(c("spMaps", "smoothr", "wesanderson"))

dat <- process_data(dirs = dirs,
                    return_data = TRUE,
                    save_data = TRUE,
                    rmNAs = TRUE)

library(spMaps)
library(wesanderson)
library(smoothr)

# maps
load(file.path(dirs$aland, "ComplexMapAland.Rdata"))
load(file.path(dirs$aland, "SpatialDataPatchesAndRoads.Rdata"))
metapops <- read.csv(file = file.path(dirs$aland, "Metapopdata2Anna.csv"))

# Aland centroid location
patch_coords <- as.data.frame(patches[,c(1,3)])
rownames(patch_coords) <- as.character(patch_coords$patch)
apply(patch_coords[as.character(unique(dat$X$pop)), 3:4], 2, mean)

## index map
png(file.path(dirs$figs, 
              "nordic_index_map.png"), 
        height = 5, 
        width = 5, 
        units = "in", 
        res = 300,
        bg = "transparent")
    par(family = "serif", mar = rep(0, times = 4))
    ref_table <- getEuropeReferenceTable()
    subs_eur <- getSpMaps(countries = c("FIN", "SWE", "NOR"))
    subs_eur <- crop(subs_eur, extent(5, 35, 56, 72))
    plot(subs_eur)
dev.off()

map_aland2_crop <- crop(map_aland2, extent(19.3, 20.9, 59.95, 60.4))
map_aland2_crop_sm <- smooth(map_aland2_crop, method = "ksmooth", smoothness = 20)
map_aland <- map_aland2_crop_sm

png(file.path(dirs$figs, 
              paste0("populations_map.png")), 
    height = 8, 
    width = 10, 
    units = "in", 
    res = 300)
    par(family = "serif", mar = rep(0, times = 4))
    plot(patches, col = "white")
    plot(map_aland, add = TRUE, col = "darkolivegreen", border = "grey75")
    for (i in 1:length(unique(dat$X$pop))) {
        ptch <- unique(dat$X$pop)[i]
        points(patches@coords[patches$patch == ptch, 1], 
               patches@coords[patches$patch == ptch, 2], 
               cex = 2, 
               pch = 19, 
               col = "darkolivegreen")
        text(patches@coords[patches$patch == ptch, 1], 
             patches@coords[patches$patch == ptch, 2], 
             labels = as.character(ptch), 
             col = "black")
    }
dev.off()

# populations
wesandcols_sel <- c(wes_palette("GrandBudapest1"),
                    wes_palette("GrandBudapest2"),
                    wes_palette("Zissou1"),
                    wes_palette("Rushmore"),
                    wes_palette("Moonrise1"),
                    wes_palette("Darjeeling2"),
                    wes_palette("IsleofDogs2"),
                    wes_palette("Royal1"),
                    wes_palette("Royal2"),
                    wes_palette("BottleRocket2"),
                    wes_palette("Chevalier1"),
                    wes_palette("IsleofDogs1"),
                    "#800080", "#32cd32", "#ff00ff", "#ffff94", "#9ABCA7", "#4D5382", "#4B7740" 
                    )
rem_cols <- c(1, 7, 11:15, 19, 21:24, 26:29, 31:33, 35, 38, 40:41, 44, 46, 48:50, 54:56)
wesandcols_sel <- wesandcols_sel[-rem_cols]
#par(mfrow=c(1,1))
#varitestikuva(paletti = wesandcols_sel)

wesandcols_sel2 <- wesandcols_sel[1:ncol(dat$Y)]
for (i in 1:length(unique(dat$X$pop))) {
    patch <- unique(dat$X$pop)[i]
    par(mar = rep(0, times = 4))
    png(file.path(dirs$figs,
                  "infections_by_pop", 
                  paste0("infections_pop_", patch,".png")), 
        height = 3, 
        width = 3, 
        bg = "transparent",
        units = "in", 
        res = 300)
        pie(colSums(dat$Y[which(dat$X$pop == patch), ]),
            col = wesandcols_sel2, 
            labels = NA)
    dev.off()
}

leg_cols <- cbind(colnames(dat$Y), wesandcols_sel2)
leg_cols_ordrd <- leg_cols[order(colSums(dat$Y), decreasing = TRUE), ]

png(file = file.path(dirs$figs, 
                     "infections_by_pop", 
                     "infections_legend.png"),
    bg = "transparent", 
    width = 2, 
    height = 6,
    units = "in", 
    res = 300)
    par(family = "serif", mar = rep(0, times = 4))
    plot(x = 1:10, 1:10, type = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
    legend("topleft", legend = leg_cols_ordrd[, 1], fill = leg_cols_ordrd[, 2], bty = 'n')
dev.off()

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
