#' Process meta17 data
#'
#' Modify meta17 data
#' @param save_data Save the data
#' @param return_data Return the data
#' @param rmNAs Remove rows containing NAs from the merged data frame
#' @export

process_data <- function(dirs,
                         return_data = TRUE,
                         save_data = TRUE,
                         rmNAs = TRUE)
{

    # environmental and host-related variables
    ######################################################################################
    
    # agricultural land area
    filename <- "agri_area.csv"
    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    agri_area <- as.data.frame(read.csv2(file.path(datapath))) #, check.names = FALSE
    colnames(agri_area) <- c("pop", "prop_nonagri_area", "prop_agri_area")
    agri_area <- apply(agri_area, 2, as.character)
    agri_area <- apply(agri_area, 2, as.numeric)

    # weather data
    filename <- "weatherdata.csv"
    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    weather <- as.data.frame(read.csv2(file.path(datapath))) #check.names = FALSE
    
    siteIDs <- as.character(unique(weather[, "site"]))

    # snow
    snow <- apply(apply(weather[which(weather$Type == "snow"), ], 
                        2, 
                        as.character), 
                  2, 
                  as.numeric)
    # snow sum year 2016 months 9-12 and year 2017 months 1-5    
    snow_prev_winter <- snow[which(snow[, "Month"] >= 9 & snow[, "Year"] == "2016" | snow[, "Month"] <= 5 & snow[, "Year"] == "2017"), ]
    snowsum_prev_winter <- sum(snow_prev_winter[snow_prev_winter[, "site"] == unique(snow_prev_winter[, "site"])[1], "value"], 
                   na.rm = TRUE)
    for (p1 in 2:length(unique(snow_prev_winter[, "site"]))) {
        p2 <- unique(snow_prev_winter[, "site"])[p1]
        sn_sm <- sum(snow_prev_winter[snow_prev_winter[, "site"] == p2, "value"], na.rm = TRUE)
        snowsum_prev_winter <- rbind(snowsum_prev_winter, sn_sm)
    }
    snowsum_prev_winter <- data.frame("snowsum" = snowsum_prev_winter[,1], 
                                      "pop" = unique(snow_prev_winter[, "site"]))
    rownames(snowsum_prev_winter) <- unique(snow_prev_winter[, "site"])
    snowsum_prev_winter <- snowsum_prev_winter[siteIDs, ]

    # precipitation
    precip <- apply(apply(weather[which(weather$Type == "prec"), ], 
                          2, 
                            as.character), 
                    2, 
                    as.numeric)
    # precipitation sum for months 4-5 2017 
    precip_spring17 <- precip[which(precip[, "Year"] == "2017" & (precip[, "Month"] == 5 | precip[, "Month"] == 4)), ]
    precsum_spring17 <- sum(precip_spring17[precip_spring17[, "site"] == unique(precip_spring17[, "site"])[1], "value"], 
                            na.rm = TRUE)
    for (p1 in 2:length(unique(precip_spring17[, "site"]))) {
        p2 <- unique(precip_spring17[, "site"])[p1]
        pr_sm <- sum(precip_spring17[precip_spring17[, "site"] == p2, "value"], na.rm = TRUE)
        precsum_spring17 <- rbind(precsum_spring17, pr_sm)
    }
    precsum_spring17 <- data.frame("precsum_spring17" = precsum_spring17[,1], 
                                      "pop" = unique(precip_spring17[, "site"]))
    rownames(precsum_spring17) <- unique(precip_spring17[, "site"])
    precsum_spring17 <- precsum_spring17[siteIDs, ]

    # precipitation sum for months 6-9 2016 
    precip_summer16 <- precip[which(precip[, "Year"] == "2016" & (precip[, "Month"] >= 6 & precip[, "Month"] <= 9)), ]
    precsum_summer16 <- sum(precip_summer16[precip_summer16[, "site"] == unique(precip_summer16[, "site"])[1], "value"], 
                            na.rm = TRUE)
    for (p1 in 2:length(unique(precip_summer16[, "site"]))) {
        p2 <- unique(precip_summer16[, "site"])[p1]
        pr_sm <- sum(precip_summer16[precip_summer16[, "site"] == p2, "value"], na.rm = TRUE)
        precsum_summer16 <- rbind(precsum_summer16, pr_sm)
    }
    precsum_summer16 <- data.frame("precsum_summer16" = precsum_summer16[,1], 
                                   "pop" = unique(precip_summer16[, "site"]))
    rownames(precsum_summer16) <- unique(precip_summer16[, "site"])
    precsum_summer16 <- precsum_summer16[siteIDs, ]

    # temperature
    temp <- apply(apply(weather[which(weather$Type == "temp"), ], 
                          2, 
                            as.character), 
                    2, 
                    as.numeric)

    # sum over frozen days Sept 2016 - Apr 2017
    temp_frozen <- temp[which(temp[, "Month"] >= 9 & temp[, "Year"] == "2016" | temp[, "Month"] <= 4 & temp[, "Year"] == "2017"), ]
    temp_frozen <- temp_frozen[which(temp_frozen[, "value"] < 0), ]

    tempsum_frozen <- sum(temp_frozen[temp_frozen[, "site"] == unique(temp_frozen[, "site"])[1], "value"], 
                   na.rm = TRUE)
    for (p1 in 2:length(unique(temp_frozen[, "site"]))) {
        p2 <- unique(temp_frozen[, "site"])[p1]
        tmp_sm <- sum(temp_frozen[temp_frozen[, "site"] == p2, "value"], na.rm = TRUE)
        tempsum_frozen <- rbind(tempsum_frozen, tmp_sm)
    }
    tempsum_frozen <- data.frame("tempsum_frozen" = tempsum_frozen[,1], 
                                   "pop" = unique(temp_frozen[, "site"]))
    rownames(tempsum_frozen) <- unique(temp_frozen[, "site"])
    tempsum_frozen <- tempsum_frozen[siteIDs, ]

    # sum over temperatures of year 2016 summer days
    temp_summer16 <- temp[which((temp[, "Month"] >= 5 & temp[, "Month"] <= 8) & temp[, "Year"] == "2016"), ]
    tempsum_summer16 <- sum(temp_summer16[temp_summer16[, "site"] == unique(temp_summer16[, "site"])[1], "value"], 
                   na.rm = TRUE)
    for (p1 in 2:length(unique(temp_summer16[, "site"]))) {
        p2 <- unique(temp_summer16[, "site"])[p1]
        tmp_sm <- sum(temp_summer16[temp_summer16[, "site"] == p2, "value"], na.rm = TRUE)
        tempsum_summer16 <- rbind(tempsum_summer16, tmp_sm)
    }
    tempsum_summer16 <- data.frame("tempsum_summer16" = tempsum_summer16[,1], 
                                   "pop" = unique(temp_summer16[, "site"]))
    rownames(tempsum_summer16) <- unique(temp_summer16[, "site"])
    tempsum_summer16 <- tempsum_summer16[siteIDs, ]

    # yearly metric of pathogen exposure to freezing 
        # the number of days between October and April 
        # in which the minimum air temperature was < 0 C and the snow depth was < 5 cm
    
    # snow & temp
    snow2 <- apply(weather[which(weather$Type == "snow"), ], 
                   2, 
                   as.character)[, c("site", "variable", "value", "Year", "Month")]
    temp2 <- apply(weather[which(weather$Type == "temp"), ], 
                   2, 
                   as.character)[, c("site", "variable", "value", "Year", "Month")]

    substrRight <- function(x, n) {
        substr(x, nchar(x) - n + 1, nchar(x))
    }

    snow2[, "variable"] <- substrRight(snow2[, "variable"], 8)    
    temp2[, "variable"] <- substrRight(temp2[, "variable"], 8)    
    snow2 <- as.data.frame(snow2)
    temp2 <- as.data.frame(temp2)
    snow2 <- apply(apply(snow2, 2, as.character), 2, as.numeric)
    temp2 <- apply(apply(temp2, 2, as.character), 2, as.numeric)
    tmp <- temp2[which(temp2[, "value"] < 0 & ((temp2[, "Month"] >= 10 & temp2[, "Year"] == 2016) | (temp2[, "Month"] <= 4 & temp2[, "Year"] == 2017))), ]
    snw <- snow2[which(snow2[, "value"] >= 5 & ((snow2[, "Month"] >= 10 & snow2[, "Year"] == 2016) | (snow2[, "Month"] <= 4 & snow2[, "Year"] == 2017))), ]
    intersect_ids <- intersect(tmp[, "variable"], snw[, "variable"])
    
    id <- intersect_ids[1]
    sitescount <- tmp[which(tmp[, "variable"] == id), ]
    sites <- sitescount[!duplicated(sitescount[, "site"]), ]
    for (i in 2:length(intersect_ids)) {
        id <- intersect_ids[i]
        sitescount <- tmp[which(tmp[, "variable"] == id), ]
        sites_tmp <- sitescount[!duplicated(sitescount[, "site"]), ]
        sites <- rbind(sites_tmp, sites)        
    }
    severe_days_count <- table(sites[, c("site")])    
    severe_days_count <- data.frame("severe_days_count" = severe_days_count)
    colnames(severe_days_count) <- c("pop", "severe_days_count")
    rownames(severe_days_count) <- severe_days_count[, "pop"]
    severe_days_count <- severe_days_count[siteIDs, ]    

    # bind weather variables together
    weathersum <- data.frame("temp_frozen_days_prev_winter" = tempsum_frozen[, "tempsum_frozen"], 
                             "temp_eff_days_summer16" = tempsum_summer16[, "tempsum_summer16"], 
                             "prec_summer16" = precsum_summer16[, "precsum_summer16"], 
                             "prec_spring17" = precsum_spring17[, "precsum_spring17"], 
                             "snow_prev_winter" = snowsum_prev_winter[, "snowsum"],
                             "severe_winter_days" = severe_days_count[, "severe_days_count"], 
                             "pop" = siteIDs)

    # viruses
    ######################################################################################
    
    # viruses
    filename <- "viruses.csv"
    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    viruses <- as.data.frame(read.csv2(file.path(datapath))) #check.names = FALSE
    spat <- viruses[, c("sampleID","pop")]
    viruses <- viruses[, -which(colnames(viruses) == "pop")]
    viruses <- viruses[, -which(colnames(viruses) == "plant")]
    viruses$sampleID <- as.character(viruses$sampleID)
    simplenames <- lapply(strsplit(colnames(viruses)[-1], "[.]"), 
                          grep, 
                          pattern = "vir|dae", 
                          value = TRUE)    
    colnames(viruses)[-1] <- paste0("sp_", simplenames)
    viruses <- viruses[, -c(which(colnames(viruses) == "sp_virus"), 
                            which(colnames(viruses) == "sp_character(0)"),
                            which(colnames(viruses) == "sp_viroids"))]

    # virus hosts
    filename <- "virus_host.csv"
    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    virus_host <- as.data.frame(read.csv(file.path(datapath))) #check.names = FALSE
    virus_host[,1] <- sub(" ", "", as.character(virus_host[,1]))
    rownames(virus_host) <- paste0("sp_", virus_host[,1])
    virus_host <- virus_host[intersect(rownames(virus_host), colnames(viruses)[-1]), ]
    virus_host <- virus_host[which(virus_host[, "host"] == "plant" | virus_host[, "host"] == "plant_fungus"), ]
    viruses <- viruses[, c("sampleID", rownames(virus_host))]

    # virus phylogenies
    filename <- "virus_phyliptree.phy"
    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    virus_phyl_tree <- ape::read.tree(file = file.path(datapath))
    virus_phyl_tree$tip.label <- paste0("sp_", virus_phyl_tree$tip.label)
    diffr <- setdiff(virus_phyl_tree$tip.label, 
                     colnames(viruses)[-1])
    virus_phyl_tree2 <- virus_phyl_tree
    for (i in 1:length(diffr)) {
        virus_phyl_tree2 <- ape::drop.tip(virus_phyl_tree2, 
                                          tip = diffr[i])
    }
    virus_phyl_tree <- virus_phyl_tree2

    # filter with phylogeny
    viruses2 <- viruses[,c("sampleID", virus_phyl_tree$tip.label)]
    rownames(viruses2) <- viruses2$sampleID
    viruses <- viruses2
    Y <- viruses[, -1]

    # surrounding vegetation, population variables
    ##########################################################################################
    
    # plants
    filename <- "plants.csv"
    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    plants <- as.data.frame(read.csv(file.path(datapath))) #check.names = FALSE
    colnames(plants)[colnames(plants) == "Sample.ID"] <- "sampleID"
    colnames(plants)[colnames(plants) == "Patch"] <- "pop"
    colnames(plants)[colnames(plants) == "Plant.ID"] <- "plant"
    plants <- plants[, -which(colnames(plants) == "pop")]
    plants <- plants[, -which(colnames(plants) == "plant")]
    plants$sampleID <- as.character(plants$sampleID)
    plants$noleaves <- as.numeric(as.character(plants$noleaves))
    plants <- plants[, -which(colnames(plants) == "vein")]
    plants$leaflenght <- as.numeric(as.character(plants$leaflenght))
    plants$leafwidth <- as.numeric(as.character(plants$leafwidth))
    plants$plant_size <- plants$noleaves * (pi * plants$leaflenght/2 * plants$leafwidth/2)
    plants$holes <- plants$roundholes + plants$largehole + plants$window + plants$roundhole
    plants$holes <- (plants$holes > 0) * 1
    plants$roundhole <- NULL
    plants$suck_or_bite <- ((plants$bite + plants$suck) > 0) * 1
    plants$virus_symptoms <- ((plants$curly + plants$yellow + plants$necrotic) > 0) * 1
    plants <- plants[, -c(2:13)]

    # vegetation
    filename <- "vegetation.csv"
    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    vegetation <- as.data.frame(read.csv2(file.path(datapath))) #check.names = FALSE
    vegetation[, c(1:2, 5:9)] <- apply(vegetation[, c(1:2, 5:9)], 2, as.numeric)
    vegetable <- table(vegetation$Patch, vegetation$Plant)
    vegediversity <- cbind(rownames(vegetable), 
                           rowSums(vegetable), 
                           vegan::diversity(vegetable, "shannon"),
                           vegan::diversity(vegetable, "simpson"),
                           vegan::diversity(vegetable, "invsimpson"))
    colnames(vegediversity) <- c("pop", "plant_rich", "shannon", "simpson", "inv_simpson")
    vegediversity <- apply(vegediversity, 2, as.numeric)
        
    # patches
    filename <- "patchdata.csv"
    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    patches <- as.data.frame(read.csv2(file.path(datapath))) #check.names = FALSE
    patches <- patches[,c(1,7:8)]
    patches$location1 <- as.numeric(as.character(patches$location1))
    patches$location2 <- as.numeric(as.character(patches$location2))

    pops <- unique(patches[, which(colnames(patches) == "Patch")])
    poplocs <- matrix(NA, ncol = 3, nrow = length(pops))
    colnames(poplocs) <- c("pop", "loc1", "loc2")
    rownames(poplocs) <- pops
    for (i in 1:length(pops)) {
        i2 <- pops[i]
        poplocs[i, "pop"] <- unique(patches[which(patches[,1] == i2), 1])
        poplocs[i, c("loc1", "loc2")] <- cluster:::pam(patches[which(patches[,1] == i2), -1], 1)$medoids
    }
    rownames(poplocs) <- NULL

    # focal plant data
    filename <- "focalplantcoordinates.csv"
    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    plantcoords <- as.data.frame(read.csv2(file.path(datapath))) #check.names = FALSE
    plantcoords <- plantcoords[c("Sample", "x", "y")]
    colnames(plantcoords)[which(colnames(plantcoords) == "Sample")] <- "sampleID"
    plantcoords$sampleID <- as.factor(plantcoords$sampleID)
    plantcoords$pop <- simplify2array(strsplit(as.character(plantcoords$sampleID), 
                                      split = "_"))[1, ]
    plantcoords <- apply(plantcoords, 2, as.character)
    plantcoords <- as.data.frame(plantcoords)
    plantcoords[, 2:4] <- apply(plantcoords[, 2:4], 2, as.numeric)

    # populations
    filename <- "populations.csv"
    datapath <- filename
    if (!file.exists(datapath)) {
        datapath <- file.path(dirs$dat, filename)
    }
    populations <- as.data.frame(read.csv2(file.path(datapath))) #check.names = FALSE
    populations$x_pop <- NA
    populations$y_pop <- NA
    populations <- apply(populations, 2, as.character)
    populations <- apply(populations, 2, as.numeric)
    for (i in 1:nrow(populations)) {
        i2 <- as.character(populations[i, "pop"])
        populations[i, c("x_pop", "y_pop")] <- cluster:::pam(plantcoords[which(plantcoords[, "pop"] == i2), c("x", "y")], 
                                                             1)$medoids
    }

    # merge
    joined_data <- merge(x = viruses, y = plants, by = "sampleID", all = FALSE, incomparables = NA)
    joined_data <- merge(x = joined_data, y = spat, by = "sampleID", all = FALSE, incomparables = NA)
    joined_data <- merge(x = joined_data, y = populations, by = "pop", all = TRUE, incomparables = NA)
    joined_data <- merge(x = joined_data, y = vegediversity, by = "pop", all = TRUE, incomparables = NA)
    joined_data <- merge(x = joined_data, y = plantcoords[, -4], by = "sampleID", all = TRUE, incomparables = NA)
    joined_data <- merge(x = joined_data, y = weathersum, by = "pop", all = TRUE, incomparables = NA)
    joined_data <- merge(x = joined_data, y = agri_area, by = "pop", all = TRUE, incomparables = NA)

    sp_ind <- grep("sp_", colnames(joined_data))

    if (rmNAs) {
        joined_data <- joined_data[-unique(which(is.na(joined_data[, -sp_ind]), 
                                                 arr.ind = TRUE)[,1]), ]
    }

    Y <- joined_data[, sp_ind]
    X <- joined_data[, -sp_ind]
    X$sampleID <- as.numeric(as.factor(X$sampleID))
    Xnum <- apply(X, 2, as.numeric)
    X <- as.data.frame(Xnum)

    if (any(colSums(Y, na.rm = TRUE) == 0)) {
        Y <- Y[, -which(colSums(Y, na.rm = TRUE) == 0)]
    }

    colnames(Y) <- sub("sp_", "", colnames(Y))
    virus_phyl_tree$tip.label <- sub("sp_", "", virus_phyl_tree$tip.label)
    
    dat <- list(Y = Y,
                X = X, 
                phylogeny = virus_phyl_tree, 
                population_locations = poplocs, 
                vegetation = vegetation)

    if (save_data) {
        saveRDS(dat, file = file.path(dirs$mod_dat, "dat.rds"))
    }
    if (return_data) {
        return(dat)
    }
}
