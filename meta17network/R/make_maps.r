#' Maps
#'
#' Make a map of Åland and an index map
#' @param dirs Directories
#' @return Saved maps as .png
#' @export

make_maps <- function(dirs, POPs) {

    library(spMaps)

    load(file.path(dirs$aland, "ComplexMapAland.Rdata"))
    load(file.path(dirs$aland, "SpatialDataPatchesAndRoads.Rdata"))
    metapops <- read.csv(file = file.path(dirs$aland, "Metapopdata2Anna.csv"))

    # index map
    png(file.path(dirs$figs, 
                  "nordic_index_map.png"), 
            height = 5, 
            width = 5, 
            units = "in", 
            res = 300,
            bg = "transparent")
            par(family = "serif", mar = rep(0, times = 4))
                ref_table <- spMaps::getEuropeReferenceTable()
                subs_eur <- spMaps::getSpMaps(countries = c("FIN", "SWE", "NOR"))
                subs_eur <- raster::crop(subs_eur, raster::extent(5, 35, 56, 72))
                plot(subs_eur)
    dev.off()


    map_aland2_crop <- raster::crop(map_aland2, raster::extent(19.3, 20.9, 59.95, 60.4))
    map_aland2_crop_sm <- smoothr::smooth(map_aland2_crop, method = "ksmooth", smoothness = 20)
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
        for (i in 1:length(unique(POPs))) {
            ptch <- unique(POPs)[i]
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

}
