# --------------------------- #
#
# R Script To Make Gene Flow Diagrams For Eunicella verrucosa
#
# Reference: Holland, Jenkins and Stevens 2017
# https://doi.org/10.1038/hdy.2017.14
#
# Tom Jenkins tom.l.jenkins@outlook.com
# 
# --------------------------- #

# Load packages
library(circlize)
library(tidyverse)


# ----------------- #
#
# Prepare data
#
# ----------------- #

# Prepare BayesAss data
baysass = tibble(Region = c("Portugal","France","Ireland","Britain"),
                 Portgual = c(0.9799,0.0085,0.0031,0.0085),
                 France   = c(0.0368,0.6674,0.0024,0.2834),
                 Ireland  = c(0.0331,0.0098,0.9647,0.0123),
                 Britain  = c(0.0046,0.0030,0.0017,0.9907)
                 )
baysass

# Convert to matrix
baysass.mat = as.matrix(baysass[, 2:5])
baysass.mat

# Gene flow matrix
dimnames(baysass.mat) = list(source = baysass$Region, sink = baysass$Region)
baysass.mat

# Prepare Migrate-n data
migrate = tibble(Region = c("Portugal","France","Ireland","Britain"),
                 Portgual = c(0, 21.4, 15.2, 48.2),
                 France   = c(33.4, 0, 15.1, 61.9),
                 Ireland  = c(51.6, 19.3, 0, 164.5),
                 Britain  = c(10, 11, 2.9, 0)
                 )
migrate

# Convert to matrix
migrate.mat = as.matrix(migrate[, 2:5])
migrate.mat

# Gene flow matrix
dimnames(migrate.mat) = list(source = migrate$Region, sink = migrate$Region)
migrate.mat

# Define region colours
cols = c("darkorange2","blue","green3","red2")


# ----------------- #
#
# BayesAss gene flow
#
# ----------------- #

# Reset the circular layout parameters
circos.clear()

# Parameters for the circular layout
circos.par(start.degree = 90, gap.degree = 6)

# Plot chord diagram
chordDiagram(x = baysass.mat, grid.col = cols, grid.border = "black", transparency = 0.25,
             order = baysass$Region, directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.1),
             annotationTrack = "grid", annotationTrackHeight = c(0.1, 0.1),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.15, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)

# Add labels to chord diagram
circos.trackPlotRegion(track.index = 1,
                       bg.border = NA,
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    # Text direction
    theta = circlize(mean(xlim), 1)[1, 1] %% 360
    dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
    circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
                niceFacing = TRUE, cex = 2, font = 1)
  }
)

# Add title
mtext("BayesAss", outer = FALSE, cex = 2, font = 2, line = -1)

# Export figure
# dev.copy2pdf(file = "geneflow_baysass.pdf", height = 10, width = 10)


# ----------------- #
#
# Migrate-n gene flow
#
# ----------------- #

# Reset the circular layout parameters
circos.clear()

# Parameters for the circular layout
circos.par(start.degree = 90, gap.degree = 6)

# Plot chord diagram
chordDiagram(x = migrate.mat, grid.col = cols, grid.border = "black", transparency = 0.25,
             order = baysass$Region, directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.1),
             annotationTrack = "grid", annotationTrackHeight = c(0.1, 0.1),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.15, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)

# Add labels to chord diagram
circos.trackPlotRegion(track.index = 1,
                       bg.border = NA,
  panel.fun = function(x, y) {
   xlim = get.cell.meta.data("xlim")
   sector.index = get.cell.meta.data("sector.index")
   # Text direction
   theta = circlize(mean(xlim), 1)[1, 1] %% 360
   dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
   circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
               niceFacing = TRUE, cex = 2, font = 1)
  }
)

# Add title
mtext("Migrate-n", outer = FALSE, cex = 2, font = 2, line = -1)

# Export figure
# dev.copy2pdf(file = "geneflow_migrate-n.pdf", height = 10, width = 10)


# ----------------- #
#
# Combined gene flow diagrams
#
# ----------------- #

# Plot chord diagrams as a single figure
# png("Eunicella_geneflow.png", width = 15, height = 10, units = "in", res = 600, type = "cairo")
par(mfrow=c(1,2))
circos.clear()
circos.par(start.degree = 90, gap.degree = 6)
chordDiagram(x = baysass.mat, grid.col = cols, grid.border = "black", transparency = 0.25,
             order = baysass$Region, directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.1),
             annotationTrack = "grid", annotationTrackHeight = c(0.1, 0.1),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.15, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)
circos.trackPlotRegion(track.index = 1,
                       bg.border = NA,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         sector.index = get.cell.meta.data("sector.index")
                         # Text direction
                         theta = circlize(mean(xlim), 1)[1, 1] %% 360
                         dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
                         circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
                                     niceFacing = TRUE, cex = 2, font = 1)
                       }
)
mtext("BayesAss", outer = FALSE, cex = 2, font = 2, line = -5)
chordDiagram(x = migrate.mat, grid.col = cols, grid.border = "black", transparency = 0.25,
             order = baysass$Region, directional = T, direction.type = "arrows",
             self.link = 1, preAllocateTracks = list(track.height = 0.1),
             annotationTrack = "grid", annotationTrackHeight = c(0.1, 0.1),
             link.border = "NA", link.sort = T, link.decreasing = T,
             link.arr.length = 0.15, link.arr.lty = 3, link.arr.col = "#252525", 
             link.largest.ontop = F)
circos.trackPlotRegion(track.index = 1,
                       bg.border = NA,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         sector.index = get.cell.meta.data("sector.index")
                         # Text direction
                         theta = circlize(mean(xlim), 1)[1, 1] %% 360
                         dd = ifelse(theta < 180 || theta > 360, "bending.inside", "bending.outside")
                         circos.text(x = mean(xlim), y = 0.8, labels = sector.index, facing = dd,
                                     niceFacing = TRUE, cex = 2, font = 1)
                       }
)
mtext("Migrate-n", outer = FALSE, cex = 2, font = 2, line = -5)
# dev.off()
