
# Project: Dispersal patterns of Urania boisduvalii (Lepidoptera: Uraniidae) 
#          driven by biotic interactions
# Authors: Claudia Nunez-Penichet
#          Jorge Soberón
#          Luis Osorio-Olvera
#-------------------------------------------------------------------------------

# Creates a plot with the number of clusters obtained with the different 
# dispersal distances

# Note: This script writes the figure in your working directory

# Data available at https://doi.org/10.6084/m9.figshare.19217592

## Loading packages
library(plotrix)
library(bam)

## Establishing working directory
setwd("Your working directory")

## Reading the data
data1 <- read.csv("CSD_data.csv")
y <- data1$Clusters
x <- data1$d

# Creating plot
jpeg("CSD_plot.jpg", width = 16.6, height = 14.5, units = "cm", res = 600)

gap.plot(y = y, x = x, gap = c(230, 570), gap.axis = "y", bgcol = "white",
         breakcol = "white", brw = 0.02, xlab = "Dispersal distance (km)",
         ylab = "Number of clusters", ytics = c(0, 50, 100, 150, 200, 600),
         xtics = c(0, 5, 10, 15, 20), pch = 19, type = "b")
box(bty = "l")

axis.break(2, 226.5, style = "gap")
axis.break(2, 226.5, breakcol = "white", style = "gap")
axis.break(2, 228 * (1 + 0.02), breakcol = "black", style = "slash")

dev.off()

# -------------------Creating maps included in figure 1---------------
km_3 <- bam_clusters(ura_omp, ngbs = 3, plot_model = T)
writeRaster(km_10@raster_map, filename = "3_km_clusters.tif", format = "GTiff")

km_15 <- bam_clusters(ura_omp, ngbs = 15, plot_model = T)
writeRaster(km_15@raster_map, filename = "15_km_clusters.tif", format = "GTiff")