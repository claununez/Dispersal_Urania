
# Project: Dispersal patterns of Urania boisduvalii (Lepidoptera: Uraniidae) 
#          driven by biotic interactions
# Authors: Claudia Nunez-Penichet
#          Jorge Soberón
#          Luis Osorio-Olvera
#-------------------------------------------------------------------------------

# Time series simulations for different toxicity and connectivity scenarios.
# The dispersal process starts from a point in one extreme of the island
# (eastern or western). The period for Omphalea to became toxic vary across
# the connectivity scenarios (10, 15,and 20 neighbors).

# Note: This script writes the results in your working directory

# Data available at https://doi.org/10.6084/m9.figshare.19217592

## Loading packages
library(bam)
library(raster)
library(purrr)
library(magrittr)
library(furrr)

set.seed(111)
rm(list = ls())
source("time_series_function.R") # this function is included in 
                                 # the repository
rm(list = ls())

# Reading the data
sp1 <- raster::raster("urania_model.tif")
sp2 <- raster::raster("omphalea_model.tif")
sparse_mod <- bam::model2sparse(sp1*sp2)

# Points for starting simulations
cc1 <- read.csv("UraniaData.csv")

## Coordinates for starting simulations in the west
occs_oriente <- c(-74.762, 20.327)

## Coordinates for starting simulations in the east
occs_poniente <- c(-84.479, 21.93)
both_occs <- matrix(c(occs_oriente, occs_poniente), ncol = 2, byrow = T)
occsList <- list(occs_oriente = occs_oriente,
                 occs_poniente = occs_poniente,
                 occs_all = both_occs)

# Defining parameters
periods_positiveL <- rep(c(3, 3, 2, 6, 3, 0), 3)
periods_negativeL <- rep(c(3, 2, 3, 3, 6, 2), 3)
ngbsL <- rev(c(10, 15, 20))

# Calculating adjacent matrices
plan(multisession(workers = length(ngbsL) + 2))
adj_matrixL <- ngbsL %>% furrr::future_map(function(x){
  adj_matrix <- bam::adj_mat(sparse_mod, ngbs = x)
}, .progress = TRUE)

plan(sequential)
gc()

# Preparing parameters for time series calculations
nprocesss <- length(adj_matrixL) *
  length(periods_negativeL) * length(occsList)

seqs_matrices <- rep(seq_along(adj_matrixL), each = length(occsList) * 
                       length(periods_negativeL))

seqs_periods <- rep(rep(seq_along(periods_negativeL), each = length(occsList)),
                    times = length(adj_matrixL))

seqs_regions <- rep(seq_along(occsList), times = length(adj_matrixL) * 
                      length(periods_negativeL))

data.frame(seqs_matrices, seqs_regions, seqs_periods)

# Creating directory for results
dir.create("Results")

# Running in parallel to obtain time series
plan(multisession(workers = 3))

options(future.globals.maxSize = 8500 * 1024^2)
df_all_tseries <- 1:nprocesss %>% furrr::future_map_dfr(function(x){
  set.seed(111)
  matrix_index <- seqs_matrices[x]
  adj_matrix <- adj_matrixL[[matrix_index]]
  periods_negative <- periods_negativeL[seqs_periods[x]]
  periods_positive <- periods_positiveL[seqs_periods[x]]
  focal_region <- NULL
  init_coords <- occsList[[seqs_regions[x]]]
  
  ts_r1 <- ts_simulte(sp1 = sp1,
                      sp2 = sp2,
                      adj_matrix = adj_matrix,
                      periods_negative = periods_negative,
                      periods_positive = periods_positive,
                      nsteps = 220,
                      npops = 1,
                      focal_region = focal_region,
                      init_pop = init_coords)
  
  rcoords <- t(sapply(7:ncol(ts_r1), function(y){
    ids <- which(ts_r1[, y] > 0)
    colMeans(ts_r1[ids, c("x", "y")])
  }))
  
  rasname <- names(occsList[seqs_regions[x]])
  
  
  df_summary <- data.frame(time = colnames(ts_r1[, -(1:6)]),
                           region_id = seqs_regions[x],
                           rasname,
                           ts_r1[1,(1:3)], rcoords,
                           mean_occ=colMeans(ts_r1[, -(1:6)]),
                           row.names = NULL)
  fname = paste0("ts_Reg_",df_summary$region_id[1],"_","ngbs_",
                 df_summary$ngbs[1],"_",df_summary$periods_parm[1],".csv")
  fpath = file.path("Results",fname)
  rio::export(df_summary,fpath)
  return(df_summary)
  
}, .progress = TRUE)

plan(sequential)

# Creating a single table
setwd("Results")

# Reading the tables to merge them
ts_1 <- list.files(path = ".", pattern = "^ts_", full.names = T)
ts_11 <- list.files(path = ".", pattern = "^ts_", full.names = F)
name <- gsub(".csv$", "", ts_11)

ts <-list()

for (i in 1:length(name)){
  ts[[i]] <- read.csv(file = ts_1[i], as.is = FALSE)
} 

spl <- do.call(rbind, ts)
write.csv(spl, "Time_series.csv", row.names = F)
