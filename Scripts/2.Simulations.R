
# Project: Dispersal patterns of Urania boisduvalii (Lepidoptera: Uraniidae) 
#          driven by biotic interactions
# Authors: Claudia Nunez-Penichet
#          Jorge Soberón
#          Luis Osorio-Olvera
#-------------------------------------------------------------------------------

# Script to perform simulations of dispersal of Urania boisduvalii among 
# Omphalea plants

# Note: This script writes the results in your working directory

# Data available at https://doi.org/10.6084/m9.figshare.19217592

## Required packages
devtools::install_github("luismurao/bam")
library(raster)
library(bam)
library(dplyr)

## Setting working directory
setwd("Your working directory")
source("Directory with the function/bam_new.R") # this function is included in 
                                                # the repository
rm(list = ls())

## Reading and plotting models of suitability (1 km resolution)
# Omphalea spp. model (from Nunez-Penichet et al. 2019)
omp_path <- "omphalea_model.tif"
omp <- raster(omp_path)
plot(omp)

# Urania boisduvalii model (from Nunez-Penichet et al. 2019)
ura_path <- "urania_model.tif"
ura <- raster(ura_path)
plot(ura)

## Calculating the product of the suitability models of U. boisduvalii and Omphalea 
#  to find the adjacent matrix
modelsparse <- model2sparse(ura*omp)

## Reading occurrence points of U. boisduvalii
cc1 <- read.csv("UraniaData.csv")

raster::plot(ura*omp)
points(cc1[, 2:3])

## Defining the coordinates to start the simulation
occs_sparse_all <- occs2sparse(modelsparse = modelsparse, occs = cc1[, 2:3])
occs_sparse_east <- occs2sparse(modelsparse = modelsparse, occs = cc1[2, 2:3])
occs_sparse_west <- occs2sparse(modelsparse = modelsparse, occs = cc1[1, 2:3])

## Calculating the adjacent matrix to x neighbors
# 10 neighbors (10 km)
set_M_10 <- adj_mat(modelsparse = modelsparse, ngbs = 10)

# 15 neighbors (15 km)
set_M_15 <- adj_mat(modelsparse = modelsparse, ngbs = 15)

# 20 neighbors (20 km)
set_M_20 <- adj_mat(modelsparse = modelsparse, ngbs = 20)

setwd("Directory for saving simulation results")


##--------------------------Performing simulations------------------------------
#----------------------------Starting West--------------------------------------
nsteps <- 300
### toxic = 0, distance = 10
sdm_ura <- bam::sdm_sim(set_A = modelsparse, set_M = set_M_10, 
                        initial_points = occs_sparse_west, nsteps = nsteps)
                        
new_sim <- bam::sim2Animation(sdm_simul = sdm_ura,
                              which_steps = 1:sdm_ura@sim_steps,
                              ani.res = 200, png_keyword = "1", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nToxic = 0 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "west_d10_tox0.html")


#-------------------------------------------------------------------------------
### toxic = 0, distance = 15
sdm_ura <- bam::sdm_sim(set_A = modelsparse, set_M = set_M_15, 
                        initial_points = occs_sparse_west, nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = sdm_ura,
                              which_steps = 1:sdm_ura@sim_steps,
                              ani.res = 200, png_keyword = "2", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nToxic = 0 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "west_d15_tox0.html")


#-------------------------------------------------------------------------------
### toxic = 0, distance = 20
sdm_ura <- bam::sdm_sim(set_A = modelsparse, set_M = set_M_20, 
                        initial_points = occs_sparse_west, nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = sdm_ura,
                              which_steps = 1:sdm_ura@sim_steps,
                              ani.res = 200, png_keyword = "3", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nToxic = 0 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "west_d20_tox0.html")


#-------------------------------------------------------------------------------
### edible = 3, toxic = 2, distance = 10
periods_toxic  <- 3 #edible
periods_suitable <- 2 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "4", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 2 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "west_d10_edi3_tox2.html")


#---------------------------------------------------
### edible = 3, toxic = 3, distance = 10
periods_toxic  <- 3 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "5", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 3 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "west_d10_edi3_tox3.html")


#---------------------------------------------------
### edible = 2, toxic = 3, distance = 10
periods_toxic  <- 2 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = 300)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "6", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 2 gen", "Toxic = 3 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "west_d10_edi2_tox3.html")


#---------------------------------------------------
nsteps <- 440
### edible = 6, toxic = 3, distance = 10
periods_toxic  <- 6 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "7", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 6 gen", "Toxic = 3 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "west_d10_edi6_tox3.html")


#---------------------------------------------------
### edible = 3, toxic = 6, distance = 10
periods_toxic  <- 3 #edible
periods_suitable <- 6 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "8", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 6 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "west_d10_eat3_tox6.html")


#---------------------------------------------------
#---------------------------------------------------
### edible = 3, toxic = 2, distance = 15
periods_toxic  <- 3 #edible
periods_suitable <- 2 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "9", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 2 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "west_d15_edi3_tox2.html")


#---------------------------------------------------
### edible = 3, toxic = 3, distance = 15
periods_toxic  <- 3 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "10", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 3 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "west_d15_edi3_tox3.html")


#---------------------------------------------------
### edible = 2, toxic = 3, distance = 15
periods_toxic  <- 2 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = 300)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "11", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 2 gen", "Toxic = 3 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "west_d15_edi2_tox3.html")


#---------------------------------------------------
nsteps <- 440
### edible = 6, toxic = 3, distance = 15
periods_toxic  <- 6 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "12", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 6 gen", "Toxic = 3 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "west_d15_edi6_tox3.html")


#---------------------------------------------------
### edible = 3, toxic = 6, distance = 15
periods_toxic  <- 3 #edible
periods_suitable <- 6 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "13", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 6 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "west_d15_eat3_tox6.html")


#---------------------------------------------------
#---------------------------------------------------
### edible = 3, toxic = 2, distance = 20
periods_toxic  <- 3 #edible
periods_suitable <- 2 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "14", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 2 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "west_d20_edi3_tox2.html")


#---------------------------------------------------
### edible = 3, toxic = 3, distance = 20
periods_toxic  <- 3 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "15", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 3 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "west_d20_edi3_tox3.html")


#---------------------------------------------------
### edible = 2, toxic = 3, distance = 20
periods_toxic  <- 2 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = 300)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "16", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 2 gen", "Toxic = 3 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "west_d20_edi2_tox3.html")


#---------------------------------------------------
nsteps <- 440
### edible = 6, toxic = 3, distance = 20
periods_toxic  <- 6 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "17", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 6 gen", "Toxic = 3 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "west_d20_edi6_tox3.html")


#---------------------------------------------------
### edible = 3, toxic = 6, distance = 20
periods_toxic  <- 3 #edible
periods_suitable <- 6 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_west,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "18", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 6 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "west_d20_eat3_tox6.html")


#-------------------------------------------------------------------------------
#----------------------------Starting East--------------------------------------
nsteps <- 300
### toxic = 0, distance = 10
sdm_ura <- bam::sdm_sim(set_A = modelsparse, set_M = set_M_10, 
                        initial_points = occs_sparse_east, nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = sdm_ura,
                              which_steps = 1:sdm_ura@sim_steps,
                              ani.res = 200, png_keyword = "19", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nToxic = 0 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "east_d10_tox0.html")


#-------------------------------------------------------------------------------
### toxic = 0, distance = 15
sdm_ura <- bam::sdm_sim(set_A = modelsparse, set_M = set_M_15, 
                        initial_points = occs_sparse_east, nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = sdm_ura,
                              which_steps = 1:sdm_ura@sim_steps,
                              ani.res = 200, png_keyword = "20", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nToxic = 0 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "east_d15_tox0.html")


#-------------------------------------------------------------------------------
### toxic = 0, distance = 20
sdm_ura <- bam::sdm_sim(set_A = modelsparse, set_M = set_M_20, 
                        initial_points = occs_sparse_east, nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = sdm_ura,
                              which_steps = 1:sdm_ura@sim_steps,
                              ani.res = 200, png_keyword = "21", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nToxic = 0 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "east_d20_tox0.html")


#-------------------------------------------------------------------------------
### edible = 3, toxic = 2, distance = 10
periods_toxic  <- 3 #edible
periods_suitable <- 2 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "22", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 2 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "east_d10_edi3_tox2.html")


#---------------------------------------------------
### edible = 3, toxic = 3, distance = 10
periods_toxic  <- 3 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "23", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 3 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "east_d10_edi3_tox3.html")


#---------------------------------------------------
### edible = 2, toxic = 3, distance = 10
periods_toxic  <- 2 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = 300)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "24", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 2 gen", "Toxic = 3 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "east_d10_edi2_tox3.html")


#---------------------------------------------------
nsteps <- 440
### edible = 6, toxic = 3, distance = 10
periods_toxic  <- 6 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "25", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 6 gen", "Toxic = 3 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "east_d10_edi6_tox3.html")


#---------------------------------------------------
### edible = 3, toxic = 6, distance = 10
periods_toxic  <- 3 #edible
periods_suitable <- 6 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "26", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 6 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "east_d10_eat3_tox6.html")


#---------------------------------------------------
#---------------------------------------------------
### edible = 3, toxic = 2, distance = 15
periods_toxic  <- 3 #edible
periods_suitable <- 2 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "27", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 2 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "east_d15_edi3_tox2.html")


#---------------------------------------------------
### edible = 3, toxic = 3, distance = 15
periods_toxic  <- 3 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "28", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 3 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "east_d15_edi3_tox3.html")


#---------------------------------------------------
### edible = 2, toxic = 3, distance = 15
periods_toxic  <- 2 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = 300)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "29", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 2 gen", "Toxic = 3 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "east_d15_edi2_tox3.html")


#---------------------------------------------------
nsteps <- 440
### edible = 6, toxic = 3, distance = 15
periods_toxic  <- 6 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "30", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 6 gen", "Toxic = 3 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "east_d15_edi6_tox3.html")


#---------------------------------------------------
### edible = 3, toxic = 6, distance = 15
periods_toxic  <- 3 #edible
periods_suitable <- 6 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "31", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 6 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "east_d15_eat3_tox6.html")


#---------------------------------------------------
#---------------------------------------------------
### edible = 3, toxic = 2, distance = 20
periods_toxic  <- 3 #edible
periods_suitable <- 2 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "32", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 2 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "east_d20_edi3_tox2.html")


#---------------------------------------------------
### edible = 3, toxic = 3, distance = 20
periods_toxic  <- 3 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "33", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 3 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "east_d20_edi3_tox3.html")


#---------------------------------------------------
### edible = 2, toxic = 3, distance = 20
periods_toxic  <- 2 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = 300)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "34", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 2 gen", "Toxic = 3 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "east_d20_edi2_tox3.html")


#---------------------------------------------------
nsteps <- 440
### edible = 6, toxic = 3, distance = 20
periods_toxic  <- 6 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "35", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 6 gen", "Toxic = 3 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "east_d20_edi6_tox3.html")


#---------------------------------------------------
### edible = 3, toxic = 6, distance = 20
periods_toxic  <- 3 #edible
periods_suitable <- 6 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_east,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "36", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 6 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "east_d20_eat3_tox6.html")


#-------------------------------------------------------------------------------
#----------------------------Starting Both--------------------------------------
nsteps <- 300
### toxic = 0, distance = 10
sdm_ura <- bam::sdm_sim(set_A = modelsparse, set_M = set_M_10, 
                        initial_points = occs_sparse_all, nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = sdm_ura,
                              which_steps = 1:sdm_ura@sim_steps,
                              ani.res = 200, png_keyword = "37", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nToxic = 0 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "all_d10_tox0.html")


#-------------------------------------------------------------------------------
### toxic = 0, distance = 15
sdm_ura <- bam::sdm_sim(set_A = modelsparse, set_M = set_M_15, 
                        initial_points = occs_sparse_all, nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = sdm_ura,
                              which_steps = 1:sdm_ura@sim_steps,
                              ani.res = 200, png_keyword = "38", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nToxic = 0 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "all_d15_tox0.html")


#-------------------------------------------------------------------------------
### toxic = 0, distance = 20
sdm_ura <- bam::sdm_sim(set_A = modelsparse, set_M = set_M_20, 
                        initial_points = occs_sparse_all, nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = sdm_ura,
                              which_steps = 1:sdm_ura@sim_steps,
                              ani.res = 200, png_keyword = "39", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nToxic = 0 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "all_d20_tox0.html")


#-------------------------------------------------------------------------------
### edible = 3, toxic = 2, distance = 10
periods_toxic  <- 3 #edible
periods_suitable <- 2 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "40", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 2 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "all_d10_edi3_tox2.html")


#---------------------------------------------------
### edible = 3, toxic = 3, distance = 10
periods_toxic  <- 3 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "41", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 3 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "all_d10_edi3_tox3.html")


#---------------------------------------------------
### edible = 2, toxic = 3, distance = 10
periods_toxic  <- 2 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = 300)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "42", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 2 gen", "Toxic = 3 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "all_d10_edi2_tox3.html")


#---------------------------------------------------
nsteps <- 440
### edible = 6, toxic = 3, distance = 10
periods_toxic  <- 6 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "43", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 6 gen", "Toxic = 3 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "all_d10_edi6_tox3.html")


#---------------------------------------------------
### edible = 3, toxic = 6, distance = 10
periods_toxic  <- 3 #edible
periods_suitable <- 6 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_10,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "44", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 6 gen",
                                               "Distance = 10 km"),
                              fmt = "HTML", filename = "all_d10_eat3_tox6.html")


#---------------------------------------------------
#---------------------------------------------------
### edible = 3, toxic = 2, distance = 15
periods_toxic  <- 3 #edible
periods_suitable <- 2 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "45", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 2 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "all_d15_edi3_tox2.html")


#---------------------------------------------------
### edible = 3, toxic = 3, distance = 15
periods_toxic  <- 3 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "46", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 3 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "all_d15_edi3_tox3.html")


#---------------------------------------------------
### edible = 2, toxic = 3, distance = 15
periods_toxic  <- 2 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = 300)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "47", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 2 gen", "Toxic = 3 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "all_d15_edi2_tox3.html")


#---------------------------------------------------
nsteps <- 440
### edible = 6, toxic = 3, distance = 15
periods_toxic  <- 6 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "48", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 6 gen", "Toxic = 3 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "all_d15_edi6_tox3.html")


#---------------------------------------------------
### edible = 3, toxic = 6, distance = 15
periods_toxic  <- 3 #edible
periods_suitable <- 6 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_15,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "49", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 6 gen",
                                               "Distance = 15 km"),
                              fmt = "HTML", filename = "all_d15_eat3_tox6.html")


#---------------------------------------------------
#---------------------------------------------------
### edible = 3, toxic = 2, distance = 20
periods_toxic  <- 3 #edible
periods_suitable <- 2 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "50", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 2 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "all_d20_edi3_tox2.html")


#---------------------------------------------------
### edible = 3, toxic = 3, distance = 20
periods_toxic  <- 3 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "51", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 3 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "all_d20_edi3_tox3.html")


#---------------------------------------------------
### edible = 2, toxic = 3, distance = 20
periods_toxic  <- 2 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = 300)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "52", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 2 gen", "Toxic = 3 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "all_d20_edi2_tox3.html")


#---------------------------------------------------
nsteps <- 440
### edible = 6, toxic = 3, distance = 20
periods_toxic  <- 6 #edible
periods_suitable <- 3 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "53", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 6 gen", "Toxic = 3 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "all_d20_edi6_tox3.html")


#---------------------------------------------------
### edible = 3, toxic = 6, distance = 20
periods_toxic  <- 3 #edible
periods_suitable <- 6 #toxic
progress_bar <- T

ura_sim <- bam_sim1(sp1 = ura, sp2 = omp, set_M = set_M_20,
                    initial_points = occs_sparse_all,
                    periods_toxic  = periods_toxic,
                    periods_suitable = periods_suitable,
                    palatable_matrices = TRUE,
                    nsteps = nsteps)

new_sim <- bam::sim2Animation(sdm_simul = ura_sim, which_steps = 1:ura_sim@sim_steps,
                              ani.res = 200, png_keyword = "54", ani.height = 700,
                              bg_color = "lightgray", suit_color = "#008837",
                              occupied_color = "#7b3294", 
                              extra_legend = c("\nEdible = 3 gen", "Toxic = 6 gen",
                                               "Distance = 20 km"),
                              fmt = "HTML", filename = "all_d20_eat3_tox6.html")


#---------------------------------------------------