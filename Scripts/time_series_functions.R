#' Project: Migratory patterns of Urania boisduvalii (Lepidoptera: Uraniidae)
#' as a function of biotic interactions.
#' Paper authors: Claudia Nunez-Penichet
#'                Jorge Soberon
#                 Luis Osorio-Olvera
#' Code author: Luis Osorio-Olvera
#' Function to estimate the time series of pixel occupation from
#' a simulation of the BAM package
#' Date: Feb-5-2021
#' @param sp1 Distribution model for species 1.
#' @param sp2 Distribution model for species 2. The raster must have the same
#' extent and resolution as the raster of sp1.
#' @param adj_matrix Adjacency matrix of an area of interest. see
#' \code{\link[bam]{adj_mat}}.
#' @param periods_negative Periods that sps2 takes to develop defense
#' mechanisms (i.e. toxic). See \code{\link[bam]{bam_sim}} for details.
#' @param periods_positive This is the time that sp2 takes to become non-toxic
#' @param nsteps Number of time steps to run the simulation.
#' @param focal_region A raster of a specific region from which estimate time series.
#' @param init_pop Coordinates of initial populations. If is null the function will take
#' a random coordinate within the focal_raster or the sp1 raster.
library(bam)
library(raster)
library(purrr)
library(magrittr)
ts_simulte <- function(sp1,sp2,adj_matrix,
                       periods_negative = 15,
                       periods_positive = 7,
                       nsteps = 220,
                       npops = 1,
                       focal_region=NULL,
                       init_pop =NULL){
  sp_inter <- sp1*sp2
  sparse_mod <- bam::model2sparse(sp_inter)

  df_coords <-  data.frame(cellIDs=sparse_mod@cellIDs,
                           coordinates(sparse_mod@bin_model)[sparse_mod@cellIDs,])

  if(!is.null(focal_region)){
    all_pchs <- focal_region[]
    nona <- which(!is.na(all_pchs))
    coords <- raster::xyFromCell(focal_region,
                                 which(all_pchs > 0))
  } else{
    coords <- df_coords
    nona <- sparse_mod@cellIDs
  }
  if(is.null(init_pop)){
    coorid <- sample(nrow(coords),1)
    init_pop <- coords[coorid,]
  }

  periods_parm <- paste0("negative_",
                         periods_negative,
                         "_positive_",periods_positive)
  source("paper/Scripts/functions_and_scripts/bam_newAlgo.R")

  init_pop_sparse  <- bam::occs2sparse(modelsparse = sparse_mod,
                                       occs = init_pop)
  sdm_s <- bam_sim1(sp1 = sp1,
                        sp2 = sp2,
                        set_M = adj_matrix,
                        initial_points = init_pop_sparse,
                        periods_toxic = periods_negative,
                        periods_suitable  = periods_positive,
                        nsteps = nsteps)

  for(i in 1:sdm_s@sim_steps){
    df_coords[[paste0("t_",i)]] <- as.vector(sdm_s@sdm_sim[[i]])
  }

  df_serie <- dplyr::inner_join(df_coords,as.data.frame(coords))
  df_serie <- data.frame(region_size =length(nona),
                         ngbs= adj_matrix@ngbs,
                         periods_parm = periods_parm,
                         df_serie)

  return(df_serie)

}

#' Code to collapse a segment of y-axis in a ggplot object
#' Author: The function is taken from Stibu at stackoverflow
#' https://stackoverflow.com/posts/35514861/revisions
#'
#'
squish_trans <- function(from, to, factor) {

  trans <- function(x) {

    if (any(is.na(x))) return(x)

    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to

    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)

    return(x)
  }

  inv <- function(x) {

    if (any(is.na(x))) return(x)

    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor

    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))

    return(x)
  }

  # return the transformation
  return(trans_new("squished", trans, inv))
}
