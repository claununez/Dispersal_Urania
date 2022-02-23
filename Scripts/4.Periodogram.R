
# Project: Dispersal patterns of Urania boisduvalii (Lepidoptera: Uraniidae) 
#          driven by biotic interactions
# Authors: Claudia Nunez-Penichet
#          Jorge Soberón
#          Luis Osorio-Olvera
#-------------------------------------------------------------------------------

# Script for Fourier analysis and time series

# Note: This script writes the results in your working directory

# Data available at https://doi.org/10.6084/m9.figshare.19217592

## Loading packages
library(pracma)

## Setting working directory
setwd("our working directory")

## Reading the data
raw_data <- read.csv("Time_series.csv")

## Defining parameters
reg <- c(1, 2, 3) # 1:simulation starting in the east, 2: starting in the west, 
                  # and 3: starting in both
relations <- c("negative_3_positive_3", "negative_2_positive_0", 
               "negative_2_positive_3", "negative_3_positive_2", 
               "negative_6_positive_3", "negative_3_positive_6")

dis_10 <- 10 #dispersal distances
dis_15 <- 15
dis_20 <- 20
ex <- c(180, 140, 50, 180, 200, 75)


for (i in reg) {
  st <- ifelse(i == 1, "East", ifelse(i == 2, "West", "all"))
  fig_name <- paste0("Your folder/", st, "_", relations, ".pdf")
  
  for (j in 1:length(relations)) {
    ## Extracting the relationship of edible/toxic of interest
    # 10 km
    suit_txc_10 <- which(raw_data$region_id == i & raw_data$ngbs == dis_10 & 
                           raw_data$periods_parm == relations[j])
    # 15 km
    suit_txc_15 <- which(raw_data$region_id == i & raw_data$ngbs == dis_15 & 
                           raw_data$periods_parm == relations[j])
    # West, 20 km
    suit_txc_20 <- which(raw_data$region_id == i & raw_data$ngbs == dis_20 & 
                           raw_data$periods_parm == relations[j])
    
    ## Extracting columns of interest
    # West, 10 km
    data_10 <- na.omit(as.matrix(raw_data[suit_txc_10, 7:9]))
    # West, 15 km
    data_15 <- na.omit(as.matrix(raw_data[suit_txc_15, 7:9]))
    # West, 20 km
    data_20 <- na.omit(as.matrix(raw_data[suit_txc_20, 7:9]))
    
    #-------------------Spectral analysis------------------
    #
    # Source: https://ms.mcmaster.ca/~bolker/eeid/2010/Ecology/Spectral.pdf
    #
    # Calculates the frequency axis in terms of cycles per sampling interval;

    # frequency = number of cycles per unit time. 
    # period = length of time for one full cycle = 1/frequency
    
    # West, 10 km
    serie_10 <- ts(data_10[, 1][ex[j]:length(data_10[, 1])], 
                   frequency = 1) #deleting transitory
    x.spec_10 <- spectrum(serie_10, span = 5, plot = FALSE)
    spx_10 <- x.spec_10$freq
    spy_10 <- 2 * x.spec_10$spec
    
    # West, 15 km
    serie_15 <- ts(data_15[, 1][ex[j]:length(data_15[, 1])], frequency = 1) 
    x.spec_15 <- spectrum(serie_15, span = 5, plot = FALSE)
    spx_15 <- x.spec_15$freq
    spy_15 <- 2 * x.spec_15$spec
    
    # West, 20 km
    serie_20 <- ts(data_20[, 1][ex[j]:length(data_20[, 1])], frequency = 1) 
    x.spec_20 <- spectrum(serie_20, span = 5, plot = FALSE)
    spx_20 <- x.spec_20$freq
    spy_20 <- 2 * x.spec_20$spec
    
    
    #Finding the peak
    # West, 10 km
    peak_10 <- which(spy_10 == max(spy_10))
    t_peak_10 <- 1/spx_10[peak_10] #time units at which there is a maximum peak
    
    # West, 15 km
    peak_15 <- which(spy_15 == max(spy_15))
    t_peak_15 <- 1/spx_15[peak_15] 
    
    # West, 20 km
    peak_20 <- which(spy_20 == max(spy_20))
    t_peak_20 <- 1/spx_20[peak_20] 
    
    ##-----------------------Plotting-----------------------------------------------
    pdf(fig_name[j], width = 7, height = 4)
    layout(matrix(c(1:12), 3, byrow = T), heights = c(1, 10, 10), 
           widths = c(1, 10, 10, 10))
    par(cex = 0.6)
    
    # Titles
    par(mar = rep(0, 4))
    plot.new() 
    
    plot.new()
    text(0.55, 0.5, "10 km dispersal distance", cex = 1.3)
    
    plot.new()
    text(0.55, 0.5, "15 km dispersal distance", cex = 1.3)
    
    plot.new()
    text(0.55, 0.5, "20 km dispersal distance", cex = 1.3)
    
    plot.new()
    text(0.5, 0.6, "Longitude", cex = 1, srt = 90)
    
    # 10 km
    par(mar = c(4, 2.3, 0.5, 0.5))
    plot(data_10[, 1][ex[j]:length(data_10[, 1])], type = "l", xlab = "", 
         ylab = "", bty = "l")
    title(xlab = "Time", line = 2.3)
    
    # 15 km
    plot(data_15[, 1][ex[j]:length(data_15[, 1])], type = "l", xlab = "", 
         ylab = "", bty = "l")
    title(xlab = "Time", line = 2.3)
    
    # 20 km
    plot(data_20[, 1][ex[j]:length(data_20[, 1])], type = "l", xlab = "", 
         ylab = "", bty = "l")
    title(xlab = "Time", line = 2.3)
    
    par(mar = rep(0, 4))
    plot.new()
    text(0.5, 0.6, "Spectral density", cex = 1, srt = 90)
    
    # 10 km
    par(mar = c(4, 2.3, 0.5, 0.5))
    plot(spy_10 ~ spx_10, xlab = "", ylab = "", type = "l", bty = "l")
    title(xlab = "Cycles per time unit", line = 2.3)
    
    # 15 km
    plot(spy_15 ~ spx_15, xlab = "", ylab = "", type = "l", bty = "l")
    title(xlab = "Cycles per time unit", line = 2.3)
    
    # 20 km
    plot(spy_20 ~ spx_20, xlab = "", ylab = "", type = "l", bty = "l")
    title(xlab = "Cycles per time unit", line = 2.3)
    dev.off()    
  }
}

