### Script to prepare the temperature data and export a csv containing the 
### starting values of the model parameters for the NLLS fitting.
rm(list=ls())
graphics.off()
load("../Data/GPDDFiltered.RData")
df <- read.csv("../Data/ThermResp.csv")

library(lattice)
require(ggplot2)
library(plyr)

prep <- function(){
######################
### Preparing data ###
######################

unique_id <- data.frame()
gen_id <- function(d){
  return(paste(d, collapse="-:-"))
}

unique_id <- apply(as.matrix(df[, c("Species_standardised", "Reference", "Latitude", "Longitude")]), 1 ,gen_id)
df3 <- data.frame(df$Species_standardised, df$Reference, df$Trait_value, df$Temperature, unique_id, stringsAsFactors = FALSE)
# Removes NAs and datasets with fewer than 5 data points.
df3 <- na.omit(df3)
df4 <<- ddply(df3, .(unique_id), subset, length(unique_id) >= 5)
names(df4) <<- c("Species_standardised", "Reference", "Trait_value", "Temperature", "unique_id")
}
###################
### Sample data ###
###################
# 
# uid <- as.character(df4$unique_id)[100]
# uid_index <- which(as.character(df3$unique_id) == uid)
# 
# samp <- df4[uid_index,]

#################
### Full data ###
#################
plot <- function(){
plot_graphs <- function(d){
  p <- ggplot(d, aes(x = Temperature, y = Trait_value)) + 
    geom_point(size = I(2), shape = I(3)) + 
    xlab("Temperature") + 
    ylab("Trait Value") + 
    labs(title = d$unique_id)
  print(p)
}
pdf("../Results/Plots_by_unique_id.pdf")
d_ply(df4, .(unique_id), plot_graphs)
dev.off()
}
###############################
### Finding E values and B0 ###
###############################
param <- function(){
find_values <- function(f){
  g <- f[with(f, order(Trait_value, Temperature)),]
  max_trait <- max(g$Trait_value)
  max_temp <- g[which(g$Trait_value == max_trait),]
  i <- subset(g, g$Trait_value <= max_trait)
  j <- subset(i, 0< g$Trait_value)
  out <- subset(j, j$Temperature <= max_temp$Temperature)
  k <- 8.617*10^(-5) ## Defining global variable, K, the Boltzmann Constant
  out$Temperature <- -1/(k*(out$Temperature + 273.15)) ## Redefining Temp as 1/kT and in kelvin
  x_10 <- -1/k*293.15 ## out_temp at 10 degrees celsius.
  mod <- lm(log(Trait_value) ~ Temperature, out) ## Creation of the linear model
  B0_value <- min(j$Trait_value)## Starting value for B0 (equal to lowest trait_value)
  E_value <- mod$coefficients[[2]][[1]] ## Activation energy ( = -gradient of 1/kT)
  ED_value <- E_value * 2 ## The Deactivation energy (some multiplier of the activation energy)
  T_pk <- max_temp$Temperature + 273.15 ## Outputs max temp in K
  if(is.na(E_value)){
  }
  else{
    if(E_value < 0){
      ED_value <- 0.0001
  }
    if(ED_value < 0){
      ED_value <- 0.0002
  }}
  data.frame(T_pk, B0_value, E_value, ED_value)
}

####################
### Write to CSV ###
####################

csvdata <- ddply(df4, .(unique_id), find_values)
csvout <- merge(df4, csvdata, by="unique_id")
csvout[is.na(csvout)] <- 0.1
write.csv(csvout,"../Data/ThermResp_startvals.csv")
}
# ### Sample data export
# samp
# y <- ddply(samp, .(unique_id), find_values)
# sampout <- merge(samp, y, by = "unique_id")
# write.csv(sampout,"../Data/samp.csv")
