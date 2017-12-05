#######################################################
# Basic Schoolfield Fitting Code
#######################################################

# set a working directory
setwd("~/Documents/Thermal_Fitting_Repository/Code/")

library(minpack.lm) # for NLLS
library(ggplot2) # for plotting

rm(list = ls())  # clear objects

# assign Boltzmann constant (units of eV * K^-1) as a global:
assign("k", 8.617 * 10^-5, envir = .GlobalEnv)  

#### Functions #####

SchoolFTpk <- function(B0, E, E_D, T_pk, T_ref, temp)
{ # Sharpe-Schoolfield model with explicit T_pk parameter
  
  # PARAMETERS/INPUTS (all temperatures in Kelvin) -
  # temp   : temperature values to evaluate function at (single, scalar or vector of values)
  # B0     : Normalisation constant (log transformed)
  # E      : Activation energy (> 0)
  # E_D    : High temperature de-activation energy (> 0) 
  # T_ref  : Standardization (reference) temperature; set to 0 if not wanted   
  # T_pk   : Temperature at which trait reaches peak value
  
  return(B0 + log(exp((-E/k) * ((1/temp)-(1/T_ref)))/(1 + (E/(E_D - E)) * exp(E_D/k * (1/T_pk - 1/temp)))))
}

#### Estimate STARTING VALUES for the nls

GetE <- function(tmp, rate, T_p, k)
{
  # Estimate starting value for E, taking linear regression using the rise part
  # of the curve only.
  # ~~~ Parameters ~~~
  # tmp  : temperature data (in K).
  # rate : rate data corresponding to temperature above.
  # T_p  : temperature at which rate peaks, used as a cutoff point.
  # k    : Boltzmann constant.
  
  tmp.w <- which(tmp <= T_p)
  if (length(tmp.w) > 1)
  {
    m <- lm(log(rate[tmp.w]) ~ I(1 / (k * (tmp[tmp.w]))))
    return(abs(summary(m)$coefficients[2, 1]))
  } else
  {
    return(0.6)
  }
}

GetB0 <- function(tmp, rate)
{
  # Estimate starting value for the normalising constant.
  # ~~~ Parameters ~~~
  # tmp   : temperature data (in K).
  # rate  : rate data corresponding to temperature above.
  # T_ref : estimate normalising constant at this temperature (in K).
  
  if (min(tmp,na.rm=TRUE) > T_ref)
  {
    return(log(min(rate[1],na.rm=TRUE)))
  } else
  {
    return(log(max(rate[which(tmp <= T_ref)],na.rm=TRUE)))
  }
}


GetTpk <- function(tmp, rate)
{
  # Temperature at which the rate is maximised (estimate of T.peak).
  # ~~~ Parameters ~~~
  # tmp  : Temperature data (in K).
  # rate : Rate data corresponding to temperature above.
  
  return(max(tmp[which.max(rate)]))
}


##-------------------------------------------------------------------------------------##

# load in the data

data <- read.csv("../Data/database.csv")

data$K <- data$ConTemp+273.15 # better make a temperature column in kelvin

# set a reference temperature (note this is in Kelvin, set here at 0C)
T_ref <- 273.15

# pull out the IDs
IDs <- as.character(unique(data$OriginalID))

# CHANGE THIS to set an alternative output directory for the fit graphs
outdir <- "../Results/fits_R/"

# Initialize empty vectors to store the parameter estimates
# that will be obtained.

species <- c()
trait <- c()
B0_sch <- c()
E_sch <- c()
E_D_sch <- c()	
T_pk_est_sch <- c()
T_pk_sch <- c()
P_pk_sch <- c()

# now loop through the IDs and fit...

for(i in 1:length(IDs)){
  # subset the data
  subs <- data[data$OriginalID == IDs[i],]
  
  species <- c(species, as.character(subs$Consumer[1]))
  trait <- c(trait, as.character(subs$StandardisedTraitName[1]))
  
  # generate starting values for the model
  T_pk_st  <- GetTpk(tmp=subs$K, rate=subs$StandardisedTraitValue)
  E_st    <- GetE(tmp=subs$K, rate=subs$StandardisedTraitValue, T_p=T_pk_st, k = k)
  B_st <- GetB0(tmp=subs$K, rate=subs$StandardisedTraitValue)
  
  # try to fit Schoolfield
  try(schoolfield_nls <- nlsLM(
    log(StandardisedTraitValue) ~ SchoolFTpk(B0, E, E_D, T_pk, T_ref, temp = K), data= subs, 
    start=list(B0 = B_st, E = E_st, E_D = 4*E_st, T_pk=T_pk_st)))
  
  
  # If fitting worked ...
  if(!is.na(schoolfield_nls[1])){ 
    
    # Collect the parameter estimates...
    B0_sch <- c(B0_sch, exp(coef(schoolfield_nls)["B0"]))
    E_sch <- c(E_sch, coef(schoolfield_nls)["E"])
    E_D_sch <- c(E_D_sch, coef(schoolfield_nls)["E_D"])
    T_pk_est_sch <- c(T_pk_est_sch, coef(schoolfield_nls)["T_pk"])
  
  
  # Calculate the peak of the curve and its 
  # corresponding temperature value.
  curr_prediction <- predict(schoolfield_nls)
  for (j in 1:length(curr_prediction))
  {
    
    # If we found the maximum performance, exit the loop.
    if (curr_prediction[j] == max(curr_prediction))
    {
      break
    }
  }
  
  T_pk_sch <- c(T_pk_sch, subs$K[j])
  P_pk_sch <- c(P_pk_sch, exp(curr_prediction[j]))
  
  
  ##############################
  # Plotting Schoolfield's fit #
  ##############################
  
  # Create a name for the output file using:
  #	- the original id number
  #   - the species name
  #   - the model
  output_name <- paste(
    i, 
    subs$Consumer[1], 
    'Schoolfield',
    sep = "_"
  )
  
  # Remove any characters that won't look good in a file name,
  # using a regular expression.
  output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
  
  # Convert spaces to underscores.
  output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
  
  # Generate predictions from the model fit...
  tmp_temps <- seq(min(
    floor(subs$K)), 
    ceiling(max(subs$K)
    ), length = 200)
  
  tmp_model <- exp(SchoolFTpk(
    coef(schoolfield_nls)["B0"],
    coef(schoolfield_nls)["E"],
    coef(schoolfield_nls)["E_D"],
    coef(schoolfield_nls)["T_pk"],
    T_ref,
    tmp_temps
  ))
  
  ModelToPlot <- data.frame(
    Temperature = tmp_temps - 273.15, 
    TraitValue = tmp_model
  )
  
  # Prepare the data points of the original values.
  DataToPlot <- data.frame(
    Temperature = subs$K - 273.15, 
    TraitValue = subs$StandardisedTraitValue
  )
  DataToPlot <- na.omit(DataToPlot)
  
  # Plot!
  p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                                    y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                             alpha = 0.7, pch = 21) + 
    geom_line(data = ModelToPlot, 
              aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
              lwd = 1.3) + 
    ggtitle(paste(subs$Consumer[1])) +
    xlab(expression(paste("Temperature (", degree, C, ")"))) + 
    ylab(subs$StandardisedTraitName[1]) +
    theme_bw() + theme(plot.title = element_text(size = 12), 
                       axis.title = element_text(size = 10))
  
  # Save it as an svg file.
  svg_file <- paste(outdir, gsub("/|#", "", output_name), ".svg", sep="")
  ggsave(filename = svg_file, plot = p, height = 4, width = 4.2)
  
} else # If fitting failed ...
{
  
  # Populate the vectors with missing values.
  B0_sch <- c(B0_sch, NA)
  E_sch <- c(E_sch, NA)
  E_D_sch <- c(E_D_sch, NA)	
  T_pk_est_sch <- c(T_pk_est_sch, NA)
  T_pk_sch <- c(T_pk_sch, NA)
  P_pk_sch <- c(P_pk_sch, NA)
}
}
  
  
  # Compile all data into a data frame.
  results <- data.frame(
    species, trait, E_sch, B0_sch, E_D_sch, 
    T_pk_est_sch, T_pk_sch, P_pk_sch
  )
  

# write the results out to a new file
write.csv(results, file = "../Results/summary_R.csv", row.names = FALSE)
