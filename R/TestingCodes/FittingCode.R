# Code to analyze intraspecific thermal response curves of metabolic 
# traits.

# Requires package minpack.lm for the NLLS (which implements the more 
# robust Levenberg-Marqualdt algorithm) and ggplot2 for plotting, neither 
# of which is not part of the R standard libraries -- Please install if 
# necessary

# TODOs:
# * Use conditional dir.create() to create missing ODir (output 
#   directory) Results (see #below) and subdirectories such as Fig 

library(ggplot2, graphics) #for plotting
library(minpack.lm) # for NLLS

rm(list = ls())  # clear objects
graphics.off()  #close all open figures and graphics objects 

################# Scroll down to main Code first ###############

##########  ##########  ##########  ##########  ##########  ####
# Functions, including thermal Response models to fit to the data # 
##########  ##########  ##########  ##########  ##########  ####

Boltz <- function(const, E, T_ref, temp)
{ # Boltzmann-Arrhenius model
	
  # PARAMETERS/INPUTS (al temperatures in Kelvin) -
  # temp   : temperature values to evaluate function at (single, scalar or vector of values)
	# const  : Normalisation constant (log transformed)
  # E      : Activation energy (> 0)

	return(const + log(exp((-E/k) * ((1/temp) - (1/T_ref)))))

}

SchoolFTpk <- function(const, E, E_D, T_pk, T_ref, temp)
{ # Sharpe-Schoolfield model with explicit T_pk parameter
	
  # PARAMETERS/INPUTS (al temperatures in Kelvin) -
  # temp   : temperature values to evaluate function at (single, scalar or vector of values)
	# const  : Normalisation constant (log transformed)
  # E      : Activation energy (> 0)
  # E_D    : High temperature de-activation energy (> 0) 
  # T_ref  : Standardization (reference) temperature; set to 0 if not wanted   
  # T_pk   : Temperature at which trait reaches peak value

	return(const + log(exp((-E/k) * ((1/temp) - (1/T_ref)))/(1 + (E/(E_D - E)) * exp(E_D/k * (1/T_pk - 1/temp)))))

}

SchoolF <- function(const, E, E_D, T_pk, T_ref, temp)
{ # Original Sharpe-Schoolfield's model with low-temp inactivation term removed, slightly modified 

  # PARAMETERS/INPUTS -
  # temp   : temperature values to evaluate function at (single, scalar or vector of values)
	# const  : Normalisation constant (log transformed)
  # E      : Activation energy (> 0)
  # E_D    : High temperature de-activation energy (> 0) 
  # T_ref  : Standardization (reference) temperature; set to 0 if not wanted   
  # T_H    : High temperature at which deactivation sets in
  
  return(const + log(exp((-E/k) * ((1/temp) - (1/T_ref)))/(1 + exp((E_D/k) * (1/T_H - 1/temp))))) 
  
}

GetStartVals <- function(tmpData){ #Generates starting values for the NLLS algorithm
  
  T_pk_strt <- max(tmpData[which(tmpData[, "OriginalTraitValue"] == mean(max(tmpData[, "OriginalTraitValue"]), na.rm = TRUE)), 
                           "ConTemp"],na.rm=TRUE)  # Use temperature of peak trait value (flux) as start value for fitted T_pk
  
  tmpData_L <- subset(tmpData, tmpData[, "ConTemp"] <= T_pk_strt)  #extract data upto peak value
  tmpData_H <- subset(tmpData, tmpData[, "ConTemp"] >= T_pk_strt)  #extract data after peak value
  
  const_strt <- log(min(tmpData$OriginalTraitValue))  #Take rate at min temp as starting values of constant
  
  if (length(unique(tmpData_L$ConTemp)) > 1 ) { #if there are enough temperature points
    x <- 1/(k * tmpData_L$ConTemp)
    y <- log(tmpData_L$OriginalTraitValue)
    fit <- summary(lm(y ~ x))
    E_strt <- abs(fit$coefficients[2]) #Take estimated from Arrhenius plot as start value of E_A
  } else{ E_strt <- 0.7}
  
  if (length(unique(tmpData_H$ConTemp)) > 1 ) { #if there are enough temperature points
    x <- 1/(k * tmpData_H$ConTemp)
    y <- log(tmpData_H$OriginalTraitValue)
    fit <- summary(lm(y ~ x))
    E_D_strt <- abs(fit$coefficients[2]) #Take estimated from Arrhenius plot as start value of E_A
  } else{ E_D_strt <- E_strt*1.1}
  
  #    E_D_strt <- E_strt * 10
  
  return(list(const_strt = const_strt, E_strt = E_strt, E_D_strt = E_D_strt, T_pk_strt = T_pk_strt))
}

########## ########## ########## ########## ########## 
########## ########## Main Code ########## ########## 
########## ########## ########## ########## ########## 

IDir <- "../Data/"
ODir <- "../Results/"
dir.create(paste(ODir,'Figs',sep = ""),showWarnings = TRUE)

FileName <- "NewDataFilled.csv"

Data <- as.data.frame(read.table(paste(IDir, FileName, sep = ""), 
									 header = TRUE, sep = ",", strip.white = TRUE))

Data[, "ConTemp"] <- Data[, "ConTemp"] + 273.15  #convert to Kelvins

# assign Boltzmann constant (units of eV * K^-1) as a global:
assign("k", 8.617 * 10^-5, envir = .GlobalEnv)  

#contro <- nls.control(maxiter = 1000, tol = 1e-20, minFactor = 1/1024, printEval = FALSE, 
#    warnOnly = TRUE) #Control parameters for the NLS fitting 

OutRows <- length(unique(Data[, "ID"]))
OutNames <- c("ID", "Trait", "lnB0", "B0", "E", "E_D", "T_pk", 
"T_pkStart", "B_pk", "B_pkStart", "OriUnits", "StanUnits", "RSS", 
"R^2", "N_lessthanT_pk", "N_greaterthanT_pk","Habitat", "Latitude", 
"Longitude", "Species", "Citation", "FigureTable", "AcclimTemp", 
"AcclimTempDur", "AcclimTempDurUnits", "Input") #Headers for output

Results <- matrix(data = "", OutRows,length(OutNames),dimnames = list(NULL,OutNames)) #Preallocate empty matrix

tries <- 1000

IDs <- unique(Data$ID)

T_ref <- 0 # Note that this this is Kelvin

for (i in 1:length(IDs)){   # loop to run analysis for each thermal response separately
  
  tmpData <- subset(Data, ID == IDs[i])
  
	#  if(!any(is.na(tmpData$StandardisedTraitValue))){
	#    tmpData$OriginalTraitValue <- (tmpData$StandardisedTraitValue)
	#    tmpData$OriginalTraitUnit <- (tmpData$StandardisedTraitUnit)
	#  } 
    
  # If there are negative values, add the lowest value to each of the others, then remove that row
  AddedVal <- NA
  if(any(tmpData$OriginalTraitValue <= 0)){
    AddedVal <- min(tmpData$OriginalTraitValue)
    tmpData$OriginalTraitValue <- (tmpData$OriginalTraitValue - min(tmpData$OriginalTraitValue))
    tmpData <- tmpData[-which(tmpData$OriginalTraitValue == min(tmpData$OriginalTraitValue)),]
  }   

  #If the added value is recorded, can this then be used to the results table, to reposition B0, Bpk and Bpkstart
    
  StartVals <- GetStartVals(tmpData)
  
  # Some possible bounds on parameters:
	#
  #    const_bds <- c(-20,20)
  #    E_bds <- c(0,2)
  #    E_D_bds <- c(2,10)
  #    T_pk_bds <- c(-50+27ODir3.15,50+273.15)
  #    const_strt <- runif(tries,min = const_bds[1], max = const_bds[2])
  #    E_strt <- runif(tries,min = E_bds[1], max = E_bds[2])
  #    E_D_strt <- runif(tries,min = E_D_bds[1], max = E_D_bds[2])
  #    T_pk_strt <- runif(tries,min = T_pk_bds[1], max = T_pk_bds[2])
         
  NLSfit <- NULL #Initialize a null variable which will be replaced if tries succeed 
  for(j in 1:tries){
    #Schoolfield model (fitting using the Levenberg-Marqualdt algorithm)
    try(NLSfit <- nlsLM(log(OriginalTraitValue) ~ 
    SchoolF(const,E,E_D,T_pk,T_ref,ConTemp), tmpData,
						start = list(const=StartVals$const_strt, E = StartVals$E_strt, E_D = StartVals$E_D_strt, T_pk = StartVals$T_pk_strt), 
						control = nls.lm.control(maxiter = 1000,ftol = .Machine$double.eps, ptol = .Machine$double.eps,maxfev = 1000),
						lower = c(-Inf, 0, 0, 0), upper = c(Inf, Inf, Inf, 273.15 + 150)))

    if (!is.null(NLSfit)){break} else{StartVals$E_D_strt <- StartVals$E_D_strt*1.1}
  }
  
  ###### PLOT THE CURRENT DATASET ######
  DataToPlot <- data.frame(Temperature = tmpData[, "ConTemp"] - 273.15, OriginalTraitValue = tmpData[, "OriginalTraitValue"])
  
  p <- ggplot() + #Plot just the raw data
    geom_point(data = DataToPlot, aes(x = Temperature, y = OriginalTraitValue), size = I(3), colour = "blue", alpha = 0.7) + 
    xlab(expression(paste("Temperature (", degree, C, ")"))) + 
    ylab(paste(tmpData[1, "StandardisedTraitName"], sep = "")) 
  #      + 
  #      geom_vline(xintercept = tmpData[1, "AcclimFixTemp"], linetype = "longdash", colour = "blue")
  
  if (!is.null(NLSfit)){ # only if the nlls fitting converged
    ########Plot fitted model/curve#########
    temps <- 273.15 + seq(min(floor(tmpData$ConTemp - 273.15)), ceiling(max(tmpData$ConTemp - 273.15)), length = 200)
    model <- exp(SchoolF(coef(NLSfit)["const"], coef(NLSfit)["E"], coef(NLSfit)["E_D"], coef(NLSfit)["T_pk"], temps))
    ModelToPlot <- data.frame(Temperature = temps - 273.15, OriginalTraitValue = model)
    # calculate R^2 in linear scale         
    mod <- exp(SchoolF(coef(NLSfit)["const"], coef(NLSfit)["E"], coef(NLSfit)["E_D"], coef(NLSfit)["T_pk"], tmpData$ConTemp))
    RSS <- sum((tmpData$OriginalTraitValue - mod)^2)
    Results[i,"R^2"] <- 1 - (RSS / sum((tmpData$OriginalTraitValue - mean(tmpData$OriginalTraitValue))^2))
    #        Results[i,"R^2"] <- 1 - ( sum(SumStat$residuals^2) / sum((log(tmpData$OriginalTraitValue) - mean(log(tmpData$OriginalTraitValue)))^2)) #log scale R^2

    p <- p +  annotate("text", x = mean(DataToPlot[, "Temperature"]), y = mean(DataToPlot[, "OriginalTraitValue"]), 
             label = paste("E =", format(coef(NLSfit)["E"], digits = 2), "eV", ",\nT_pk =", format(coef(NLSfit)["T_pk"] - 273.15, digits = 2), "C", ",\nE_D =", format(coef(NLSfit)["E_D"], digits = 2), "eV", ",\nB0 =", format(exp(coef(NLSfit)["const"]), digits = 2), tmpData[1, "OriginalTraitUnit"], ",\nR^2 =", format(as.numeric(Results[i,"R^2"]), digits = 3))) 
             
    if (as.numeric(Results[i, "R^2"]) > 0 ){ # add the curve only if R^2 is positive
      p <- p + geom_line(data = ModelToPlot, aes(x = Temperature, y = OriginalTraitValue), colour = "red") 
    }            
    
    ###### Add to results table ######
    Results[i,"lnB0"] <- coef(NLSfit)["const"]
    Results[i,"B0"] <- exp(coef(NLSfit)["const"])
    Results[i,"E"] <- coef(NLSfit)["E"]
    Results[i,"T_pk"] <- coef(NLSfit)["T_pk"] - 273.15
    Results[i,"E_D"] <- coef(NLSfit)["E_D"]
    SumStat <- summary(NLSfit)
    Results[i,"RSS"] <- sum(SumStat$residuals^2)
  }
  
  ggsave(filename = paste(ODir,'../Results/Figs/', IDs[i], ".svg", sep = ""), plot = p, height = 4, width = 4.2)
  
  Results[i,"ID"] <- as.character(tmpData$ID[1])
  Results[i,"Habitat"] <- as.character(tmpData$Habitat[1])
  Results[i,"Species"] <- as.character(tmpData$Consumer[1])
  Results[i,"N_lessthanT_pk"] <- nrow(tmpData[tmpData[, "ConTemp"] < StartVals$T_pk_strt,])
  Results[i,"N_greaterthanT_pk"] <- nrow(tmpData[tmpData[, "ConTemp"] > StartVals$T_pk_strt,])
  Results[i,"Trait"] <- as.character(tmpData$StandardisedTraitName[1])s
  Results[i,"AcclimTemp"] <- as.character(tmpData$Labtemp[1])
  Results[i,"AcclimTempDur"] <- as.character(tmpData$Labtime[1])
  Results[i,"AcclimTempDurUnits"] <- as.character(tmpData$Labtimeunit[1])
  Results[i,"Latitude"] <- as.character(tmpData$Latitude[1])
  Results[i,"Longitude"] <- as.character(tmpData$Longitude[1])
  Results[i,"Citation"] <- as.character(tmpData$Citation[1])
  Results[i,"FigureTable"] <- as.character(tmpData$FigureTable[1])
  Results[i,"T_pkStart"] <- StartVals$T_pk_strt - 273.15
  Results[i,"OriUnits"] <- as.character(tmpData$OriginalTraitUnit[1])
  Results[i,"StanUnits"] <- as.character(tmpData$StandardisedTraitUnit[1])
  Results[i,"Input"] <- as.character(tmpData$Input[1])
  Results[i,"B_pkStart"] <- exp(tmpData[which.max(tmpData$ConTemp),which(colnames(tmpData)=="OriginalTraitValue")])
  Results[i,"B_pk"] <- max(ModelToPlot$OriginalTraitValue)
  
}

write.csv(Results, paste(ODir,"Results.csv",sep = ""), row.names = FALSE)

return
