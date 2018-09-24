# This code is meant to fit thermal performance curves using various 
# thermal performance curve (TPC) models. 

# The main function starts after model functions have been loaded. You 
# we need to specify the Data file, the Model you want to run 
# ("Boltzmann", "Schoolfield" or "all" if we want to run both). You also 
# need to specify if you want to save plots (PLOT=TRUE) (one figure will be 
# plot for each model) or if you are running both models we can choose the 
# Overplot=TRUE selection in case that you want a unique figure with both 
# models overlaid.

# Depending on if you want to use the Schoolfield model with explicit
# T_pk parameter or not you have to set the SchoolTPk as TRUE or FALSE

# Requires package minpack.lm for the NLLS (which implements the more 
# robust Levenberg-Marqualdt algorithm) and ggplot2 for plotting, neither 
# of which is not part of the R standard libraries -- Please install if 
# necessary

#### Tref is specified at the beginning as GlobalEnvironment, so change 
# the value here in case you need it.

# Required packages...
library(ggplot2, graphics)
library(lattice)
library(minpack.lm)  ### The nls is now run under this package using the nlsLM instead of nls
library(sme)
library(truncnorm)

#############################
# F  U  N  C  T  I  O  N  S #
#############################

#### assign Tref as GlobalEnv
# Tref is the standardization temperature (in K). 
# This needs to be any value below the peak of the curve.
assign("Tref", 0 + 273.15, envir = .GlobalEnv) ## Set to Inf if you want to remove the normalization, so 1/Tref will be 0 

#### Estimate STARTING VALUES for the nls

GetE <- function(tmp, rate, T.p, k=8.62e-5)
    {
    # Estimate starting value for E, taking linear regression using the rise part
    # of the curve only.
    # ~~~ Parameters ~~~
    # tmp  : temperature data (in K).
    # rate : rate data corresponding to temperature above.
    # T.p  : temperature at which rate peaks, used as a cutoff point.
    # k    : Boltzmann constant.
    

  tmp.w <- which(tmp <= T.p)
  if (length(tmp.w) > 1 & length(unique(tmp[tmp.w]))>1)
  {
    m <- lm(log(rate[tmp.w]) ~ I(1 / (k * (tmp[tmp.w]))))
     return(abs(summary(m)$coefficients[2, c('Estimate', 'Std. Error')]))
  } else
  {
    return(c(0.7,2))  # Arbitrary estimate if we can't do regression.
  }

}

GetB0 <- function(tmp, rate)
{
    # Estimate starting value for the normalising constant.
    # ~~~ Parameters ~~~
    # tmp   : temperature data (in K).
    # rate  : rate data corresponding to temperature above.
    # T.ref : estimate normalising constant at this temperature (in K).

    if (min(tmp,na.rm=TRUE) > Tref)
        {
    return(log(min(rate[1],na.rm=TRUE)))
} else
    {
    return(log(max(rate[which(tmp <= Tref)],na.rm=TRUE)))
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




###################### Boltzmann - Arrhenius model.
Boltzmann.Arrhenius <- function(B0, E, temp) {
    
    # Boltzmann's constant. Units imply that E is in eV.
    k <- 8.62e-5 

    # B0 is the normalization constant.  
    # E is the activation energy.
    # Tref is the standardization temperature (in K).
    
    calc <- B0 - E/k * (1/temp - 1/Tref)

    return(calc)
}

###################### Schoolfield type models ######################

# Schoolfield function runs two different schoolfield models, one with
# explicit T_pk parameter (then we need to set SchoolTpk as TRUE),
#or the original one.

Schoolfield <- function(B0, E, E_D, T_h, temp, SchoolTpk=TRUE)
{ 
    # PARAMETERS/INPUTS (all temperatures in Kelvin) -
    # Boltzmann's constant. Units imply that E is in eV.
    k <- 8.62e-5 
    # temp   : temperature values to evaluate function at (single, scalar or vector of values)
    # B0     : Normalisation constant (log transformed)
    # E      : Activation energy (> 0)
    # E_D    : High temperature de-activation energy (> 0) 
    # Tref   : Standardization (reference) temperature; set to 0 if not wanted. To do    # this the whole term 1/Tref should be 0, so Tref has to be set to Inf
    # T_h (if SchoolTpk=TRUE): Temperature at which trait reaches peak value (Tpk)
    # T_h (if not Tpk)   : High temperature at which have of the enzyme units, on 
    # 				average are inactiveated.

    if (SchoolTpk==TRUE) # Sharpe-Schoolfield model with explicit T_pk parameter
        {
    return(B0 + log(exp((-E/k) * ((1/temp)-(1/Tref)))/(1 + (E/(E_D - E)) * exp(E_D/k * (1/T_h - 1/temp)))))

} else { # Original Sharpe-Schoolfield's model with low-temp inactivation term removed, slightly modified                                                                                                                     
      return(B0 + log(exp((-E/k) * ((1/temp) - (1/Tref)))/(1 + exp((E_D/k) * (1/T_h - 1/temp))))) }

}

####################################################################################
############################# M  A  I  N    C  O  D  E #############################
####################################################################################

TPCFit <- function(Data,PLOT=TRUE,OverPLOT=FALSE,Model="Schoolfield",SchoolTpk=TRUE,rand.st=FALSE, n.rand=20){ 
    # MODEL can be also "Boltzmann" or "all" for both
    # School can be TRUE (explicit Tpk parameter) or FALSE (Original Schoolfield)
    # rand.st    : boolean - use random starting values?
    # n.rand     : numeric - number of random starting values used (if
    #              `rand.st=TRUE`).

    if (OverPLOT==TRUE) PLOT <- FALSE


    # Loads Data
    curvespt <- Data
    curvespt$OriginalTraitValue <- as.numeric(curvespt$OriginalTraitValue)

    # Merges ConTemp and AmbientTemp in the same column
    NATemp <- which(is.na(curvespt$ConTemp))
    NAAmbientTemp <- which(!is.na(curvespt$AmbientTemp[NATemp]))
    curvespt$ConTemp[NATemp][NAAmbientTemp] <- curvespt$AmbientTemp[NATemp][NAAmbientTemp]
    curvespt$ConTemp <- as.numeric(curvespt$ConTemp)

    # Transform temperatures to Kelvin and log-transform the
    # trait values.
    curvespt$K <- curvespt$ConTemp + 273.15


    ###################################################################
    # Create unique species/individual IDs, with a series of vectors. #
    ###################################################################

    # Initialize the current ID at 0. You'll see why later...
    current_ID <- 0

    # Initialize an empty vector of species IDs (as a number).
    id <- c()

    # Initialize an empty vector of species names.
    id_spp <- c()

    # Initialize an empty vector of processes (e.g., photosynthesis, respiration).
    id_process <- c()

    # Read each row of the data frame.
    for (k in 1:nrow(curvespt))
        {

    # If the ID of this row is different from the current ID...
    if ( current_ID != curvespt$FinalID[k])
	{
    
    # Change the current ID to the one found in this row.
    current_ID <- curvespt$FinalID[k]
    
    # Add information for this species to the 3 vectors 
    # that we initialized above.
    id <- c(id, current_ID)
    id_spp <- c(id_spp, curvespt$Consumer[k])
    id_process <- c(id_process, curvespt$StandardisedTraitName[k])
}
}

    # Initialize empty vectors to store the parameter estimates
    # that will be obtained.

    if (Model=="Boltzmann" | Model=="all")
        {
    E_boltz <- c()
    B0_boltz <- c()
    T_pk_boltz <- c()
    P_pk_boltz <- c()
    AIC_boltz <- c()
    r_sq_boltz <- c()
}

    if (Model=="Schoolfield" | Model=="all")
        {
    B0_sch <- c()
    E_sch <- c()
    E_D_sch <- c()	
    T_h_sch <- c()
    T_pk_sch <- c()
    P_pk_sch <- c()
    AIC_sch <- c()
    r_sq_sch <- c()
    B0exp1_sch <- c()
    B0exp2_sch<- c()
    B0exp5_sch<- c()
    B0exp10_sch<- c()
}

    if (Model=="all"){
        selected_model<- c()}

    P_pkBug <- c()
    B0Bug <- c()
    AftPk <-c()
    BefPk <- c()


    # Go through every single species ID and try to fit the the selected models.
    for(i in 1:length(id)) 
        { 

    # Get only the part of the data that correspond to that particular ID.
    current_dataset <- curvespt[curvespt$FinalID == id[i],]

    ## ##  # If there are NA, remove
    ## if(any(is.na(current_dataset$OriginalTraitValue)) | any(is.na(current_dataset$K))){
    ## IS.NA <- which(is.na(current_dataset$OriginalTraitValue) | is.na(current_dataset$K))
    ## current_dataset <- current_dataset[-IS.NA,]}

    
    ## If there are negative values, substract the minimum value
    MinVal <- NA
    if (min(current_dataset$OriginalTraitValue,na.rm=TRUE)<=0){
        MinVal <- min(current_dataset$OriginalTraitValue)
        current_dataset$OriginalTraitValue <-current_dataset$OriginalTraitValue - MinVal
        current_dataset <-current_dataset[-which(current_dataset$OriginalTraitValue==0),]}

    ## Calculates number of data points after and before the Tpk
    Order <- order(current_dataset$K)
    Max <- which.max(current_dataset$OriginalTraitValue[Order])
    if (length(Max)>0){
    AftPk <- c(AftPk,(length(current_dataset$OriginalTraitValue)-Max))
    BefPk <- c(BefPk,(Max-1))} else {
    AftPk <- c(AftPk,NA)
    BefPk <- c(BefPk,NA)} 

    # Runs only if we have at least 5 data points and more than 5 diff temperatures
    if (length(unique(current_dataset$OriginalTraitValue))>=5 && length(unique(current_dataset$K))>=5)
        {

    # Estimate T.h as being approximately T.peak.
    T.h.st  <- GetTpk(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)
    E.st    <- GetE(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue, T.p=T.h.st)
    B.st <- GetB0(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)

    if (Model=="Boltzmann" | Model=="all"){
        ###############################
        # Boltzmann - Arrhenius model #
        ###############################

   

        # Initialize the fitting variable to NA (empty).
        boltzmann_nls <- NA

        
        # Try and fit the model.
        try( 
            boltzmann_nls <- nlsLM(
                log(OriginalTraitValue) ~ Boltzmann.Arrhenius(B0, E, temp = K),
                start = c(B0 = B.st, E = E.st),
                lower=c(B0=-Inf, E=0),
                upper=c(B0=Inf,  E=Inf),
                control=list(minFactor=1 / 2^16, maxiter=1e4),
                data = current_dataset, 
                na.action=na.omit),
            silent=TRUE
            )
        
        # If fitting worked ...
        if(!is.na(boltzmann_nls[1])) 
            { 

    # Collect the parameter estimates...
    if (!is.na(MinVal)){ ## Add MinVal if it was substracted
        B0_boltz <- c(B0_boltz, (coef(boltzmann_nls)["B0"]+MinVal))
        if (Model=="Boltzmann") B0Bug <- c(B0Bug,TRUE) ## When both models are fit, B0Bug is attached in Schoolfield
    }else {
         B0_boltz <- c(B0_boltz, coef(boltzmann_nls)["B0"])
         if (Model=="Boltzmann") B0Bug <- c(B0Bug,FALSE)
     }
    E_boltz <- c(E_boltz, coef(boltzmann_nls)["E"])
    AIC_boltz<- c(AIC_boltz, AIC(boltzmann_nls))

    # Calculate the R squared value as: 1 - (rss/tss)
    rss <- sum((exp(predict(boltzmann_nls)) - 
                    current_dataset$OriginalTraitValue)^2, 
               na.rm = TRUE)
    tss <- sum(
        (current_dataset$OriginalTraitValue - 
             mean(current_dataset$OriginalTraitValue, na.rm = TRUE))^2, 
        na.rm = TRUE)
    
    if ( tss != 0 )
        {
    r_sq_boltz <- c(r_sq_boltz, 1 - (rss/tss))
} else
    {
    r_sq_boltz <- c(r_sq_boltz, 1)
}
    
    # Calculate the peak of the curve and its 
    # corresponding temperature value.
    curr_prediction <- predict(boltzmann_nls)
    for (j in 1:length(curr_prediction))
        {
    
    # If we found the maximum performance, exit the loop.
    if (curr_prediction[j] == max(curr_prediction))
        {
    break
}
}
    
    T_pk_boltz <- c(T_pk_boltz, current_dataset$K[j])
    if (!is.na(MinVal)){ ## Add MinVal if it was substracted
        P_pk_boltz <- c(P_pk_boltz, (curr_prediction[j]+MinVal))
        if (Model=="Boltzmann") P_pkBug <- c(P_pkBug,TRUE) ## When both models are fit, PpkBug is attached in Schoolfield

    }else {
         P_pk_boltz <- c(P_pk_boltz, curr_prediction[j])
         if (Model=="Boltzmann") P_pkBug <- c(P_pkBug,FALSE)
     }


    #######################################
    # Plotting Boltzmann - Arrhenius' fit #
    #######################################
    
    # Create a name for the output file using:
    #	- the original id number
    #   - the species name
    #   - the model
    output_name <- paste(
        current_dataset$FinalID[1], 
        current_dataset$Consumer[1], 
        'Boltzmann_Arrhenius',
        sep = "_"
        )
    
    # Remove any characters that won't look good in a file name,
    # using a regular expression.
    output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
    
    # Convert spaces to underscores.
    output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
    
    # CHANGE THIS to set an alternative output directory.
    outdir <- "./"
    
    # Generate predictions from the model fit...
    tmp_temps <- seq(min(
        floor(current_dataset$K)), 
                     ceiling(max(current_dataset$K)
                             ), length = 200)
    
    tmp_model <- exp(Boltzmann.Arrhenius(
        coef(boltzmann_nls)["B0"],
        coef(boltzmann_nls)["E"],
        tmp_temps
        ))
    
    ModelToPlotB <- data.frame(
        Temperature = tmp_temps - 273.15, 
        TraitValue = tmp_model
        )
    
    # Prepare the data points of the original values.
    DataToPlot <- data.frame(
        Temperature = current_dataset$K - 273.15, 
        TraitValue = current_dataset$OriginalTraitValue
        )
    DataToPlot <- na.omit(DataToPlot)

    #### If we want individual plots
    if (PLOT==TRUE) {
        # Plot!
        p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                       y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                                   alpha = 0.7, pch = 21) + 
            geom_line(data = ModelToPlotB, 
                      aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
                      lwd = 1.3) +
                
                ggtitle(paste(current_dataset$Consumer[1])) +
                    xlab(expression(paste("Temperature (", degree, C, ")"))) + 
                        ylab(current_dataset$StandardisedTraitName[1]) +
                            theme_bw() + theme(plot.title = element_text(size = 12), 
                                               axis.title = element_text(size = 10)) +
                                annotate("text", size = 3, label=             
                                             paste("R^2 boltz=", sprintf("%.2f", r_sq_boltz[i]), "\nE boltz=", format(coef(boltzmann_nls)["E"], digits = 3),"\nAIC boltz=",format(AIC(boltzmann_nls),digits=3)), 
                                         x = min(DataToPlot[, "Temperature"]),
                                         y = mean(DataToPlot[, "TraitValue"]),
                                         hjust=0,
                                         fontface = 3)
        
        # Save it as an svg file.
        svg_file <- paste(outdir, gsub("/|#", "", output_name), ".svg", sep="")
            ggsave(filename = svg_file, plot = p, height = 4, width = 4.2)
    }
    
    
} else # If fitting failed ...
    {
    # Populate the vectors with missing values.
    E_boltz <- c(E_boltz, NA)
    B0_boltz <- c(B0_boltz, NA)
    T_pk_boltz <- c(T_pk_boltz, NA)
    P_pk_boltz <- c(P_pk_boltz, NA)
    AIC_boltz <- c(AIC_boltz, NA)
    r_sq_boltz <- c(r_sq_boltz, NA)
    if (Model=="Boltzmann") { ## When both models are fit, this is attached in Schoolfield
        B0Bug <- c(B0Bug,NA)
        P_pkBug <- c(P_pkBug,NA)}
}
    }

    if (Model=="Schoolfield" | Model=="all"){
        
    #####################
    # Schoolfield model #
    #####################
  

        if (rand.st)
            {
    # Create randomised starting points.
    E.st.pe <- E.st[1]  # Slope value.
    T.h.st  <- c(T.h.st, rnorm(n.rand-1, mean=T.h.st, sd=15))
    # We need truncated normal to ensure we don't get negative values of E.
    E.st   <- c(E.st.pe, rtruncnorm(n.rand-1, a=0, b=Inf, mean=E.st[1], sd=2 * E.st[2]))
    B.st   <- exp(B.st)
    # Randomise on linear scale. Again, we don't want negative rates.
    B.st <- c(B.st, log(rtruncnorm(n.rand-1, a=0, b=Inf, mean=B.st, sd=B.st / 2)))

    # We'll select the best model using AICc. Many of these turn out to be
    # similar.
    aics.out <- rep(NA, n.rand)

    for (h in 1:n.rand)
        {
    schoolfield_nls <- try(nlme(
        log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K,SchoolTpk=SchoolTpk),
       fixed = list(B0 + E + E_D + T_h ~ 1),
	random = pdDiag(B0 + E + E_D + T_h  ~ 1),
        method='ML',
        start=c(B0=B.st[h], E=E.st[h], E_D=4*E.st[h],T_h=T.h.st[h]),
        lower=c(B0=-Inf, E=0,  E_D=0,  T_h=250),
        upper=c(B0=Inf,  E=30, E_D=50, T_h=350),
        data=current_dataset,control=list(minFactor=1 / 2^16, maxiter=1024)),
                           silent=TRUE)
  

    if (class(schoolfield_nls) != 'try-error')
        {
    aics.out[h] <- AICc(schoolfield_nls)
}
}
    w <- which.min(aics.out)

    if (length(w) > 0)
        {
    schoolfield_nls <- try(nlme(
        log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K,SchoolTpk=SchoolTpk),
        fixed = list(B0 + E + E_D + T_h ~ 1),
	random = pdDiag(B0 + E + E_D + T_h  ~ 1),
        method='ML',
        start=c(B0=B.st[w], E=E.st[w], E_D=4*E.st[w],T_h=T.h.st[w]),
        lower=c(B0=-Inf, E=0,  E_D=0,  T_h=250),
        upper=c(B0=Inf,  E=30, E_D=50, T_h=350),
        data=current_dataset,control=list(minFactor=1 / 2^16, maxiter=1024)),
                          silent=TRUE)
    
} else
    {
    schoolfield_nls <- NA
}
} else
    {


    schoolfield_nls <- NA
    try( 
        schoolfield_nls <- nlme(
            log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K,SchoolTpk=SchoolTpk),
           fixed = list(B0 + E + E_D + T_h ~ 1),
            random = pdDiag(B0 + E + E_D + T_h  ~ 1),
            method='ML',
            start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
            lower=c(B0=-Inf,   E=0,    E_D=0, T_h=0),
            upper=c(B0=Inf,    E=Inf,  E_D=Inf, T_h=Inf),
            data=current_dataset, control=list(minFactor=1 / 2^16, maxiter=1024)),
        silent=TRUE)
    
}
          
 
    
    # If fitting worked ...
    if(!is.na(schoolfield_nls[1])) 
	{ 

    # Collect the parameter estimates..
    if (!is.na(MinVal)){ ## Add MinVal if it was substracted
        B0Bug <- c(B0Bug,TRUE)
        B0_sch <- c(B0_sch, (coef(schoolfield_nls)["B0"]+MinVal))
    }else {
         B0Bug <- c(B0Bug,FALSE)
         B0_sch <- c(B0_sch, coef(schoolfield_nls)["B0"])
     }
    E_sch <- c(E_sch, coef(schoolfield_nls)["E"])
    E_D_sch <- c(E_D_sch, coef(schoolfield_nls)["E_D"])
    T_h_sch <- c(T_h_sch, coef(schoolfield_nls)["T_h"])
    AIC_sch<- c(AIC_sch, AIC(schoolfield_nls))
    
    # Calculate the R squared value as: 1 - (rss/tss)
    rss <- sum((exp(predict(schoolfield_nls)) - 
                    current_dataset$OriginalTraitValue)^2, 
               na.rm = TRUE)
    tss <- sum(
        (current_dataset$OriginalTraitValue - 
             mean(current_dataset$OriginalTraitValue, na.rm = TRUE))^2, 
        na.rm = TRUE)
    
    if ( tss != 0 )
        {
    r_sq_sch <- c(r_sq_sch, 1 - (rss/tss))
} else
    {
    r_sq_sch <- c(r_sq_sch, 1)
}
    
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
    
    T_pk_sch <- c(T_pk_sch, current_dataset$K[j])
    if (!is.na(MinVal)){ ## Add MinVal if it was substracted
        P_pkBug <- c(P_pkBug,TRUE)
        P_pk_sch <- c(P_pk_sch, (curr_prediction[j]+MinVal))
    }else {
         P_pkBug <- c(P_pkBug,FALSE)
         P_pk_sch <- c(P_pk_sch, curr_prediction[j])
     }


    # Generate predictions for B0 at diff TºC from the model fit...
    tmp_temps <- seq(min(
        floor(1+273.15)), 
                     ceiling(10+273.15)
                             , by=1)
    
    tmp_model <- exp(Schoolfield(
        coef(schoolfield_nls)["B0"],
        coef(schoolfield_nls)["E"],
        coef(schoolfield_nls)["E_D"],
        coef(schoolfield_nls)["T_h"],
        tmp_temps
        ))

    B0exp1_sch <- c(B0exp1_sch,tmp_model[1])
    B0exp2_sch <- c(B0exp2_sch,tmp_model[2])
    B0exp5_sch <- c(B0exp5_sch,tmp_model[5])
    B0exp10_sch <- c(B0exp10_sch,tmp_model[10])    


    ##############################
    # Plotting Schoolfield's fit #
    ##############################
    
    # Create a name for the output file using:
    #	- the original id number
    #   - the species name
    #   - the model
    output_name <- paste(
        current_dataset$FinalID[1], 
        current_dataset$Consumer[1], 
        'Schoolfield',
        sep = "_"
        )
    
    
    # Remove any characters that won't look good in a file name,
    # using a regular expression.
    output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
    
    # Convert spaces to underscores.
    output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
    
    # CHANGE THIS to set an alternative output directory.
    outdir <- "./"
    
    # Generate predictions from the model fit...
    tmp_temps <- seq(min(
        floor(current_dataset$K)), 
                     ceiling(max(current_dataset$K)
                             ), length = 200)
    
    tmp_model <- exp(Schoolfield(
        coef(schoolfield_nls)["B0"],
        coef(schoolfield_nls)["E"],
        coef(schoolfield_nls)["E_D"],
        coef(schoolfield_nls)["T_h"],
        tmp_temps
        ))
    
    ModelToPlotS <- data.frame(
        Temperature = tmp_temps - 273.15, 
        TraitValue = tmp_model
        )
    
    # Prepare the data points of the original values.
    DataToPlot <- data.frame(
        Temperature = current_dataset$K - 273.15, 
        TraitValue = current_dataset$OriginalTraitValue
        )
    DataToPlot <- na.omit(DataToPlot)

    #### If we want individual plots
    if (PLOT==TRUE) {
        # Plot!
        p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                       y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                                   alpha = 0.7, pch = 21) + 
            geom_line(data = ModelToPlotS, 
                      aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
                      lwd = 1.3) +                           
                ggtitle(paste(current_dataset$Consumer[1])) +
                    xlab(expression(paste("Temperature (", degree, C, ")"))) + 
                        ylab(current_dataset$StandardisedTraitName[1]) +
                            theme_bw() + theme(plot.title = element_text(size = 12), 
                                               axis.title = element_text(size = 10)) +
                                annotate("text", size = 3, label=             
                                             paste("R^2","sch=", sprintf("%.2f", r_sq_sch[i]),"\nE sch=", format(coef(schoolfield_nls)["E"], digits = 3),"\nAIC sch=",format(AIC(schoolfield_nls),digits=3)), 
                                         x = min(DataToPlot[, "Temperature"]),
                                         y = mean(DataToPlot[, "TraitValue"]),
                                         hjust=0,
                                         fontface = 3)
        
        # Save it as an svg file.
        svg_file <- paste(outdir, gsub("/|#", "", output_name), ".svg", sep="")
            ggsave(filename = svg_file, plot = p, height = 4, width = 4.2)

    }
} else # If fitting failed ...
    {
    # Populate the vectors with missing values.
    B0_sch <- c(B0_sch, NA)
    B0Bug <- c(B0Bug,NA)
    E_sch <- c(E_sch, NA)
    E_D_sch <- c(E_D_sch, NA)	
    T_h_sch <- c(T_h_sch, NA)
    T_pk_sch <- c(T_pk_sch, NA)
    P_pk_sch <- c(P_pk_sch, NA)
    P_pkBug <- c(P_pkBug,NA)
    AIC_sch <- c(AIC_sch, NA)
    r_sq_sch <- c(r_sq_sch, NA)
    B0exp1_sch <- c(B0exp1_sch,NA)
    B0exp2_sch<- c(B0exp2_sch,NA)
    B0exp5_sch<- c(B0exp5_sch,NA)
    B0exp10_sch<- c(B0exp10_sch,NA)
}

}

    if (Model=="all")
        {

    if (OverPLOT==TRUE) { ## In case we want both models in the same figure with the overlaid fits.
        ##############################
        # Plotting both fits #
        ##############################
        
        # Create a name for the output file using:
        #	- the original id number
        #   - the species name
        #   - the model
        output_name <- paste(
            current_dataset$FinalID[1], 
            current_dataset$Consumer[1], 
            'Schoolfield_Boltz',
            sep = "_"
            )
        
        
        # Remove any characters that won't look good in a file name,
        # using a regular expression.
        output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
        
        # Convert spaces to underscores.
        output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
        
        # CHANGE THIS to set an alternative output directory.
        outdir <- "./"
        
        
        # Prepare the data points of the original values.
        DataToPlot <- data.frame(
            Temperature = current_dataset$K - 273.15, 
            TraitValue = current_dataset$OriginalTraitValue
            )
        DataToPlot <- na.omit(DataToPlot)

        ## Trait Units for plot
        Unit <- current_dataset$StandardisedTraitUnit[1]
        if (is.na(Unit)) Unit <- current_dataset$OriginalTraitUnit[1]
        
        # Plot!
        p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                       y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                                   alpha = 0.7, pch = 21) + 
            geom_line(data = ModelToPlotB, 
                      aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
                      lwd = 1.3) +            
                geom_line(data = ModelToPlotS, 
                          aes(x = Temperature, y = TraitValue), colour = "red", 
                          lwd = 1.3) +             
                    ggtitle(paste(current_dataset$Consumer[1])) +
                        xlab(expression(paste("Temperature (", degree, C, ")"))) + 
                            ylab(paste(current_dataset$StandardisedTraitName[1],"\n",Unit)) +
                                theme_bw() + theme(plot.title = element_text(size = 12), 
                                                   axis.title = element_text(size = 10)) +
                                    annotate("text", size = 3, label=             
                                                 paste("R^2","sch=", sprintf("%.2f", r_sq_sch[i]),"\nE sch=", format(coef(schoolfield_nls)["E"], digits = 3),"\nAIC sch=",format(AIC(schoolfield_nls),digits=3),"\nR^2 boltz=", sprintf("%.2f", r_sq_boltz[i]), "\nE boltz=", format(coef(boltzmann_nls)["E"], digits = 3),"\nAIC boltz=",format(AIC(boltzmann_nls),digits=3)), 
                                             x = min(DataToPlot[, "Temperature"]),
                                             y = mean(DataToPlot[, "TraitValue"]),
                                             hjust=0,
                                             fontface = 3)
        
        # Save it as an svg file.
        svg_file <- paste(outdir, gsub("/|#", "", output_name), ".svg", sep="")
            ggsave(filename = svg_file, plot = p, height = 4, width = 4.2)

    }     

    ##################################################################
    # Compare the two models using the Akaike Information Criterion. #
    ##################################################################
    
    # If both models failed to fit, add NA. 
    if (is.na(AIC_sch[i]) && is.na(AIC_boltz[i]))
	{
    selected_model <- c(selected_model, NA)
    
    # If only one of the two models could be fit, that is 
    # automatically the winner!
} else if (is.na(AIC_sch[i]) && !is.na(AIC_boltz[i]))
      {
    selected_model <- c(selected_model, 'boltzmann')
} else if (is.na(AIC_boltz[i]) && !is.na(AIC_sch[i]))
      {
    selected_model<- c(selected_model, "schoolfield")	
    
    # If both models were able to fit and Schoolfield's AIC
    # was lower, then that is the better model for this curve.
} else if (AIC_sch[i] < AIC_boltz[i])
      {
    selected_model<- c(selected_model,  "schoolfield")
    
    # And the opposite for the Boltzmann - Arrhenius model.
} else
    {
    selected_model<- c(selected_model,  "boltzmann")
}
}

}  else # If there are not enough values
    {
    
    # Populate the vectors with missing values.
    if (Model=="Schoolfield" | Model=="all"){
        B0_sch <- c(B0_sch, NA)
        B0Bug <- c(B0Bug,NA)
        E_sch <- c(E_sch, NA)
        E_D_sch <- c(E_D_sch, NA)	
        T_h_sch <- c(T_h_sch, NA)
        T_pk_sch <- c(T_pk_sch, NA)
        P_pk_sch <- c(P_pk_sch, NA)
        P_pkBug <- c(P_pkBug,NA)
        AIC_sch <- c(AIC_sch, NA)
        r_sq_sch <- c(r_sq_sch, NA)
        B0exp1_sch <- c(B0exp1_sch,NA)
        B0exp2_sch<- c(B0exp2_sch,NA)
        B0exp5_sch<- c(B0exp5_sch,NA)
        B0exp10_sch<- c(B0exp10_sch,NA)
    }    
    # Populate the vectors with missing values.
    if (Model=="Boltzmann" | Model=="all"){
        E_boltz <- c(E_boltz, NA)
        B0_boltz <- c(B0_boltz, NA)	
        T_pk_boltz <- c(T_pk_boltz, NA)
        P_pk_boltz <- c(P_pk_boltz, NA)
        AIC_boltz <- c(AIC_boltz, NA)
        r_sq_boltz <- c(r_sq_boltz, NA)}
    
    if (Model=="Boltzmann") { ## When both models are fit, this is attached in Schoolfield
        B0Bug <- c(B0Bug,NA)
        P_pkBug <- c(P_pkBug,NA)}
}
}


    ##################################################################
    # RESULTS FILE #
    ##################################################################
  

    if (Model=="all"){
        # Compile all data into a data frame.
        results <- data.frame(
            id, id_spp, id_process,E_boltz, E_sch, B0_boltz, B0_sch, B0Bug, E_D_sch, 
            T_h_sch, T_pk_boltz, T_pk_sch, P_pk_boltz, P_pk_sch, P_pkBug,AIC_boltz, 
            AIC_sch, r_sq_boltz, r_sq_sch, AftPk,BefPk,selected_model,B0exp1_sch,
            B0exp2_sch,B0exp5_sch,B0exp10_sch
            )}

    if (Model=="Schoolfield"){
        results <- data.frame(
            id, id_spp, id_process,E_sch, B0_sch, B0Bug, E_D_sch, 
            T_h_sch,  T_pk_sch, P_pk_sch, P_pkBug,
            AIC_sch,  r_sq_sch, AftPk,BefPk,B0exp1_sch,
            B0exp2_sch,B0exp5_sch,B0exp10_sch
            )}

    if (Model=="Boltzmann"){
        results <- data.frame(
            id, id_spp, id_process,E_boltz, B0_boltz, B0Bug,
            T_pk_boltz, P_pk_boltz, P_pkBug,
            AIC_boltz,  r_sq_boltz, AftPk,BefPk    )}


    # Write the results as a CSV file.
    write.csv(results, file = "results.csv", row.names = FALSE)

    return()
}

