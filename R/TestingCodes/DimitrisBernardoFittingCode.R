# Required packages...
library(ggplot2)
library(lattice)

#############################
# F  U  N  C  T  I  O  N  S #
#############################

# The Boltzmann - Arrhenius model.
Boltzmann.Arrhenius <- function(B0, E, temp) {
    
    # Boltzmann's constant. Units imply that E is in eV.
    k <- 8.62e-5 

    # B0 is the normalization constant.
    
    # E is the activation energy.

    # T_ref is the standardization temperature (in K). 
    # This needs to be any value below the peak of the 
    # curve.
    T_ref <- 273.15 + 10

    calc <- B0 * exp(-E*((1/(k*temp)) - (1/(k*T_ref))))

    return(log(calc))
}

# A version of the Schoolfield model which ignores 
# the inactivation of the rate-limiting enzyme at 
# low temperature.
Schoolfield <- function(B0, E, E_D, T_h, temp) {

    # T_ref is the standardization temperature (in K). 
    # This needs to be any value below the peak of the 
    # curve.
    T_ref <- 273.15 + 10
    
    # Boltzmann's constant. Units imply that E and E_D are in eV.
    k <- 8.62e-5

    # B0 is the normalization constant.
    
    # E is the activation energy.
    
    # E_D is the de-activation energy.
    
    # T_h is the temperature at which the rate-limiting enzyme 
    # is 50% active and 50% denatured due to high temperature.

    calc <- B0 * exp(-E * ((1/(k * temp)) - (1/(k * T_ref)))) / 
        (1 + exp(E_D/k * ((1/T_h) - (1/temp))))

    return(log(calc))
}


#### Estimate starting values for the nls

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
    m <- lm(log(rate[tmp.w]) ~ I(1 / (k * (tmp[tmp.w]))))
    return(abs(summary(m)$coefficients[2, 1]))
}

GetB0 <- function(tmp, rate, T.ref)
{
    # Estimate starting value for the normalising constant.
    # ~~~ Parameters ~~~
    # tmp   : temperature data (in K).
    # rate  : rate data corresponding to temperature above.
    # T.ref : estimate normalising constant at this temperature (in K).

    if (min(tmp) > T.ref)
        {
    return(min(rate[1]))
} else
    {
    return(max(rate[which(tmp < T.ref)]))
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


############################
# M  A  I  N    C  O  D  E #
############################

# CHANGE THIS LINE for R to read your input CSV file.
#curves<-read.csv("path/to/your/input/file.csv", stringsAsFactors = FALSE)
curvespt<-read.csv("/tmp/sheryl.csv", stringsAsFactors = FALSE)
curvespt <- read.csv("RichardSherylNewTemplate.csv")

# Remove negative and zero values before transformation.
# To take a logarithm you need positive values.
curvespt$OriginalTraitValue[curvespt$OriginalTraitValue <= 0] <- NA

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
    if ( current_ID != curvespt$OriginalID[k])
	{
    
    # Change the current ID to the one found in this row.
    current_ID <- curvespt$OriginalID[k]
    
    # Add information for this species to the 3 vectors 
    # that we initialized above.
    id <- c(id, current_ID)
    id_spp <- c(id_spp, curvespt$Consumer[k])
    id_process <- c(id_process, curvespt$StandardisedTraitName[k])
}
}

# Initialize empty vectors to store the parameter estimates
# that will be obtained.
E_boltz <- c()
B0_boltz <- c()
T_pk_boltz <- c()
P_pk_boltz <- c()
B0_sch <- c()
E_sch <- c()
E_D_sch <- c()	
T_h_sch <- c()
T_pk_sch <- c()
P_pk_sch <- c()
AIC_boltz <- c()
AIC_sch <- c()
r_sq_boltz <- c()
r_sq_sch <- c()
selected_model<- c()

# Go through every single species ID and try to fit the two models.
for(i in 1:length(id)) 
{ 

    # Get only the part of the data that correspond to that particular ID.
    current_dataset <- curvespt[curvespt$OriginalID == i,]

    # Estimate T.h as being approximately T.peak.
    T.h.st  <- GetTpk(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)
    E.st    <- GetE(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue, T.p=T.h.st)
    B.st <- GetB0(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue, T.ref=273.15 + 10)


    ###############################
    # Boltzmann - Arrhenius model #
    ###############################
    
    # Initialize the fitting variable to NA (empty).
    boltzmann_nls <- NA

    
    # Try and fit the model.
    try( 
        boltzmann_nls <- nls(
            log(OriginalTraitValue) ~ Boltzmann.Arrhenius(B0, E, temp = K), 
            start = c(B0 = B.st, E = E.st), 
            data = current_dataset, 
            na.action=na.omit
            )
	)
    
    # If fitting worked ...
    if(!is.na(boltzmann_nls[1])) 
	{ 

    # Collect the parameter estimates...   	
    E_boltz <- c(E_boltz, coef(boltzmann_nls)["E"])
    B0_boltz <- c(B0_boltz, coef(boltzmann_nls)["B0"])
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
    P_pk_boltz <- c(P_pk_boltz, curr_prediction[j])
    
    #######################################
    # Plotting Boltzmann - Arrhenius' fit #
    #######################################
    
    # Create a name for the output file using:
    #	- the original id number
    #   - the species name
    #   - the model
    output_name <- paste(
        current_dataset$OriginalID[1], 
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
    
    ModelToPlot <- data.frame(
        Temperature = tmp_temps - 273.15, 
        TraitValue = tmp_model
        )
    
    # Prepare the data points of the original values.
    DataToPlot <- data.frame(
        Temperature = current_dataset$K - 273.15, 
        TraitValue = current_dataset$OriginalTraitValue
        )
    DataToPlot <- na.omit(DataToPlot)
    
    # Plot!
    p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                   y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + 
        geom_line(data = ModelToPlot, 
                  aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
                  lwd = 1.3) + 
            ggtitle(paste(current_dataset$Consumer[1])) +
                xlab(expression(paste("Temperature (", degree, C, ")"))) + 
                    ylab(current_dataset$StandardisedTraitName[1]) +
                        theme_bw() + theme(plot.title = element_text(size = 12), 
                                           axis.title = element_text(size = 10)) +
                            annotate("text", size = element_text(size = 10), label = 
                                         paste("italic(R^2==", sprintf("%.2f", 1 - (rss/tss)), ")", sep = ""), 
                                     y = max(DataToPlot$TraitValue)  - 0.05 * 
                                         (max(DataToPlot$TraitValue) - 
                                              min(DataToPlot$TraitValue)), 
                                     x = max(DataToPlot$Temperature) - 0.85 * 
                                         (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
                                     parse = TRUE, fontface = 3)
    
    # Save it as an svg file.
    svg_file <- paste(outdir, gsub("/|#", "", output_name), ".svg", sep="")
        ggsave(filename = svg_file, plot = p, height = 4, width = 4.2)

    
} else # If fitting failed ...
    {
    
    # Populate the vectors with missing values.
    E_boltz <- c(E_boltz, NA)
    B0_boltz <- c(B0_boltz, NA)	
    T_pk_boltz <- c(T_pk_boltz, NA)
    P_pk_boltz <- c(P_pk_boltz, NA)
    AIC_boltz <- c(AIC_boltz, NA)
    r_sq_boltz <- c(r_sq_boltz, NA)
}

    #####################
    # Schoolfield model #
    #####################
    
    schoolfield_nls <- NA
    try( 
        schoolfield_nls <- nls(
            log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K), 
            start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
            
            data = current_dataset, 
            na.action=na.omit
            )
	)  
    
    # If fitting worked ...
    if(!is.na(schoolfield_nls[1])) 
	{ 

    # Collect the parameter estimates...
    B0_sch <- c(B0_sch, coef(schoolfield_nls)["B0"])
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
    P_pk_sch <- c(P_pk_sch, curr_prediction[j])
    
    ##############################
    # Plotting Schoolfield's fit #
    ##############################
    
    # Create a name for the output file using:
    #	- the original id number
    #   - the species name
    #   - the model
    output_name <- paste(
        current_dataset$OriginalID[1], 
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
    
    ModelToPlot <- data.frame(
        Temperature = tmp_temps - 273.15, 
        TraitValue = tmp_model
        )
    
    # Prepare the data points of the original values.
    DataToPlot <- data.frame(
        Temperature = current_dataset$K - 273.15, 
        TraitValue = current_dataset$OriginalTraitValue
        )
    DataToPlot <- na.omit(DataToPlot)
    
    # Plot!
    p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                   y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + 
        geom_line(data = ModelToPlot, 
                  aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
                  lwd = 1.3) + 
            ggtitle(paste(current_dataset$Consumer[1])) +
                xlab(expression(paste("Temperature (", degree, C, ")"))) + 
                    ylab(current_dataset$StandardisedTraitName[1]) +
                        theme_bw() + theme(plot.title = element_text(size = 12), 
                                           axis.title = element_text(size = 10)) +
                            annotate("text", size = element_text(size = 10), label = 
                                         paste("italic(R^2==", sprintf("%.2f", 1 - (rss/tss)), ")", sep = ""), 
                                     y = max(DataToPlot$TraitValue)  - 0.05 * 
                                         (max(DataToPlot$TraitValue) - min(DataToPlot$TraitValue)), 
                                     x = max(DataToPlot$Temperature) - 0.85 * 
                                         (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
                                     parse = TRUE, fontface = 3)
    
    # Save it as an svg file.
    svg_file <- paste(outdir, gsub("/|#", "", output_name), ".svg", sep="")
        ggsave(filename = svg_file, plot = p, height = 4, width = 4.2)
    
} else # If fitting failed ...
    {

    # Populate the vectors with missing values.
    B0_sch <- c(B0_sch, NA)
    E_sch <- c(E_sch, NA)
    E_D_sch <- c(E_D_sch, NA)	
    T_h_sch <- c(T_h_sch, NA)
    T_pk_sch <- c(T_pk_sch, NA)
    P_pk_sch <- c(P_pk_sch, NA)
    AIC_sch <- c(AIC_sch, NA)
    r_sq_sch <- c(r_sq_sch, NA)
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

# Compile all data into a data frame.
results <- data.frame(
    id, id_spp, id_process, E_boltz, E_sch, B0_boltz, B0_sch, E_D_sch, 
    T_h_sch, T_pk_boltz, T_pk_sch, P_pk_boltz, P_pk_sch, AIC_boltz, 
    AIC_sch, r_sq_boltz, r_sq_sch, selected_model
    )

# Write the results as a CSV file.
write.csv(results, file = "results.csv", row.names = FALSE)
