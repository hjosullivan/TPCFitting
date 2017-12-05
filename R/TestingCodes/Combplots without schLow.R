#Creating PDFs with the regular plots (either Schoolfield or Boltzmann depending on which was selected)
#next to plots of the log of the trait vs 1/kT. Plots are created using both the original database and the
#results of the original analysis.


setwd("/Users/adamkhwaja/Dropbox")  #change to your wd

library(ggplot2)
library(lattice)
library(minpack.lm)  ### The nls is now run under this package using the nlsLM instead of nls
library(gridExtra)


curvespt<-read.csv("RichardSherylEditedTemplate.csv", stringsAsFactors = FALSE, fileEncoding = "latin1")
#reading original database, from which temperature and trait values will be used to construct plots


# Remove negative and zero values before transformation.
# To take a logarithm you need positive values.
curvespt$OriginalTraitValue[curvespt$OriginalTraitValue <= 0] <- NA

# Transform temperatures to Kelvin and log-transform the
# trait values.
curvespt$K <- curvespt$ConTemp + 273.15

k <- 8.62e-5


#Boltzmann function for use in later plotting

Boltzmann.Arrhenius <- function(B0,E,temp) {
  

  k <- 8.62e-5 
  
 
  T_ref <- 273.15 + 10
  
  calc <- B0 * exp(-E*((1/(k*temp)) - (1/(k*T_ref))))
  
  return(log(calc))
}


# A version of the Schoolfield model which ignores 
# the inactivation of the rate-limiting enzyme at 
# low temperature.
Schoolfield <- function(B0, E, E_D, T_h, temp) {
  
 
  T_ref <- 273.15 + 10
  
  
  k <- 8.62e-5
  
  calc <- B0 * exp(-E * ((1/(k * temp)) - (1/(k * T_ref)))) / 
    (1 + exp(E_D/k * ((1/T_h) - (1/temp))))
  
  return(log(calc))
}


# A version of the Schoolfield model which includes 
# the inactivation of the rate-limiting enzyme at 
# low temperature.
SchoolfieldLow <- function(B0, E, E_Dh, E_Dl, T_l,T_h, temp) {
  
  
  T_ref <- 273.15 + 10
  
  
  k <- 8.62e-5
  
  calc <- B0 * exp(-E * ((1/(k * temp)) - (1/(k * T_ref)))) / 
    (1 + exp(E_Dh/k * ((1/T_h) - (1/temp))) + exp(E_Dl/k * ((1/T_l) - (1/temp))))
  
  
  return(log(calc))
}




results <- read.csv("results.csv")  #reading results of previous analysis with estimates for the parameters
#and outcomes of model selection



for(i in 1:nrow(results))
  
{
  
  current_dataset <- results[results$id == i,]
  current_dataset_temp <- curvespt[curvespt$OriginalID == i,]  #loading both original database data
  #and the results data
  
  
  
  
  if(is.na(current_dataset$selected_model))  #if no model has been chosen, plot data with no model
    
  { 
    
    # Prepare the data points of the original values.
    DataToPlot <- data.frame(
      Temperature = current_dataset_temp$K - 273.15, 
      TraitValue = current_dataset_temp$OriginalTraitValue
    )
    DataToPlot <- na.omit(DataToPlot)
    
    
    
    
    # Plot!
    p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                                      y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + ggtitle(paste(current_dataset$id_spp)) +
      xlab(expression(paste("Temperature (", degree, C, ")"))) + 
      ylab(current_dataset_temp$StandardisedTraitName) +
      theme_bw() + theme(plot.title = element_text(size = 12), 
                         axis.title = element_text(size = 10))
    
    
    
    # CHANGE THIS to set an alternative output directory.
    outdir <- "/Users/adamkhwaja/"
    
    DataToPlot_log <- data.frame(
      InverseTemp = (1/((current_dataset_temp$K - 273.15)*k)), 
      LogTraitValue = (log(current_dataset_temp$OriginalTraitValue))
    )
    DataToPlot_log <- na.omit(DataToPlot_log)
    
    
    
    
    # Plot!
    q <- ggplot() + geom_point(data = DataToPlot_log, aes(x = InverseTemp, 
                                                          y = LogTraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + ggtitle(paste(current_dataset$id_spp)) +
      xlab(expression(paste("1/kT"))) + 
      ylab(expression(paste("log(trait value)"))) +
      theme_bw() + theme(plot.title = element_text(size = 12), 
                         axis.title = element_text(size = 10))
    
    
    output_name <- paste(
      current_dataset$id, 
      current_dataset$id_spp, 
      'Combplot_no_model',
      sep = "_"
    )
    
    # Remove any characters that won't look good in a file name,
    # using a regular expression.
    output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
    
    # Convert spaces to underscores.
    output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
    
    # CHANGE THIS to set an alternative output directory.
    outdir <- "/Users/adamkhwaja/"
    
    # Save it as an pdf file.
    pdf_file <- paste(outdir, gsub("/|#", "", output_name), ".pdf", sep="")
    
    
    pdf(pdf_file)
    
    grid.arrange(p, q, ncol = 2)
    
    dev.off()
    
  }
  
  
  
  ####BOLTZMANN####
  
  
  else if (current_dataset$selected_model == "boltzmann" )
    
  { 
    
    tmp_temps <- seq(min(
      floor(current_dataset_temp$K)), 
      ceiling(max(current_dataset_temp$K)
      ), length = 200)
    
    tmp_model <- exp(Boltzmann.Arrhenius(
      current_dataset$B0_boltz,
      current_dataset$E_boltz,
      tmp_temps
    ))
    
    ModelToPlot <- data.frame(
      Temperature = tmp_temps - 273.15, 
      TraitValue = tmp_model
    )
    
    # Prepare the data points of the original values.
    DataToPlot <- data.frame(
      Temperature = current_dataset_temp$K - 273.15, 
      TraitValue = current_dataset_temp$OriginalTraitValue
    )
    DataToPlot <- na.omit(DataToPlot)
    
    # Plot!
    p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                                      y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + 
      geom_line(data = ModelToPlot, 
                aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
                lwd = 1.3) + 
      ggtitle(paste(current_dataset$id_spp)) +
      xlab(expression(paste("Temperature (", degree, C, ")"))) + 
      ylab(current_dataset_temp$StandardisedTraitName) +
      theme_bw() + theme(plot.title = element_text(size = 12), 
                         axis.title = element_text(size = 10)) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(R^2==", sprintf("%.2f", current_dataset$r_sq_boltz), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.05 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(E==", sprintf("%.2f", current_dataset$E_boltz), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.1 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(B[0]==", sprintf("%.2f", current_dataset$B0_boltz), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.15 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(T[pk]==", sprintf("%.2f", current_dataset$T_pk_boltz), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.2 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.82 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) 
    
    
    
    #####creating the log temp vs 1/kT graph#####
    
    
    
    ModelToPlot_log <- data.frame(
      InverseTemp = (1/(k*(tmp_temps - 273.15))), 
      LogTraitValue = log(tmp_model)
    )
    
    DataToPlot_log <- data.frame(
      InverseTemp = (1/((current_dataset_temp$K)*k)), 
      LogTraitValue = (log(current_dataset_temp$OriginalTraitValue))
    )
    DataToPlot_log <- na.omit(DataToPlot_log)
    
    
    q <- ggplot() + geom_point(data = DataToPlot_log, aes(x = InverseTemp, 
                                                          y = LogTraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + 
      geom_line(data = ModelToPlot_log, 
                aes(x = InverseTemp, y = LogTraitValue), colour = "#1b9e77", 
                lwd = 1.3) + 
      ggtitle(paste(current_dataset$id_spp)) +
      xlab(expression(paste("1/kT"))) + 
      ylab(expression(paste("log(trait value)"))) +
      theme_bw() + theme(plot.title = element_text(size = 12), 
                         axis.title = element_text(size = 10))
    
    
    
    output_name <- paste(
      current_dataset$id, 
      current_dataset$id_spp, 
      'Combplot_boltzmann',
      sep = "_"
    )
    
    # Remove any characters that won't look good in a file name,
    # using a regular expression.
    output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
    
    # Convert spaces to underscores.
    output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
    
    # CHANGE THIS to set an alternative output directory.
    outdir <- "/Users/adamkhwaja/"
    
    # Save it as an pdf file.
    pdf_file <- paste(outdir, gsub("/|#", "", output_name), ".pdf", sep="")
    
    
    pdf(pdf_file)
    
    grid.arrange(p, q, ncol = 2)
    
    dev.off()
    
    
  }
  
  
  #####SCHOOLFIELD#####
  
  
  
  else if (current_dataset$selected_model == "schoolfield")
    
  {
    
    tmp_temps <- seq(min(
      floor(current_dataset_temp$K)), 
      ceiling(max(current_dataset_temp$K)
      ), length = 200)      
    
    tmp_model <- exp(Schoolfield(
      current_dataset$B0_sch,
      current_dataset$E_sch,
      current_dataset$E_D_sch,
      current_dataset$T_h_sch,
      tmp_temps
    ))
    
    
    ModelToPlot <- data.frame(
      Temperature = tmp_temps - 273.15, 
      TraitValue = tmp_model
    )
    
    
    DataToPlot <- data.frame(
      Temperature = current_dataset_temp$K - 273.15, 
      TraitValue = current_dataset_temp$OriginalTraitValue
    )
    DataToPlot <- na.omit(DataToPlot)
    
    # Plot!
    p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                                      y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + 
      geom_line(data = ModelToPlot, 
                aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
                lwd = 1.3) + 
      ggtitle(paste(current_dataset$id_spp)) +
      xlab(expression(paste("Temperature (", degree, C, ")"))) + 
      ylab(current_dataset_temp$StandardisedTraitName) +
      theme_bw() + theme(plot.title = element_text(size = 12), 
                         axis.title = element_text(size = 10)) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(R^2==", sprintf("%.2f", current_dataset$r_sq_sch), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.05 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(E==", sprintf("%.2f", current_dataset$E_sch), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.1 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(E[D]==", sprintf("%.2f", current_dataset$E_D_sch), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.15 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(B[0]==", sprintf("%.2f", current_dataset$B0_sch), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.2 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(T[pk]==", sprintf("%.2f", current_dataset$T_pk_sch), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.25 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.82 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) 
    
    
    
    #####temp vs 1/kT#####
    
    ModelToPlot_log <- data.frame(
      InverseTemp = (1/(k*(tmp_temps - 273.15))), 
      LogTraitValue = log(tmp_model)
    )
    
    DataToPlot_log <- data.frame(
      InverseTemp = (1/((current_dataset_temp$K)*k)), 
      LogTraitValue = (log(current_dataset_temp$OriginalTraitValue))
    )
    DataToPlot_log <- na.omit(DataToPlot_log)
    
    
    q <- ggplot() + geom_point(data = DataToPlot_log, aes(x = InverseTemp, 
                                                          y = LogTraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + 
      geom_line(data = ModelToPlot_log, 
                aes(x = InverseTemp, y = LogTraitValue), colour = "#1b9e77", 
                lwd = 1.3) + 
      ggtitle(paste(current_dataset$id_spp)) +
      xlab(expression(paste("1/kT"))) + 
      ylab(expression(paste("log(trait value)"))) +
      theme_bw() + theme(plot.title = element_text(size = 12), 
                         axis.title = element_text(size = 10))
    
    
    ###name and save as two plots###
    
    
    output_name <- paste(
      current_dataset$id, 
      current_dataset$id_spp, 
      'Combplot_schoolfield',
      sep = "_"
    )
    
    # Remove any characters that won't look good in a file name,
    # using a regular expression.
    output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
    
    # Convert spaces to underscores.
    output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
    
    # CHANGE THIS to set an alternative output directory.
    outdir <- "/Users/adamkhwaja/"
    
    # Save it as an pdf file.
    pdf_file <- paste(outdir, gsub("/|#", "", output_name), ".pdf", sep="")
    
    
    pdf(pdf_file)
    
    grid.arrange(p, q, ncol = 2)
    
    dev.off()
    
    
  }
  
  
  ###schoolfieldLow###
  
  
  else if(current_dataset$selected_model == "schoolfieldLow")
    
  { 
    
    
    tmp_temps <- seq(min(
      floor(current_dataset_temp$K)), 
      ceiling(max(current_dataset_temp$K)
      ), length = 200)      
    
    tmp_model <- exp(Schoolfield(
      current_dataset$B0_sch,
      current_dataset$E_sch,
      current_dataset$E_D_sch,
      current_dataset$T_h_sch,
      tmp_temps
    ))
    
    
    ModelToPlot <- data.frame(
      Temperature = tmp_temps - 273.15, 
      TraitValue = tmp_model
    )
    
    
    DataToPlot <- data.frame(
      Temperature = current_dataset_temp$K - 273.15, 
      TraitValue = current_dataset_temp$OriginalTraitValue
    )
    DataToPlot <- na.omit(DataToPlot)
    
    # Plot!
    p <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
                                                      y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + 
      geom_line(data = ModelToPlot, 
                aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
                lwd = 1.3) + 
      ggtitle(paste(current_dataset$id_spp)) +
      xlab(expression(paste("Temperature (", degree, C, ")"))) + 
      ylab(current_dataset_temp$StandardisedTraitName) +
      theme_bw() + theme(plot.title = element_text(size = 12), 
                         axis.title = element_text(size = 10)) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(R^2==", sprintf("%.2f", current_dataset$r_sq_sch), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.05 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(E==", sprintf("%.2f", current_dataset$E_sch), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.1 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(E[D]==", sprintf("%.2f", current_dataset$E_D_sch), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.15 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(B[0]==", sprintf("%.2f", current_dataset$B0_sch), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.2 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.85 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) +
      annotate("text", size = element_text(size = 10), label = 
                 paste("italic(T[pk]==", sprintf("%.2f", current_dataset$T_pk_sch), ")", sep = ""), 
               y = max(DataToPlot$TraitValue)  - 0.25 * 
                 (max(DataToPlot$TraitValue) - 
                    min(DataToPlot$TraitValue)), 
               x = max(DataToPlot$Temperature) - 0.82 * 
                 (max(DataToPlot$Temperature) - min(DataToPlot$Temperature)),
               parse = TRUE, fontface = 3) 
    
    
    
    #####temp vs 1/kT#####
    
    ModelToPlot_log <- data.frame(
      InverseTemp = (1/(k*(tmp_temps - 273.15))), 
      LogTraitValue = log(tmp_model)
    )
    
    DataToPlot_log <- data.frame(
      InverseTemp = (1/((current_dataset_temp$K)*k)), 
      LogTraitValue = (log(current_dataset_temp$OriginalTraitValue))
    )
    DataToPlot_log <- na.omit(DataToPlot_log)
    
    
    q <- ggplot() + geom_point(data = DataToPlot_log, aes(x = InverseTemp, 
                                                          y = LogTraitValue), size = 3, col = "black", bg = "lightcyan2", 
                               alpha = 0.7, pch = 21) + 
      geom_line(data = ModelToPlot_log, 
                aes(x = InverseTemp, y = LogTraitValue), colour = "#1b9e77", 
                lwd = 1.3) + 
      ggtitle(paste(current_dataset$id_spp)) +
      xlab(expression(paste("1/kT"))) + 
      ylab(expression(paste("log(trait value)"))) +
      theme_bw() + theme(plot.title = element_text(size = 12), 
                         axis.title = element_text(size = 10))
    
    
    ###name and save as two plots###
    
    
    output_name <- paste(
      current_dataset$id, 
      current_dataset$id_spp, 
      'Combplot_schoolfield',
      sep = "_"
    )
    
    # Remove any characters that won't look good in a file name,
    # using a regular expression.
    output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
    
    # Convert spaces to underscores.
    output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
    
    # CHANGE THIS to set an alternative output directory.
    outdir <- "/Users/adamkhwaja/"
    
    # Save it as an pdf file.
    pdf_file <- paste(outdir, gsub("/|#", "", output_name), ".pdf", sep="")
    
    
    pdf(pdf_file)
    
    grid.arrange(p, q, ncol = 2)
    
    dev.off()
    
  }
  
} 



