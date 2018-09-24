### Script to produce the plots of the fits and analyse
### how good the fits are.
rm(list=ls())
graphics.off()
library(ggplot2)
library(reshape2)

res_df <- read.csv("../Results/results.csv")
ids <- unique(as.character(res_df$Unique_id))
K <- 8.617 * 10^(-5)
## REQUIRED MULTIPLOT FUNCTION
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## Plots graph for each dataset
thermal_plots <- function(){
  pdf("../Results/Model_fits.pdf",10,10)
  for(i in ids){
    # Plot thermal curves for each dataset to a pdf.
    gaugo <- function(x){
      max_trait * exp((-E_gaugo * (x - T_pk_gaugo) * (x - T_pk_gaugo) - exp((E_D_gaugo * (x - T_pk_gaugo) - theta))))
    }
    # Define the school function with the parameters
    school <- function(x){
      B0 * exp(-E * ((1/(K*(x+273.15)) - (1/(K*283.15))))) / (1 + (E/(E_D - E)) * exp((E_D / K) * (1 / T_pk - 1 / (x+273.15))))
    }
    # Define the cubic function with the parameters
    cubic <- function(x){
      alpha + beta * x + gamma * (x^2) + epsilon * (x^3)
    }
    # Subset the data to contain all those of the same id
    sub <- subset(res_df,res_df$Unique_id==i)
    # Define the school parameters
    B0 <- sub$B0_school[1]
    E <- sub$E_school[1]
    E_D <- sub$E_D_school[1]
    T_pk <- sub$T_pk_school[1]
    # Define the gaugo parameters
    max_trait <- max(sub$Trait_Vals)
    E_gaugo <- sub$E_gaugo[1]
    E_D_gaugo <- sub$E_D_gaugo[1]
    T_pk_gaugo <- sub$T_pk_gaugo[1]
    theta <- sub$theta[1]
    alpha <- sub$alpha[1]
    beta <- sub$beta[1]
    gamma <- sub$gamma[1]
    epsilon <- sub$epsilon[1]
    x <- data.frame(seq(from=min(sub$Temp_Vals),to=max(sub$Temp_Vals),by=0.01))
    out <-data.frame(x,apply(X = x,MARGIN = 1,FUN = school),apply(X = x,MARGIN = 1,FUN = gaugo),apply(X = x,MARGIN = 1,FUN = cubic))
    names(out) <- c("x","school","gaugo","cubic")
    out <- melt(out,id="x")
    p <- ggplot()+
      geom_point(data=sub,aes(x=Temp_Vals,y=Trait_Vals))+
      geom_line(data = out,aes(x=x,y=value,colour=variable))+
      theme_bw()+
      scale_colour_discrete(name="Model",labels=c("Schoolfield Model","Gaussian-Gompertz Model","Cubic Model"))+
      ylab("Trait Value")+
      xlab("Temperature Value")+
      ggtitle(paste("Model fitting for",sub$Species_stand,", ",sub$Reference))
    print(p)
  }
  dev.off()
}


## Producing bar charts to give an indication of how good the different models were.
prod_bar <- function(){
  pdf("../Results/Figures/bar.pdf",5,4)
  school_tot <- length(unique(subset(res_df,res_df$School_success=="True")$Unique_id))
  school_good <- length(unique(subset(res_df,res_df$School_success=="True" & res_df$School_qual=="True")$Unique_id))
  both_total <- length(unique(subset(res_df,res_df$School_success=="True" & res_df$Gaugo_success=="True")$Unique_id)) 
  gaugo_tot <- length(unique(subset(res_df,res_df$Gaugo_success=="True")$Unique_id))
  gaugo_good <- length(unique(subset(res_df,res_df$Gaugo_success=="True" & res_df$Gaug_qual=="True")$Unique_id))
  cubic_tot <- length(unique(subset(res_df,res_df$Cubic_success=="True")$Unique_id))
  cubic_good <-length(unique(subset(res_df,res_df$Cubic_success=="True" & res_df$Cubic_qual=="True")$Unique_id))
  cubic_choose <- length(unique(subset(res_df,res_df$Choose_cubic=="True")$Unique_id))
  v <- data.frame(c(school_tot-school_good,school_good,gaugo_tot-gaugo_good,gaugo_good,cubic_tot-cubic_good,cubic_good),
                  c("Poor Fit","Good Fit","Poor Fit","Good Fit","Poor Fit","Good Fit"),
                  c("Schoolfield","Schoolfield","GauGo","GauGo","Cubic","Cubic"))
  names(v)<-c("Total","Level_of_Fit","Model")
  w <- data.frame(c(length(unique(subset(res_df,res_df$Model_name=="school")$Unique_id)),length(unique(subset(res_df,res_df$Model_name=="gaugo")$Unique_id)),cubic_choose))
  names(w) <- c("Total")
  p <- ggplot(v)+geom_bar(aes(x=Model,weight=Total,fill=Level_of_Fit)) + 
    scale_fill_manual(name="Level of Fit",guide=guide_legend(reverse=TRUE),
                      values=c("lightgreen","indianred","dodgerblue")) + 
    ylab("Number of Datasets") +
    theme_bw()
  print(p)
  dev.off()
  pdf("../Results/Figures/pie.pdf",5,5)
  pie(w$Total,labels=c("Schoolfield Model","Gaussian\n-Gompertz Model","Cubic Model"),col=c("indianred","lightgreen","dodgerblue"))
  dev.off()
}
r_distr <- function(){
  pdf("../Results/Figures/r_distr.pdf",7,3)
  data <- unique(subset(res_df,select=-c(Trait_Vals,Temp_Vals)))
  data <- melt(data = data,id.vars="Unique_id",
               measure.vars=c("R_Squared_gaugo","R_Squared_school","R_Squared_cubic"))
  p <- ggplot(data)+geom_density(aes(x = value,col=variable))+
    scale_colour_discrete(name="Model",labels=c("Gaussian-Gompertz Model","Schoolfield Model","Cubic Model"))+
    xlim(c(-1,1))+theme_bw()+xlab("R Squared Value")+ylab("Density")
  print(p)
  dev.off()
}

prep_example <- function(){
  pdf("../Results/Figures/prep_example.pdf",7,4)
  x <- data.frame(c(0,0.0013,0.022,0.3,2.5,5),c(0,5,10,15,20,25))
  n <- data.frame(c(30,35,40,45),c(4.9,4.2,1.1,0))
  names(x) <- c("y","x")
  names(n) <- c("x","y")
  samp <- function(x){
    exp(0.45*x-9)
  }
  p <- ggplot(x)+geom_point(aes(x=x,y=y,col="black"))+theme_bw()+
    geom_point(data=n,aes(x=x,y=y,col="red"))+
    stat_function(data=data.frame(seq(from=0,to=25,by=0.01)),fun=samp)+
    ylab("Log of Trait Value")+xlab("Temperature Value")+
    scale_y_log10(limits=c(0.00001,10))+
    scale_colour_manual(name="",labels=c("Included points","Excluded points"),values=c("black","red"))+
    annotate(geom="text",x=10,y=1.2,label="Gradient is E")+
    annotate(geom="text",x=8,y=0.00008,label="X intercept is B0")
  print(p)
  dev.off()
}
# Define the gaugo function with the parameters

poor_cubic <-function(){    # Plot thermal curves for each dataset to a pdf.
    gaugo <- function(x){
      max_trait * exp((-E_gaugo * (x - T_pk_gaugo) * (x - T_pk_gaugo) - exp((E_D_gaugo * (x - T_pk_gaugo) - theta))))
    }
    # Define the school function with the parameters
    school <- function(x){
      B0 * exp(-E * ((1/(K*(x+273.15)) - (1/(K*283.15))))) / (1 + (E/(E_D - E)) * exp((E_D / K) * (1 / T_pk - 1 / (x+273.15))))
    }
    # Define the cubic function with the parameters
    cubic <- function(x){
      alpha + beta * x + gamma * (x^2) + epsilon * (x^3)
    }
    pdf("../Results/Figures/poor_cubic.pdf",6,3.5)
    i <- "Ostreopsis siamensis-:-Morton et al. (1992)-:- 24.708-:- -81.125"
    # Subset the data to contain all those of the same id
    sub <- subset(res_df,res_df$Unique_id==i)
    
    # Define the school parameters
    B0 <- sub$B0_school[1]
    E <- sub$E_school[1]
    E_D <- sub$E_D_school[1]
    T_pk <- sub$T_pk_school[1]
    
    # Define the gaugo parameters
    max_trait <- max(sub$Trait_Vals)
    E_gaugo <- sub$E_gaugo[1]
    E_D_gaugo <- sub$E_D_gaugo[1]
    T_pk_gaugo <- sub$T_pk_gaugo[1]
    theta <- sub$theta[1]
    # Define the cubic parameters
    alpha <- sub$alpha[1]
    beta <- sub$beta[1]
    gamma <- sub$gamma[1]
    epsilon <- sub$epsilon[1]
    x <- data.frame(seq(from=min(sub$Temp_Vals),to=max(sub$Temp_Vals),by=0.01))
    out <-data.frame(x,apply(X = x,MARGIN = 1,FUN = school),apply(X = x,MARGIN = 1,FUN = gaugo),apply(X = x,MARGIN = 1,FUN = cubic))
    names(out) <- c("x","school","gaugo","cubic")
    out <- melt(out,id="x")
    p <- ggplot()+
      geom_point(data=sub,aes(x=Temp_Vals,y=Trait_Vals))+
      geom_line(data = out,aes(x=x,y=value,colour=variable))+
      theme_bw()+
      scale_colour_discrete(name="Model",labels=c("Schoolfield Model","Gaussian-Gompertz Model","Cubic Model"))+
      ylab("Trait Value")+
      xlab("Temperature Value")
    print(p)
    dev.off()
  }
good_curve <-function(){
    # Plot thermal curves for each dataset to a pdf.
    gaugo <- function(x){
      max_trait * exp((-E_gaugo * (x - T_pk_gaugo) * (x - T_pk_gaugo) - exp((E_D_gaugo * (x - T_pk_gaugo) - theta))))
    }
    # Define the school function with the parameters
    school <- function(x){
      B0 * exp(-E * ((1/(K*(x+273.15)) - (1/(K*283.15))))) / (1 + (E/(E_D - E)) * exp((E_D / K) * (1 / T_pk - 1 / (x+273.15))))
    }
    # Define the cubic function with the parameters
    cubic <- function(x){
      alpha + beta * x + gamma * (x^2) + epsilon * (x^3)
    }
    pdf("../Results/Figures/interest.pdf",7.4,3.3)
    
    i <- "Synechococcus-:-Mackey et al. (2013)-:- 15.310-:-  72.590"
    # Subset the data to contain all those of the same id
    sub <- subset(res_df,res_df$Unique_id==i)
    
    # Define the school parameters
    B0 <- sub$B0_school[1]
    E <- sub$E_school[1]
    E_D <- sub$E_D_school[1]
    T_pk <- sub$T_pk_school[1]
    
    # Define the gaugo parameters
    max_trait <- max(sub$Trait_Vals)
    E_gaugo <- sub$E_gaugo[1]
    E_D_gaugo <- sub$E_D_gaugo[1]
    T_pk_gaugo <- sub$T_pk_gaugo[1]
    theta <- sub$theta[1]
    # Define the cubic parameters
    alpha <- sub$alpha[1]
    beta <- sub$beta[1]
    gamma <- sub$gamma[1]
    epsilon <- sub$epsilon[1]
    x <- data.frame(seq(from=min(sub$Temp_Vals),to=max(sub$Temp_Vals),by=0.01))
    out <-data.frame(x,apply(X = x,MARGIN = 1,FUN = school),apply(X = x,MARGIN = 1,FUN = gaugo),apply(X = x,MARGIN = 1,FUN = cubic))
    names(out) <- c("x","school","gaugo","cubic")
    out <- melt(out,id="x")
    q <- ggplot()+
      geom_point(data=sub,aes(x=Temp_Vals,y=Trait_Vals))+
      geom_line(data = out,aes(x=x,y=value,colour=variable))+
      theme_bw()+ggtitle("A")+
      ylab("Trait Value")+
      xlab("Temperature Value")+theme(legend.position="none")
    j <- "Prorocentrum mexicanum-:-Morton et al. (1992)-:- 24.712-:- -81.125"
    # Subset the data to contain all those of the same id
    sub <- subset(res_df,res_df$Unique_id==j)
    
    # Define the school parameters
    B0 <- sub$B0_school[1]
    E <- sub$E_school[1]
    E_D <- sub$E_D_school[1]
    T_pk <- sub$T_pk_school[1]
    
    # Define the gaugo parameters
    max_trait <- max(sub$Trait_Vals)
    E_gaugo <- sub$E_gaugo[1]
    E_D_gaugo <- sub$E_D_gaugo[1]
    T_pk_gaugo <- sub$T_pk_gaugo[1]
    theta <- sub$theta[1]
    # Define the cubic parameters
    alpha <- sub$alpha[1]
    beta <- sub$beta[1]
    gamma <- sub$gamma[1]
    epsilon <- sub$epsilon[1]
    y <- data.frame(seq(from=min(sub$Temp_Vals),to=max(sub$Temp_Vals),by=0.01))
    out <-data.frame(y,apply(X = y,MARGIN = 1,FUN = school),apply(X = y,MARGIN = 1,FUN = gaugo),apply(X = y,MARGIN = 1,FUN = cubic))
    names(out) <- c("x","school","gaugo","cubic")
    out <- melt(out,id="x")
    p <- ggplot()+
      geom_point(data=sub,aes(x=Temp_Vals,y=Trait_Vals))+
      geom_line(data = out,aes(x=x,y=value,colour=variable))+
      theme_bw()+
      scale_colour_discrete(name="Model",labels=c("Schoolfield","Gaussian\n-Gompertz","Cubic"))+
      ylab("Trait Value")+ggtitle("B")+
      xlab("Temperature Value")
    multiplot(q,p,layout=matrix(c(1,2,2),nrow=1))
    dev.off()
  }
jumps <- function(){
  # Plot thermal curves for each dataset to a pdf.
  gaugo <- function(x){
    max_trait * exp((-E_gaugo * (x - T_pk_gaugo) * (x - T_pk_gaugo) - exp((E_D_gaugo * (x - T_pk_gaugo) - theta))))
  }
  # Define the school function with the parameters
  school <- function(x){
    B0 * exp(-E * ((1/(K*(x+273.15)) - (1/(K*283.15))))) / (1 + (E/(E_D - E)) * exp((E_D / K) * (1 / T_pk - 1 / (x+273.15))))
  }
  # Define the cubic function with the parameters
  cubic <- function(x){
    alpha + beta * x + gamma * (x^2) + epsilon * (x^3)
  }
  pdf("../Results/Figures/jumps.pdf",7.4,3.3)
  
  i <- "Skeletonema costatum-:-Suzuki and Takahashi (1995)-:- 35.500-:- 140.000"
  # Subset the data to contain all those of the same id
  sub <- subset(res_df,res_df$Unique_id==i)
  
  # Define the school parameters
  B0 <- sub$B0_school[1]
  E <- sub$E_school[1]
  E_D <- sub$E_D_school[1]
  T_pk <- sub$T_pk_school[1]
  
  # Define the gaugo parameters
  max_trait <- max(sub$Trait_Vals)
  E_gaugo <- sub$E_gaugo[1]
  E_D_gaugo <- sub$E_D_gaugo[1]
  T_pk_gaugo <- sub$T_pk_gaugo[1]
  theta <- sub$theta[1]
  # Define the cubic parameters
  alpha <- sub$alpha[1]
  beta <- sub$beta[1]
  gamma <- sub$gamma[1]
  epsilon <- sub$epsilon[1]
  x <- data.frame(seq(from=min(sub$Temp_Vals),to=max(sub$Temp_Vals),by=0.01))
  out <-data.frame(x,apply(X = x,MARGIN = 1,FUN = school),apply(X = x,MARGIN = 1,FUN = gaugo),apply(X = x,MARGIN = 1,FUN = cubic))
  names(out) <- c("x","school","gaugo","cubic")
  out <- melt(out,id="x")
  q <- ggplot()+
    geom_point(data=sub,aes(x=Temp_Vals,y=Trait_Vals))+
    geom_line(data = out,aes(x=x,y=value,colour=variable))+
    theme_bw()+scale_colour_discrete(name="Model",labels=c("Schoolfield","Gaussian\n-Gompertz","Cubic"))+
    ylab("Trait Value")+
    xlab("Temperature Value")
  print(q)
  dev.off()
}

## Function to run all plots for figures
run_all <- function(){
  print("Plotting figures ...")
  print("running poor_cubic")
  poor_cubic()
  print("running good_curve ...")
  good_curve()
  print("running jumps ...")
  jumps()
  print("running prep_example ...")
  prep_example()
  print("running prod_bar ...")
  prod_bar()
  print("running r_distr...")
  r_distr()
  print("running high_counts ...")
  high_counts()
  print("running remove_cubic ...")
  remove_cubic()
  print("Completed plotting of figures")
}
get_percentage <- function(){
  w <- data.frame(c(length(unique(subset(res_df,res_df$Model_name=="school")$Unique_id)),length(unique(subset(res_df,res_df$Model_name=="gaugo")$Unique_id)),cubic_choose))
  school <- w[1,]/661*100
  gaugo <- w[2,]/661*100
  cubic <- w[3,]/661*100
  print("School is :")
  print(school)
  print("Gaugo is :")
  print(gaugo)
  print("Cubic is :")
  print(cubic)
}

high_counts <- function(){
  y <- data.frame(res_df,as.data.frame(apply(X = res_df,1, FUN = function(x) sum(as.character(res_df$Unique_id)==x[1]))))
  names(y) <- c(names(res_df),"count")
  res_df2<-subset(y,y$count>18)
  cubic_choose <- length(unique(subset(res_df2,res_df2$Model_name=="cubic")$Unique_id))
  w <- data.frame(c(length(unique(subset(res_df2,res_df2$Model_name=="school")$Unique_id)),
                    length(unique(subset(res_df2,res_df2$Model_name=="gaugo")$Unique_id)),cubic_choose))
  names(w) <- c("Total")
  w
  pdf("../Results/Figures/high.pdf",5,5)
  pie(w$Total,labels=c("Schoolfield Model","Gaussian\n-Gompertz Model","Cubic Model"),col=c("indianred","lightgreen","dodgerblue"))
  dev.off()
}

remove_cubic <- function(){
  y <- length(unique(subset(res_df,res_df$AIC_GauGo>res_df$AIC_Schoolf)$Unique_id))
  z <- length(unique(subset(res_df,res_df$AIC_GauGo<res_df$AIC_Schoolf)$Unique_id))
  w <- data.frame(c(y,z))
  names(w)<-"Total"
  pdf("../Results/Figures/remove_cubic.pdf",5,5)
  pie(w$Total,labels=c("Schoolfield Model","Gaussian\n-Gompertz Model"),col=c("indianred","lightgreen"))
  dev.off()
}