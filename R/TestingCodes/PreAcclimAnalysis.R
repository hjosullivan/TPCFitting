#########################################################################################################################
##################### Preliminar acclimation analysis for Richard&Sheryl  ###############################
#########################################################################################################################
library(stringr)


setwd("~/ImperialCollege/Ecoinformatics/Database/working-datasets/GlobalDataset/")
load("GlobalDataset.Rdata")

Data <- GlobalDataset[which(GlobalDataset$SubmittedBy=="Alex" | GlobalDataset$SubmittedBy=="Richard" | GlobalDataset$SubmittedBy=="Sheryl"),]


###### Selects subset of data using unique species, unique citation and unique trait

Rate <- paste(Data$ConSpecies,Data$StandardisedTraitName,Data$Citation,sep="_")

UniqueSp <- unique(Rate)

Size <- rep(NA,length(UniqueSp))
TempData <- rep(NA,length(UniqueSp))
Acclimation <- rep(NA,length(UniqueSp))
LabGrowthTemp <- rep(NA,length(UniqueSp))
AmbientTemp <- rep(NA,length(UniqueSp))
UniqueCurvePerSp <- rep(NA,length(UniqueSp))
for (i in 1:length(UniqueSp)){

    WHICH <- which(Rate==UniqueSp[i])
    if (length(WHICH>0)) {
        if (any(!is.na(Data$ConSize[WHICH]))) Size[i] <- TRUE else Size[i] <- FALSE
        Curves <- Data$OriginalID[WHICH]
        if (length(unique(Curves))>1) {
            UniqueCurvePerSp[i] <- length(unique(Curves))
            DifTemp <- rep(NA,length(unique(Curves)))
            Acc <- rep(NA,length(unique(Curves)))
            for (j in 1:length(unique(Curves))){
                Acc1 <- unique(Data$ConAccTemp[which(Data$OriginalID==unique(Curves)[j])])
                Acc2 <- unique(Data$ConAccFixTemp[which(Data$OriginalID==unique(Curves)[j])])
                if (length(which(!is.na(unique(Acc1))))>0) {Acc[j] <- TRUE} else if  (length(which(!is.na(unique(Acc2))))>0) {Acc[j] <- TRUE}  else {Acc[j] <- FALSE}

                Temp <- Data$ConTemp[which(Data$OriginalID==unique(Curves)[j])]
                if (length(unique(Temp))>3) DifTemp[j] <- TRUE else DifTemp[j] <- FALSE
            }
            if(length(which(Acc!=FALSE))>1) Acclimation[i] <- TRUE else Acclimation[j] <- FALSE
            LabGr <- rep(NA,length(unique(Curves)))
            for (j in 1:length(unique(Curves))){
                LabGr1 <- unique(Data$LabGrowthTemp[which(Data$OriginalID==unique(Curves)[j])])
                
                if (any(!is.na(unique(LabGr1))))  {LabGr[j] <- TRUE} else {LabGr[j] <- FALSE}
            }
            
            if(length(which(LabGr==TRUE))>1) LabGrowthTemp[i] <- TRUE else LabGrowthTemp[j] <- FALSE
            AmbTemp <- rep(NA,length(unique(Curves)))
            for (j in 1:length(unique(Curves))){
                AmbTemp1 <- unique(Data$AmbientTemp[which(Data$OriginalID==unique(Curves)[j])])
                
                if (any(!is.na(unique(AmbTemp1))))  {AmbTemp[j] <- TRUE} else {AmbTemp[j] <- FALSE}
            }
            
            if(length(which(AmbTemp==TRUE))>1) AmbientTemp[i] <- TRUE else AmbientTemp[j] <- FALSE
            
            if(length(which(DifTemp==TRUE))>1) TempData[i] <- TRUE else TempData[j] <- FALSE
        }} else {
            TempData[i] <- NA
            UniqueCurvePerSp[i] <- NA
            Size[i] <- NA
            Acclimation[i] <- NA
            LabGrowthTemp[i] <- NA
            AmbientTemp[i] <- NA
        }
}


DataForTemp <- cbind(Acclimation,LabGrowthTemp,AmbientTemp)
IsData <- rep(NA,length=nrow(DataForTemp))
for (i in 1:nrow(DataForTemp))
    {
        IsData[i] <- length(which(DataForTemp[i,]==TRUE))}

Accli <- which(TempData==TRUE & !is.na(UniqueCurvePerSp) & IsData>1)

AccliSubSet <- Data[which(!is.na(match(Rate,UniqueSp[Accli]))),]


OKData <- rep(NA,length=length(UniqueSp))
for (i in Accli)
    {
        WHICH <- which(Rate==UniqueSp[i])
        Lab <- unique(Data$LabGrowthTemp[WHICH])
        LabT <- unique(Data$LabGrowthDur[WHICH])
        Acc <- unique(Data$ConAccTemp[WHICH])
        AccT <- unique(Data$ConAccTime[WHICH])
        AccFix <- unique(Data$ConAccFixTemp[WHICH])
        AccFixT <- unique(Data$ConAccFixTime[WHICH])
        AmbientTemp <- unique(Data$AmbientTemp[WHICH])
        if (!is.na(Acc)){            
            if (!is.na(Lab)){
                if (Lab!=Acc) {OKData[i] <- TRUE}
           } else if (!is.na(AmbientTemp)){
                 if (AmbientTemp!=Acc) {OKData[i] <- TRUE}}

        } else if (!is.na(AccFix)){            
            if (!is.na(Lab)){
                if (Lab!=AccFix) {OKData[i] <- TRUE}
           } else if (!is.na(AmbientTemp)){
            if (AmbientTemp!=AccFix) {OKData[i] <- TRUE}}

        } else if (!is.na(AmbientTemp)) {
              if (!is.na(Lab)) {
                  if (AmbientTemp!=Lab) {OKData[i] <- TRUE}}
 
        } else  {OKData[i] <- FALSE}
    }


AccliDataToFit <- Data[which(!is.na(match(Rate,UniqueSp[which(OKData==TRUE)]))),]
save(AccliDataToFit,file="~/ImperialCollege/Ecoinformatics/Database/working-datasets/Richard-Alex/AccliDataToFit.Rdata")


#### Loads results from SchoolArrheniusFit
results <- read.csv("~/ImperialCollege/Ecoinformatics/Database/working-datasets/Richard-Alex/resultsNorm10.csv")

SpData <- unique(Data$ConSpecies[which(!is.na(match(Rate,UniqueSp[which(OKData==TRUE)])))])
results$id_spp <- SpData[match( as.character(results$id_spp),unique(as.character(results$id_spp)))]  ## Renames species name according to the original data

ID <- rep(NA,length=length(SpData))
Tpk2 <- rep(NA,length=length(SpData))
for (i in 1:length(SpData))
    {
    WHICH <- which(!is.na(match(results$id_spp,SpData[i])))
    ## We have several curves for the same species, but we take only the curves for which we have AccFixTemp
    Is.Na <- which(!is.na(AccliDataToFit$ConAccFixTemp[match(results$id[WHICH],AccliDataToFit$OriginalID)]))
    if (length(Is.Na)>0){WHICH <- WHICH[Is.Na]}
    Model <- c(results$AIC_boltz[WHICH],results$AIC_sch[WHICH])
    Row <- WHICH[which(results[WHICH,]==Model[which.min(Model)],arr.ind=TRUE)[1]]
    Col <- which(results[WHICH,]==Model[which.min(Model)],arr.ind=TRUE)[2]
    if (Col==16) Tpk2[i] <- results$T_pk_sch[Row]
    if (Col==15) Tpk2[i] <- results$T_pk_boltz[Row]
    ID[i] <- results$id[Row]
}

### Calculates Tpk2 (AccliTemp)
Tpk2 <- Tpk2-273.15

MatchID <- match(ID,AccliDataToFit$OriginalID)
AccliDataToFit$AmbientTemp[MatchID]
AccliDataToFit$LabGrowthTemp[MatchID]
AccliDataToFit$ConAccFixTemp[MatchID]

### Calculates Tpk1 (LabGrowthTemp)
Tpk1 <- AccliDataToFit$LabGrowthTemp[MatchID]

### T1 is equal to Tpk (LabGrowthTemp) and T2 is the ConAccFixTemp
T2 <- AccliDataToFit$ConAccFixTemp[MatchID]
T1 <- Tpk1

### When we have Ambient Temperature but not Acclimation Temperature, I've seen that Ambient Temperature could be use as the initial temperature (Effects of root temperature on growth and photosynthesisin conifer seedlingsduring shoot elongation) and we use LabGrowthTemp as Acclimation Temperature(T2) because the curve was measured at this temperature so...

AmbTemp <- which(!is.na(AccliDataToFit$AmbientTemp[MatchID]))
T1[AmbTemp] <- AccliDataToFit$AmbientTemp[MatchID][AmbTemp]
T2[AmbTemp] <- AccliDataToFit$LabGrowthTemp[MatchID][AmbTemp]

Tpk1[AmbTemp] <- T1[AmbTemp]

### Acclimation Magnitude assuming that Tpk1 is the closest temperature to T1
AccliMagnitude <- (Tpk2-Tpk1)/(T2-T1)


### Acclimation Time (Assuming that t1 is time at growth temp and if this is NA, we take 0 instead)
t1 <- AccliDataToFit$LabGrowthDur[MatchID]
t1[which(AccliDataToFit$LabGrowthDurUnit[MatchID]=="weeks")] <- t2[which(AccliDataToFit$LabGrowthDurUnit[MatchID]=="weeks")]*7

t2 <- AccliDataToFit$ConAccFixTime[MatchID]


### t2 is LabGrowthTime when the curve is measured at the growth temp and t1 is 0
t2[AmbTemp] <- t1[AmbTemp]
t1[AmbTemp] <- 0

AccliRate <- AccliMagnitude/(t2-t1)


###### Calulates acclimation just selecting the thermal curve with the lower AIC; now we can compare different curves without taking into account the AIC
for (i in 1:length(SpData))
    {
    WHICH <- which(results$id_spp==SpData[i])
    Is.Na <- which(!is.na(AccliDataToFit$ConAccFixTemp[match(results$id[WHICH],AccliDataToFit$OriginalID)]))
    if (length(Is.Na)>0){WHICH <- WHICH[Is.Na]}

    T2 <- rep(NA,length=length(WHICH))
    T1 <- rep(NA,length=length(WHICH))
    t1 <- rep(NA,length=length(WHICH))
    t2 <- rep(NA,length=length(WHICH))
    Tpk2 <- rep(NA,length=length(WHICH))
    Tpk1 <- rep(NA,length=length(WHICH))
    for (j in 1:length(WHICH)){
        if (!is.na(results$selected_model[WHICH][j])) {
            if (as.character(results$selected_model)[WHICH][j]=="schoolfield") Tpk2[j] <- results$T_pk_sch[WHICH][j] else Tpk2[j] <- results$T_pk_boltz[WHICH][j]
        } else Tpk2[j] <- NA

        
        T2[j] <- unique(AccliDataToFit$LabGrowthTemp[which(!is.na(match(AccliDataToFit$OriginalID,results$id[WHICH][j])))])
        t2[j] <- unique(AccliDataToFit$LabGrowthDur[which(!is.na(match(AccliDataToFit$OriginalID,results$id[WHICH][j])))])
        if (!is.na(t2[j])){
        if (unique(AccliDataToFit$LabGrowthDur[which(!is.na(match(AccliDataToFit$OriginalID,results$id[WHICH][j])))])=="weeks") t2[j] <- t2[j]*7
        if (unique(AccliDataToFit$LabGrowthDur[which(!is.na(match(AccliDataToFit$OriginalID,results$id[WHICH][j])))])=="months") t2[j] <- t2[j]*30}

        
        T1[j] <- unique(AccliDataToFit$ConAccFixTemp[which(!is.na(match(AccliDataToFit$OriginalID,results$id[WHICH][j])))])
        t1[j] <- unique(AccliDataToFit$ConAccFixTime[which(!is.na(match(AccliDataToFit$OriginalID,results$id[WHICH][j])))])
        if (!is.na(t1[j])){
                 if (unique(AccliDataToFit$ConAccFixTime[which(!is.na(match(AccliDataToFit$OriginalID,results$id[WHICH][j])))])=="weeks") t1[j] <- t1[j]*7
                 if (unique(AccliDataToFit$ConAccFixTime[which(!is.na(match(AccliDataToFit$OriginalID,results$id[WHICH][j])))])=="months") t1[j] <- t1[j]*30}

            if (is.na(T1[j])) T1[j] <- unique(AccliDataToFit$AmbientTemp[which(!is.na(match(AccliDataToFit$OriginalID,results$id[WHICH][j])))])
        
        
    }


#### Loads results from SchoolArrheniusFit
results <- read.csv("~/ImperialCollege/Ecoinformatics/Database/working-datasets/Richard-Alex/resultsNorm5.csv")

Data$OriginalID <- as.numeric(Data$OriginalID)
    
    #### Converts time to seconds

    Data$LabGrowthDur <- as.numeric(Data$LabGrowthDur)

    Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="hours")] <-     Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="hours")]/24
    Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="month")] <-     Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="month")]*30
    Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="months")] <-     Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="months")]*30
    Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="weeks")] <-     Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="weeks")]*7
    Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="week")] <-     Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="week")]*7    
    Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="year")] <-     Data$LabGrowthDur[which(Data$LabGrowthDurUnit=="year")]*356    

    Data$ConAccTime[which(Data$ConAccTimeUnit=="hours")] <-     Data$ConAccTime[which(Data$ConAccTimeUnit=="hours")]/24
    Data$ConAccTime[which(Data$ConAccTimeUnit=="h")] <-     Data$ConAccTime[which(Data$ConAccTimeUnit=="h")]/24
    Data$ConAccTime[which(Data$ConAccTimeUnit=="month")] <-     Data$ConAccTime[which(Data$ConAccTimeUnit=="month")]*30
    Data$ConAccTime[which(Data$ConAccTimeUnit=="week")] <-     Data$ConAccTime[which(Data$ConAccTimeUnit=="week")]*7

    Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="hour")] <-     Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="hour")]/24
    Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="hours")] <-     Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="hours")]/24
    Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="minutes")] <-     (Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="minutes")]/24)/60
    Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="months")] <-     Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="months")]*30
    Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="weeks")] <-     Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="weeks")]*7
    Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="week")] <-     Data$ConAccFixTime[which(Data$ConAccFixTimeUnit=="week")]*7    


    ##### For each species calculate acclimation

    UniqueSp <- unique(Data$ConSpecies)
    DataP <- paste(Data$ConSpecies,Data$StandardisedTraitName,Data$Citation,Data$Figure)
    UniqueDataP <- unique(DataP)

    E1 <-c()
    Tpk1 <- c()
    T01 <- c()
    Ppk1 <- c()
    B01 <-c()
    T1 <- c()
    ID1 <- c()
    I1 <- c()
    t1 <- c()
    
    for (i in 1:length(UniqueDataP))
        {
    WHICH <- which(DataP==UniqueDataP[i])
    IS.NA <- which(is.na(Data$ConAccTemp[WHICH]))
    if (length(IS.NA)>0) WHICH <- WHICH[-IS.NA]
    if(length(unique(Data$ConAccTemp[WHICH]))>1) {

        IDx <- Data$OriginalID[WHICH][match(unique(Data$ConAccTemp[WHICH]),Data$ConAccTemp[WHICH])]
        Ix <- rep(i,length(Data$OriginalID[WHICH][match(unique(Data$ConAccTemp[WHICH]),Data$ConAccTemp[WHICH])]))
        Tx  <- Data$ConAccTemp[WHICH][match(IDx,Data$OriginalID[WHICH])]
        tx  <- Data$ConAccTime[WHICH][match(IDx,Data$OriginalID[WHICH])]
        if(any(!is.na(tx))) tx[1] <- 0 
        
        Tpk <- rep(NA,length(IDx))
        Ppk <-rep(NA,length(IDx))
        B0 <- rep(NA,length(IDx))
        E0 <- rep(NA,length(IDx))
        
        for (j in 1:length(IDx)){
            if (!is.na(as.character(results$selected_model[IDx][j]))){
                
                if (as.character(results$selected_model[IDx][j])=="schoolfield"){
                    if (results$r_sq_sch[IDx][j]>0.5){ ## excludes low R^2
                    Tpk[j] <- results$T_pk_sch[IDx][j]
                    Ppk[j] <- results$P_pk_sch[IDx][j]
                    B0[j] <- results$B0_sch[IDx][j]
                    E0[j] <- results$E_sch[IDx][j]}}

                if (as.character(results$selected_model[IDx][j])=="boltzmann"){
                    if (results$r_sq_boltz[IDx][j]>0.5){
                    Tpk[j] <- results$T_pk_boltz[IDx][j]
                    Ppk[j] <- results$P_pk_boltz[IDx][j]
                    B0[j] <- results$B0_boltz[IDx][j]
                    E0[j] <- results$E_boltz[IDx][j]}}
            }
        }
        Tpk1 <- c(Tpk1,Tpk)
        Ppk1 <- c(Ppk1,Ppk)
        B01 <- c(B01,B0)
        E1 <- c(E1,E0)
        
        ID1 <- c(ID1,IDx)
        I1 <- c(I1,Ix)
        T1 <- c(T1,Tx)
        t1 <- c(t1,tx)
        
    }
}



           WholeDiff <- function(x){
               x <- as.vector(x)
               Dif<- x-x[1]
               return(Dif)
           }

    ### Removes those TPC for which B0 is too high (this is because we don't have enough data points before the Tpk if we look at the plots)
    BadFit <- which(B01>150)
    I1 <-I1[-BadFit]
    ID1 <-ID1[-BadFit]
    Tpk1 <-Tpk1[-BadFit]
    Ppk1 <-Ppk1[-BadFit]
    E1 <- E1[-BadFit]        
    B01 <- B01[-BadFit]
    T1 <-T1[-BadFit]
    t1 <- t1[-BadFit]
    

    ### Calculates acclimation as the difference between values at each different acclimation temperature and the base temperature
    At1 <- as.vector(unlist(lapply(split(t1,I1),WholeDiff)))
    ATa1 <- as.vector(unlist(lapply(split(T1,I1),WholeDiff)))
    ATpk1 <- as.vector(unlist(lapply(split(Tpk1,I1),WholeDiff)))
    APpk1 <-as.vector(unlist(lapply(split(Ppk1,I1),WholeDiff)))
    AB01 <-as.vector(unlist(lapply(split(B01,I1),WholeDiff)))
    AE1 <- as.vector(unlist(lapply(split(E1,I1),WholeDiff)))
   
    Remove <- which(ATa1==0)
    ATa1 <-ATa1[-Remove]
    ATpk1 <-ATpk1[-Remove]
    APpk1 <-APpk1[-Remove]
    AB01 <-AB01[-Remove]
    AE1 <-AE1[-Remove]
    IdNumber1 <- I1[-Remove]
    OriginalID1 <- ID1[-Remove]
    At1 <- At1[-Remove]

    

    
    ### For those UniqueDataP we don't have ConAccTemp, search for LabGrowthTemp
    UniqueDataP2 <- 1:length(UniqueDataP)
    UniqueDataP2 <- UniqueDataP2[-unique(IdNumber1)]

    E2 <-c()
    Tpk2 <- c()
    T02 <- c()
    Ppk2 <- c()
    B02 <-c()
    T2 <- c()
    ID2 <- c()
    I2 <- c()
    t2 <- c()
    
    for (i in UniqueDataP2)
        {
    WHICH <- which(DataP==UniqueDataP[i])
    IS.NA <- which(is.na(Data$LabGrowthTemp[WHICH]))
    if (length(IS.NA)>0) WHICH <- WHICH[-IS.NA]
    if(length(unique(Data$LabGrowthTemp[WHICH]))>1) {

        IDx <- Data$OriginalID[WHICH][match(unique(Data$LabGrowthTemp[WHICH]),Data$LabGrowthTemp[WHICH])]
        Ix <- rep(i,length(Data$OriginalID[WHICH][match(unique(Data$LabGrowthTemp[WHICH]),Data$LabGrowthTemp[WHICH])]))
        Tx  <- Data$LabGrowthTemp[WHICH][match(IDx,Data$OriginalID[WHICH])]
        tx  <- Data$LabGrowthDur[WHICH][match(IDx,Data$OriginalID[WHICH])]
        if(any(!is.na(tx))) tx[1] <- 0 
        
        Tpk <- rep(NA,length(IDx))
        Ppk <-rep(NA,length(IDx))
        B0 <- rep(NA,length(IDx))
        E0 <- rep(NA,length(IDx))
        
        for (j in 1:length(IDx)){
            if (!is.na(as.character(results$selected_model[IDx][j]))){
                
                if (as.character(results$selected_model[IDx][j])=="schoolfield"){
                    if (results$r_sq_sch[IDx][j]>0.5){ ## excludes low R^2
                        Tpk[j] <- results$T_pk_sch[IDx][j]
                        Ppk[j] <- results$P_pk_sch[IDx][j]
                        B0[j] <- results$B0_sch[IDx][j]
                        E0[j] <- results$E_sch[IDx][j]}}

                if (as.character(results$selected_model[IDx][j])=="boltzmann"){
                    if (results$r_sq_boltz[IDx][j]>0.5){
                        Tpk[j] <- results$T_pk_boltz[IDx][j]
                        Ppk[j] <- results$P_pk_boltz[IDx][j]
                        B0[j] <- results$B0_boltz[IDx][j]
                        E0[j] <- results$E_boltz[IDx][j]}}

            }
        }
        Tpk2 <- c(Tpk2,Tpk)
        Ppk2 <- c(Ppk2,Ppk)
        B02 <- c(B02,B0)
        E2 <- c(E2,E0)
        
        ID2 <- c(ID2,IDx)
        I2 <- c(I2,Ix)
        T2 <- c(T2,Tx)
        t2 <- c(t2,tx)
        
    }
}

    ### Removes those TPC for which B0 is too high (this is because we don't have enough data points before or after the Tpk if we look at the plots)
    BadFit <- which(B02>150)
    I2 <-I2[-BadFit]
    ID2 <-ID2[-BadFit]
    Tpk2 <-Tpk2[-BadFit]
    Ppk2 <-Ppk2[-BadFit]
    E2 <- E2[-BadFit]        
    B02 <- B02[-BadFit]
    T2 <-T2[-BadFit]
    t2 <- t2[-BadFit]
    

    ### Calculates acclimation as the difference between values at each different acclimation temperature and the base temperature
    At2 <- as.vector(unlist(lapply(split(t2,I2),WholeDiff)))
    ATa2 <- as.vector(unlist(lapply(split(T2,I2),WholeDiff)))
    ATpk2 <- as.vector(unlist(lapply(split(Tpk2,I2),WholeDiff)))
    APpk2 <-as.vector(unlist(lapply(split(Ppk2,I2),WholeDiff)))
    AB02 <-as.vector(unlist(lapply(split(B02,I2),WholeDiff)))
    AE2 <- as.vector(unlist(lapply(split(E2,I2),WholeDiff)))
   
    Remove <- which(ATa2==0)
    ATa2 <-ATa2[-Remove]
    ATpk2 <-ATpk2[-Remove]
    APpk2 <-APpk2[-Remove]
    AB02 <-AB02[-Remove]
    AE2 <-AE2[-Remove]
    IdNumber2 <- I2[-Remove]
    OriginalID2 <- ID2[-Remove]
    At2 <- At2[-Remove]


    ### For those UniqueDataP we don't have ConAccTemp neither LabGrowthTemp, search for ConAccFixTemp
    UniqueDataP3 <- 1:length(UniqueDataP)
    UniqueIdNumber <- sort(c(unique(IdNumber1),unique(IdNumber2)))
    UniqueDataP3 <- UniqueDataP3[-UniqueIdNumber]

    

    E3 <-c()
    Tpk3 <- c()
    T03 <- c()
    Ppk3 <- c()
    B03 <-c()
    T3 <- c()
    ID3 <- c()
    I3 <- c()
    t3 <- c()
    
    for (i in UniqueDataP3)
        {
    WHICH <- which(DataP==UniqueDataP[i])
    IS.NA <- which(is.na(Data$ConAccFixTemp[WHICH]))
    if (length(IS.NA)>0) WHICH <- WHICH[-IS.NA]
    if(length(unique(Data$ConAccFixTemp[WHICH]))>1) {

        IDx <- Data$OriginalID[WHICH][match(unique(Data$ConAccFixTemp[WHICH]),Data$ConAccFixTemp[WHICH])]
        Ix <- rep(i,length(Data$OriginalID[WHICH][match(unique(Data$ConAccFixTemp[WHICH]),Data$ConAccFixTemp[WHICH])]))
        Tx  <- Data$ConAccFixTemp[WHICH][match(IDx,Data$OriginalID[WHICH])]
        tx  <- Data$ConAccFixTime[WHICH][match(IDx,Data$OriginalID[WHICH])]
        if(any(!is.na(tx))) tx[1] <- 0 
        
        
        Tpk <- rep(NA,length(IDx))
        Ppk <-rep(NA,length(IDx))
        B0 <- rep(NA,length(IDx))
        E0 <- rep(NA,length(IDx))
        
        for (j in 1:length(IDx)){
            if (!is.na(as.character(results$selected_model[IDx][j]))){
                
                if (as.character(results$selected_model[IDx][j])=="schoolfield"){
                    if (results$r_sq_sch[IDx][j]>0.5){ ## excludes low R^2
                        Tpk[j] <- results$T_pk_sch[IDx][j]
                        Ppk[j] <- results$P_pk_sch[IDx][j]
                        B0[j] <- results$B0_sch[IDx][j]
                        E0[j] <- results$E_sch[IDx][j]}}

                if (as.character(results$selected_model[IDx][j])=="boltzmann"){
                    if (results$r_sq_boltz[IDx][j]>0.5){
                        Tpk[j] <- results$T_pk_boltz[IDx][j]
                        Ppk[j] <- results$P_pk_boltz[IDx][j]
                        B0[j] <- results$B0_boltz[IDx][j]
                        E0[j] <- results$E_boltz[IDx][j]}}

            }
        }
        Tpk3 <- c(Tpk3,Tpk)
        Ppk3 <- c(Ppk3,Ppk)
        B03 <- c(B03,B0)
        E3 <- c(E3,E0)
        
        ID3 <- c(ID3,IDx)
        I3 <- c(I3,Ix)
        T3 <- c(T3,Tx)
        t3 <- c(t3,tx)
    }
}


     ### Calculates acclimation as the difference between values at each different acclimation temperature and the base temperature
    At3 <- as.vector(unlist(lapply(split(t3,I3),WholeDiff)))
    ATa3 <- as.vector(unlist(lapply(split(T3,I3),WholeDiff)))
    ATpk3 <- as.vector(unlist(lapply(split(Tpk3,I3),WholeDiff)))
    APpk3 <-as.vector(unlist(lapply(split(Ppk3,I3),WholeDiff)))
    AB03 <-as.vector(unlist(lapply(split(B03,I3),WholeDiff)))
    AE3 <- as.vector(unlist(lapply(split(E3,I3),WholeDiff)))
   
    Remove <- which(ATa3==0)
    ATa3 <-ATa3[-Remove]
    ATpk3 <-ATpk3[-Remove]
    APpk3 <-APpk3[-Remove]
    AB03 <-AB03[-Remove]
    AE3 <-AE3[-Remove]
    IdNumber3 <- I3[-Remove]
    OriginalID3 <- ID3[-Remove]
    At3 <- At3[-Remove]

    
    

    ATa <- c(ATa1,ATa2,ATa3)
    ATpk <- c(ATpk1,ATpk2,ATpk3)
    APpk <- c(APpk1,APpk2,APpk3)
    AB0 <- c(AB01,AB02,AB03)
    AE <- c(AE1,AE2,AE3)
    At <- c(At1,At2,At3)
    TotalID <- c(OriginalID1,OriginalID2,OriginalID3)
    
    Sp <- Data$Consumer[match(UniqueDataP[sort(c(IdNumber1,IdNumber2,IdNumber3))],DataP)]
    Habitat <-Data$Habitat[match(UniqueDataP[sort(c(IdNumber1,IdNumber2,IdNumber3))],DataP)]
    Trait <- Data$StandardisedTraitName[match(UniqueDataP[sort(c(IdNumber1,IdNumber2,IdNumber3))],DataP)]

    lapply(lapply(split(Sp,Trait),unique),length)
   
    

    ############################### PLOTS ########################
    par(mfrow=c(2,2))

    ## ATPK~ATA
    Is.NA <- which(is.na(ATpk) | is.na(ATa))
    x <- ATa[-Is.NA]
    y <- ATpk[-Is.NA]
    plot(y~x,ylab="",xlab="",cex.axis=0.75,pch=21,cex=.4)
    Modelo <- lm(y~x)
    A <-round(summary(Modelo)$coefficients[2,1],3)
    B <-round(summary(Modelo)$coefficients[1,1],3)
    C <-round(summary(Modelo)$r.squared,3)
    D <- round(summary(Modelo)$coefficients[2,4],3)
    Eq <- paste("y= ",A,"x +(",B,")",sep="")
    
    xx <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=length(x))
    x1Predict <- xx
    Predict <- predict(Modelo,newdata=data.frame(x1=x1Predict))
    lines(Predict~x,col="red",lwd=2)

    legend("bottomright",c(Eq,substitute(paste(R^{2},"= ",C,sep=""),list(C=C)),expression(italic(p)*"-value= 0.001")),cex=0.75,horiz=FALSE,bty="n")
    mtext(expression(paste(Delta,"Tpk", sep="")),2,line=2.4,cex=0.75,font=1)
    mtext(expression(paste(Delta,"Ta", sep="")),1,line=2.4,cex=0.75,font=1)

    ## AE~ATa
    Is.NA <- which(is.na(AE) | is.na(ATa))
    x <- ATa[-Is.NA]
    y <- AE[-Is.NA]
    plot(y~x,ylab="",xlab="",cex.axis=0.75,pch=21,cex=.4)
    Modelo <- lm(y~x)
    A <-round(summary(Modelo)$coefficients[2,1],3)
    B <-round(summary(Modelo)$coefficients[1,1],3)
    C <-round(summary(Modelo)$r.squared,3)
    D <- round(summary(Modelo)$coefficients[2,4],3)
    Eq <- paste("y= ",A,"x +(",B,")",sep="")
    
    xx <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=length(x))
    x1Predict <- xx
    Predict <- predict(Modelo,newdata=data.frame(x1=x1Predict))
    lines(Predict~x,col="red",lwd=2)

    legend("bottomright",c(Eq,substitute(paste(R^{2},"= ",C,sep=""),list(C=C)),expression(italic(p)*"-value= 0.402")),cex=0.75,horiz=FALSE,bty="n")
    mtext(expression(paste(Delta,"E", sep="")),2,line=2.4,cex=0.75,font=1)
    mtext(expression(paste(Delta,"Ta", sep="")),1,line=2.4,cex=0.75,font=1)


    ## AB0~ATa
    Is.NA <- which(is.na(AB0) | is.na(ATa))
    x <- ATa[-Is.NA]
    y <- AB0[-Is.NA]
    plot(y~x,ylab="",xlab="",cex.axis=0.75,pch=21,cex=.4)
     Modelo <- lm(y~x)
    A <-round(summary(Modelo)$coefficients[2,1],3)
    B <-round(summary(Modelo)$coefficients[1,1],3)
    C <-round(summary(Modelo)$r.squared,3)
    D <- round(summary(Modelo)$coefficients[2,4],3)
    Eq <- paste("y= ",A,"x +(",B,")",sep="")
    
    xx <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=length(x))
    x1Predict <- xx
    Predict <- predict(Modelo,newdata=data.frame(x1=x1Predict))
    lines(Predict~x,col="red",lwd=2)

    legend("bottomright",c(Eq,substitute(paste(R^{2},"= ",C,sep=""),list(C=C)),expression(italic(p)*"-value= 0.016")),cex=0.75,horiz=FALSE,bty="n")
    mtext(expression(paste(Delta,"B0", sep="")),2,line=2.4,cex=0.75,font=1)
    mtext(expression(paste(Delta,"Ta", sep="")),1,line=2.4,cex=0.75,font=1)

    ## APpk~ATa
    Is.NA <- which(is.na(APpk) | is.na(ATa))
    x <- ATa[-Is.NA]
    y <- APpk[-Is.NA]
    plot(y~x,ylab="",xlab="",cex.axis=0.75,pch=21,cex=.4)
    Modelo <- lm(y~x)
    A <-round(summary(Modelo)$coefficients[2,1],3)
    B <-round(summary(Modelo)$coefficients[1,1],3)
    C <-round(summary(Modelo)$r.squared,3)
    D <- round(summary(Modelo)$coefficients[2,4],3)
    Eq <- paste("y= ",A,"x +(",B,")",sep="")
    
    xx <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=length(x))
    x1Predict <- xx
    Predict <- predict(Modelo,newdata=data.frame(x1=x1Predict))
    lines(Predict~x,col="red",lwd=2)

    legend("bottomright",c(Eq,substitute(paste(R^{2},"= ",C,sep=""),list(C=C)),expression(italic(p)*"-value= 0.471")),cex=0.75,horiz=FALSE,bty="n")
    mtext(expression(paste(Delta,"Ppk", sep="")),2,line=2.4,cex=0.75,font=1)
    mtext(expression(paste(Delta,"Ta", sep="")),1,line=2.4,cex=0.75,font=1)



    #### WE NEED AT LEAST 3 DATA POINTS BEFORE Tpk (ALEX TAKES 2), but see TotalID[which(AB0>50)]

    ### PROBLEM with TIME, what should we use as T0??? Sometimes T0 is > than T1. WHICH should be used as initial temperature? the lowest? That is 0 or the value ??

    ### We have less data than Alex, look at the diary


   ############################### PLOTS ACCLIMATION RATE########################
    par(mfrow=c(2,2))

    ## ATPK~ATA
    x <- At
    y <- ATpk/ATa
    Is.NA <- which(is.na(x) | is.na(y) )
    x <- x[-Is.NA]
    y <- y[-Is.NA]

    plot(y~x,ylab="",xlab="",cex.axis=0.75,pch=21,cex=.4)
    Modelo <- lm(y~x)
    A <-round(summary(Modelo)$coefficients[2,1],3)
    B <-round(summary(Modelo)$coefficients[1,1],3)
    C <-round(summary(Modelo)$r.squared,3)
    D <- round(summary(Modelo)$coefficients[2,4],3)
    Eq <- paste("y= ",A,"x +(",B,")",sep="")
    
    xx <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=length(x))
    x1Predict <- xx
    Predict <- predict(Modelo,newdata=data.frame(x1=x1Predict))
    lines(Predict~x,col="red",lwd=2)

    legend("bottomright",c(Eq,substitute(paste(R^{2},"= ",C,sep=""),list(C=C)),expression(italic(p)*"-value= 0.906")),cex=0.75,horiz=FALSE,bty="n")
    mtext(expression(paste(Delta,"Tpk/",Delta,"Ta", sep="")),2,line=2.4,cex=0.75,font=1)
    mtext(expression(paste(Delta,"t", sep="")),1,line=2.4,cex=0.75,font=1)

    ## AE~ATa
    x <- At
    y <- AE/ATa
    Is.NA <- which(is.na(x) | is.na(y) )
    x <- x[-Is.NA]
    y <- y[-Is.NA]

    plot(y~x,ylab="",xlab="",cex.axis=0.75,pch=21,cex=.4)
    Modelo <- lm(y~x)
    A <-round(summary(Modelo)$coefficients[2,1],3)
    B <-round(summary(Modelo)$coefficients[1,1],3)
    C <-round(summary(Modelo)$r.squared,3)
    D <- round(summary(Modelo)$coefficients[2,4],3)
    Eq <- paste("y= ",A,"x +(",B,")",sep="")
    
    xx <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=length(x))
    x1Predict <- xx
    Predict <- predict(Modelo,newdata=data.frame(x1=x1Predict))
    lines(Predict~x,col="red",lwd=2)

    legend("bottomright",c(Eq,substitute(paste(R^{2},"= ",C,sep=""),list(C=C)),expression(italic(p)*"-value= 0.869")),cex=0.75,horiz=FALSE,bty="n")
    mtext(expression(paste(Delta,"E/",Delta,"Ta", sep="")),2,line=2.4,cex=0.75,font=1)
    mtext(expression(paste(Delta,"t", sep="")),1,line=2.4,cex=0.75,font=1)


    ## AB0~ATa
    x <- At
    y <- AB0/ATa
    Is.NA <- which(is.na(x) | is.na(y) )
    x <- x[-Is.NA]
    y <- y[-Is.NA]
    
    plot(y~x,ylab="",xlab="",cex.axis=0.75,pch=21,cex=.4)
     Modelo <- lm(y~x)
    A <-round(summary(Modelo)$coefficients[2,1],3)
    B <-round(summary(Modelo)$coefficients[1,1],3)
    C <-round(summary(Modelo)$r.squared,3)
    D <- round(summary(Modelo)$coefficients[2,4],3)
    Eq <- paste("y= ",A,"x +(",B,")",sep="")
    
    xx <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=length(x))
    x1Predict <- xx
    Predict <- predict(Modelo,newdata=data.frame(x1=x1Predict))
    lines(Predict~x,col="red",lwd=2)

    legend("bottomright",c(Eq,substitute(paste(R^{2},"= ",C,sep=""),list(C=C)),expression(italic(p)*"-value= 0.016")),cex=0.75,horiz=FALSE,bty="n")
    mtext(expression(paste(Delta,"B0/",Delta,"Ta", sep="")),2,line=2.4,cex=0.75,font=1)
    mtext(expression(paste(Delta,"t", sep="")),1,line=2.4,cex=0.75,font=1)

    ## APpk~ATa
    x <- At
    y <- APpk/ATa
    Is.NA <- which(is.na(x) | is.na(y) )
    x <- x[-Is.NA]
    y <- y[-Is.NA]
    
    plot(y~x,ylab="",xlab="",cex.axis=0.75,pch=21,cex=.4)
    Modelo <- lm(y~x)
    A <-round(summary(Modelo)$coefficients[2,1],3)
    B <-round(summary(Modelo)$coefficients[1,1],3)
    C <-round(summary(Modelo)$r.squared,3)
    D <- round(summary(Modelo)$coefficients[2,4],3)
    Eq <- paste("y= ",A,"x +(",B,")",sep="")
    
    xx <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=length(x))
    x1Predict <- xx
    Predict <- predict(Modelo,newdata=data.frame(x1=x1Predict))
    lines(Predict~x,col="red",lwd=2)

    legend("bottomright",c(Eq,substitute(paste(R^{2},"= ",C,sep=""),list(C=C)),expression(italic(p)*"-value= 0.471")),cex=0.75,horiz=FALSE,bty="n")
    mtext(expression(paste(Delta,"Ppk/",Delta,"Ta", sep="")),2,line=2.4,cex=0.75,font=1)
    mtext(expression(paste(Delta,"t", sep="")),1,line=2.4,cex=0.75,font=1)

