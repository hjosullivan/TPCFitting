load("GlobalDataset/GlobalDataset.Rdata")
source("../../Analysis/code/TPCFitting.R")
Data <- GlobalDataset
Data <- TerData
TPCFit(Data,PLOT=FALSE,OverPLOT=FALSE,Model="Schoolfield",SchoolTpk=TRUE,rand.st=TRUE, n.rand=100)
