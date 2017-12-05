"""
Created on Tue Mar 10 10:46:11 2015

@author: sam
"""
import ThermalResponse_model_fitting
import subprocess
import rpy2.robjects as robjects
r=robjects.r
r.setwd("/home/sam/Documents/CMEECoursework/CMEELongPrac/Code")
r.source("NLLS_data_preparation.R")
print("Starting data preparation ...")
r.prep()
print("Completed data preparation")
print("Starting plotting ...")
r.plot()
print("Completed plotting")
print(" Starting initial parameter values generation ...")
r.param()
print("Finished data preparation")
print("Starting model fitting")
ThermalResponse_model_fitting.main("../Data/ThermResp_startvals.csv")
r.source("Thermal_plotting.R")
r.run_all()
print("Plotting models for datasets. \n This may take some time...")
r.thermal_plots()
print("Completed plotting models for datasets")
print("Cropping figures...")
subprocess.call(["bash", "cropall.sh", "../Results/Figures/"])
print("Completed plotting figures")
print("Starting LaTeX pdf generation")
subprocess.call(["bash", "LatexCompile.sh", "../Code/","Report","../Report/"])
print("Completed LaTeX pdf generation")