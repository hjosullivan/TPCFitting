#!/usr/bin/env python
""" This code performs non-linear least squares fitting of different
unimodal functions to experimental thermal response curves."""
__author__ = 'Samuel Thompson (samuel.thompson14@imperial.ac.uk)'
__version__ = '0.0.1'
# This code performs non-linear least squares fitting of different
# unimodal functions to experimental thermal response curves.
# Scroll down to the section called "MAIN CODE" first.

from math import log, exp, pi
import scipy
from lmfit import Minimizer, minimize, Parameters, Parameter, report_fit, fit_report
import csv
import numpy
from numpy import *
import re
import sys



#############################
# F  U  N  C  T  I  O  N  S #
#############################

def GauGo_eq(temps, E, E_D, T_pk, theta, max_trait):
	"""Gaussian-Gompertz model, used for calculating trait values at 
	a given temperature"""
	function = max_trait * exp((-E * (temps - T_pk) * (temps - T_pk) - exp((E_D * (temps - T_pk) - theta))))                 
	return numpy.array(function, dtype=numpy.float64) 

def schoolf_eq(temps, B0, E, E_D, T_pk):
	"""Schoolfield model, used for calculating trait values at a given temperature"""
	function = B0 * exp(-E * ((1/(K*temps)) - (1/(K*283.15)))) / (1 + (E/(E_D - E)) * exp((E_D / K) * (1 / T_pk - 1 / temps)))
	return numpy.array(map(log,function), dtype=numpy.float64)

def GauGo(params, temps, traits):
	"""Gaussian - Gompertz model to be called by the gaugo_model() 
	function"""
	parvals = params.valuesdict()
	E = parvals['E']
	E_D = parvals['E_D']
	T_pk = parvals['T_pk']
	theta = parvals['theta']
	max_trait = parvals['max_trait']
	traits_pred = GauGo_eq(temps, E, E_D, T_pk, theta, max_trait)
	traits_diff = traits_pred - traits
	return(traits_diff)

def schoolf(params, temps, traits):
	"""Schoolfield model, to be called by schoolfield_model()"""
	parvals = params.valuesdict()
	B0 = parvals['B0_start']
	E = parvals['E']
	E_D = parvals['E_D']
	T_pk = parvals['T_pk']
	traits_pred = schoolf_eq(temps, B0, E, E_D, T_pk)
	traits_diff = (e**traits_pred) - traits+0.00001
	return(traits_diff)
	#~ 
def gaugo_model(datatest):
	"""NLLS fitting to the Gaussian-Gompertz model; this 
	function will contain the lmfit.minimize calls to the GauGo() 
	function. It will return the results of the fitting"""
	T_pk_start = float(datatest[0][6]) - 273.15
	E_start = float(datatest[0][8])
	E_D_start = float(datatest[0][9])
	traits = numpy.array(datatest[0:,4], dtype = float)+0.00001
	temps = numpy.array(datatest[0:,5], dtype = float)
	max_trait = max(traits)
#	max_trait2 = max_trait
	params = Parameters()
	params.add('theta', vary = True, value=7)
	params.add('E', value=E_start, vary= True, min=-1000, max=1000)
	params.add('E_D', value=E_D_start, vary = True, min=-1000, max=1000)
	params.add('T_pk', value=T_pk_start, vary = True, min=0.0001, max=70)
	params.add('max_trait', value=max_trait, vary = False)#,min=max_trait2,max=1000)
	out = minimize(GauGo, params, args=(temps, traits),method="nelder")
	par_out_gaugo = out.params
	## Calculates the r squared value to determine how well the model fits the data.
	r_squared_gaugo = 1-out.residual.var()/numpy.var(traits)
	nvarys_gaugo= out.nvarys
	ndata_gaugo = out.ndata
	return(par_out_gaugo, r_squared_gaugo,nvarys_gaugo, ndata_gaugo,out.chisqr)

def schoolfield_model(datatest):
	"""NLLS fitting to the Schoolfield model; this function will 
	contain the lmfit.minimize calls to the schoolf() function. This is 
	where you can constrain the parameters."""
#	 Prepare the parameters and their bounds:
#	 Importing the data
	T_pk_start = float(datatest[0][6])
	B0_start = float(datatest[0][7])
	E_start = float(datatest[0][8])
	E_D_start = float(datatest[0][9])
#	 The dataset which contains temps and trait values.
	traits = numpy.array(datatest[0:,4], dtype = float)
	temps = numpy.array(datatest[0:,5], dtype = float) + 273.15
#	 Defining parameters
	params = Parameters()
	params.add('theta', vary = False, value=7)
	params.add('B0_start', value = B0_start, vary = True, min = -10, max = 1000)
	params.add('E', value=E_start, vary= True, min=0.0000000000000001, max=10)
	params.add('E_D', value=E_D_start, vary = True, min=0.0000000000000001, max=10)
	params.add('T_pk', value=T_pk_start, vary = True, min=270, max=330)
#	 Minimising the Model
	out = minimize(schoolf, params, args=(temps, traits),method="leastsq")
	par_out_school = out.params
#	 Calculates the r squared value to determine how well the model fits the data.
	r_squared_school = 1-out.residual.var()/numpy.var(traits)
	nvarys_school= out.nvarys
	ndata_school = out.ndata
	return(par_out_school, r_squared_school, nvarys_school, ndata_school,out.chisqr)
	
### Cubic function to compare against

def cubic_eq(temps, alpha, beta, gamma, epsilon):
	"""Cubic equation, used for calculating trait values at 
	a given temperature"""
	function = alpha + beta * temps + gamma * (temps ** 2) + epsilon * (temps **3)
	return(numpy.array(function,dtype=numpy.float64))	
	
def cubic(params, temps, traits):
	"""Cubic model, to be called by cubic_model()"""
	parvals = params.valuesdict()
	alpha = parvals['alpha']
	beta = parvals['beta']
	gamma = parvals['gamma']
	epsilon = parvals['epsilon']
	traits_pred = cubic_eq(temps, alpha, beta, gamma, epsilon)
	traits_diff = traits_pred - traits
	return(traits_diff)	
	
	
def cubic_model(datatest):
	"""NLLS fitting to the cubic model; this function will 
	contain the lmfit.minimize calls to the cubicf() function. This is 
	where you can constrain the parameters."""
#	 Prepare the parameters and their bounds:
#	 Importing the data
#	 The dataset which contains temps and trait values.
	traits = numpy.array(datatest[0:,4], dtype = float)
	temps = numpy.array(datatest[0:,5], dtype = float)
#	 Defining parameters
	params = Parameters()
	params.add('alpha', vary = True, value=1, min = -1000, max = 1000)
	params.add('beta', vary = True, value=0, min = -1000, max = 1000)
	params.add('gamma', vary = True, value=0, min = -1000, max = 1000)
	params.add('epsilon', vary = True, value=0, min = -1000, max = 1000)
#	 Minimising the Model
	out = minimize(cubic, params, args=(temps, traits),method="leastsq")
	par_out_school = out.params
#	 Calculates the r squared value to determine how well the model fits the data.
	r_squared_school = 1-out.residual.var()/numpy.var(traits)
	nvarys_school= out.nvarys
	ndata_school = out.ndata
	return(par_out_school, r_squared_school, nvarys_school, ndata_school,out.chisqr)



def AICrss(n, k, rss):
	"""Calculate the Akaike Information Criterion value, using:

	- n:   number of observations
	- k:   number of parameters
	- rss: residual sum of squares
	"""
	return n * log((2 * pi) / n) + n + 2 + n * log(rss) + 2 * k

def BICrss(n, k, rss):
	"""Calculate the Bayesian Information Criterion value, using:
	
	- n:   number of observations
	- k:   number of parameters
	- rss: residual sum of squares
	"""

	return n + n * log(2 * pi) + n * log(rss / n) + (log(n)) * (k + 1)

### Function to find the best model, given the AIC and BIC values.
### The AIC value will be preferred over the BIC value if they contradict.

def best_mod(AIC_school = None,BIC_school = None, AIC_gaugo = None, BIC_gaugo = None, model_success = [False, False]):
	"""Decides on the best model for the dataset"""
	if(model_success[1]==False):
		if(model_success[0]==True):
			return("school")
		else:
			return(None)
	else:
		if(model_success[0]==True):
			if(AIC_school < AIC_gaugo):
				return("school")
			else:
				if(AIC_school == AIC_gaugo):
					if(BIC_school < BIC_gaugo):
						return("school")
					else:
						return("gaugo")
				else:
					return("gaugo")
		else:
			return("gaugo")

def success_mod(rss,rsg,rsc,threshold):
	"""Defines the different models as successful or not dependent upon a 
	specific threshold residual sum of squares value."""
	ret = [True,True,True]
	if(rsg < threshold):
		ret[0]=False
	if(rss < threshold):
		ret[1]=False
	if(rsc < threshold):
		ret[2]=False
	return(ret)

#~ ############################
#~ # M  A  I  N    C  O  D  E #
#~ ############################
def main(argv):
	"""Performs fitting to the Gaussian-Gompertz, Schoolfield and Cubic model,
	and returns the best fits as a csv file to ../Results/results.csv"""
	#Produce an error is there is no dataset provided.
	data = numpy.genfromtxt(argv,dtype = None,delimiter = ',',deletechars='"')
	#input file "../Data/ThermResp_startvals.csv"
	# Define the Boltzmann constant (units of eV * K^-1).
	global K
	K = 8.617 * 10 ** (-5)
	#sampdata = numpy.genfromtxt("../Data/samp.csv",dtype = None,delimiter = ',')
	#Open the csv file to write the output to.
	ids = list(set(data[1:,1]))
	results = open("../Results/results.csv", 'w')
	results_csv = csv.writer(results, delimiter=",")
	results_csv.writerow(
				['Unique_id','Species_stand', 'Reference', 
				'Latitude', 'Longitude', 'Trait', 'Trait_Vals', 'Temp_Vals',
				'E_gaugo', 'E_stderr_gaugo', 'T_pk_gaugo', 'T_pk_stderr_gaugo', 'E_D_gaugo', 
				'E_D_stderr_gaugo', 'theta', 'theta_stderr', 'R_Squared_gaugo', 'B0_school', 
				'B0_stderr_school', 'E_school', 'E_stderr_school', 'T_pk_school', 'T_pk_stderr_school', 'E_D_school', 
				'E_D_stderr_school', 'R_Squared_school', 
				'Model_name', 'DataPoints_rise', 'DataPoints_fall', 
				'AIC_GauGo', 'AIC_Schoolf', 'BIC_GauGo','BIC_Schoolf', "AIC_cubic",
				'School_success','Gaugo_success',"Cubic_success","School_qual","Gaug_qual","Cubic_qual",
				"R_Squared_cubic","alpha","beta","gamma","epsilon","Choose_cubic"])
	num = 0
	sc = 0
	gg = 0
	for i in ids:
		res = Parameters()
		res.add('theta', value=None)
		res.add('B0_start', value=None)
		res.add('E', value=None)
		res.add('E_D', value=None)
		res.add('T_pk', value=None)
		res=(res,None,None,None)
		res2 = (res,None,None,None)
		AIC_school = None
		AIC_gaugo = None
		BIC_gaugo = None
		BIC_school = None
		x = data[data[:,1] == i]
		model_success = [False,False,False]
		try:
			res = schoolfield_model(x)
			AIC_school = AICrss(res[3],res[2],res[4])
			BIC_school = BICrss(res[3],res[2],res[4])
			model_success[0] = True
			sc = sc + 1
		except:
			print("\nCannot produce schoolfield model for data for " + i)
		try:
			res2 = gaugo_model(x)
			AIC_gaugo = AICrss(res2[3],res2[2],res2[4])
			BIC_gaugo = BICrss(res2[3],res2[2],res2[4])
			model_success[1] = True
			gg = gg + 1
		except:
			print("\nCannot produce gaugo model for data for " + i)
		try:
			res3 = cubic_model(x)
			AIC_cubic = AICrss(res3[3],res3[2],res3[4])
			model_success[2] = True
		except:
			print("\nCannot produce cubic model for data for " + i)
		model_choice =  best_mod(AIC_school,BIC_school,AIC_gaugo,BIC_school, model_success)
		l = 0
		suc = success_mod(res[1],res2[1],res3[1],0.75)
		choose_cubic=False	
		if(AIC_cubic<AIC_gaugo):
			if(AIC_cubic<AIC_school):
				choose_cubic=True
				model_choice="cubic"
		for j in x:
			params_s = res[0]
			params_g = res2[0]
			params_c = res3[0]
			n  = x[l]
			ref = n[1].split("-:-")
			ref[0] = re.sub('"',"",ref[0])
			ref[3] = re.sub('"',"",ref[3])
			results_csv.writerow([re.sub('"',"",n[1]),ref[0],ref[1],ref[2],ref[3],"growth rate",n[4],n[5],params_g['E'].value,params_g['E'].stderr, 
				params_g['T_pk'].value,params_g['T_pk'].stderr,params_g['E_D'].value,params_g['E_D'].stderr, params_g['theta'].value, 
				params_g['theta'].stderr,res2[1],params_s['B0_start'].value, params_s['B0_start'].stderr, 
				params_s['E'].value, params_s['E'].stderr, params_s['T_pk'].value,params_s['T_pk'].stderr, params_s['E_D'].value, 
				params_s['E_D'].stderr, res[1], model_choice, 'DataPoints_rise','DataPoints_fall',
				AIC_gaugo, AIC_school, BIC_gaugo, BIC_school,AIC_cubic,
				model_success[0], model_success[1],model_success[2],suc[1],suc[0],suc[2],
				res3[1],params_c["alpha"].value,params_c["beta"].value,params_c["gamma"].value,params_c["epsilon"].value,choose_cubic])
			l = l + 1			
		num = num + 1
		sys.stdout.write("\r[" + "=" * (num / 20) +  " " * ((len(ids) - num)/ 20) + "]" +  str(num * 100 / len(ids)) + "%" + " Completed model for " + str(num) + "      ")
		sys.stdout.flush()
	results.close()	
	print("\nNumber of Schoolfield Models is:")
	print(sc)
	print("Number of Gaugo Models is:")
	print(gg)
				
if __name__ == "__main__":
	main(sys.argv[1])		

