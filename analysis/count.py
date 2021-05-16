

#import argparse
#import pickle
#import scipy
#from scipy.stats import norm
import numpy
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
#==========load txt file
#data1 = numpy.loadtxt('../global_fit_new/result/R11_xsfit_RE_GMfromRatio_Model8.txt', float)#load the fitted radius
#data1 = numpy.loadtxt('../global_fit_new/result/R11_xsfit_RE_GMisBernauer_1GeV_1e3_Model8.txt', float)#load the fitted radius
#data1 = numpy.loadtxt('../global_fit_new/result/R11_xsfit_RE_GMisDipole_1GeV_1e3_Model8.txt', float)#load the fitted radius
#data1 = numpy.loadtxt('../global_fit_new/result/R11_xsfit_RE_GMfromRatio_Cut_03GeV_chi2cut_Model1.txt', float)#load the fitted radius
#data1 = numpy.loadtxt('../global_fit_new/result/poly10_xsfit_RE_GMfromRatio_check_Model1.txt', float)#load the fitted radius
#data1 = numpy.loadtxt('../global_fit_new/result/R11_xsfit_RE_GMRandom_Model1.txt', float)#load the fitted radius
#data1 = numpy.loadtxt('../fit/fit_result/Rfiterr_R11_bootstrap_MainzinPRad_60bins.txt', float)#load the fitted radius
data1 = numpy.loadtxt('../fit/loop_9models/fit_result/CF3_Cut02GeV_Model2_1e4.txt', float)#load the fitted radius
#data1 = numpy.loadtxt('../fit_PRad/fit_result/PRad_chi2_boostrap_20_25.txt', float)#load the fitted radius
#data1 = numpy.loadtxt('../PRadII_rebin/fit/fit_result/PRadII_rebin_R11.txt', float)#load the fitted radius
#radius = data1[:,0]
#chi2 = data1[0,:]
#r = data1[1,:]
r = data1[:]
#radius = data1[:]
#print(len(radius))
i = 0
radius = []
while (i<len(r)):
    #if (chi2[i]<2000):
    if (r[i]<1.0 and r[i]>0.79):
        radius.append(r[i])
    i += 1
R_mean = numpy.mean(radius)
R_rms = numpy.std(radius)
#B_max = numpy.amax(Bias)
print(R_mean)
print(R_rms)
plt.hist(radius,bins=90)
plt.savefig("graph.png")
