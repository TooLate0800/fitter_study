

import argparse
#import pickle
#import scipy
#from scipy.stats import norm
import numpy
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
#==========load txt file

x = []
rms = []
RMSE = []
y = [1, 2, 3, 4, 5, 6, 7, 8, 9]
R_input=[0.8300, 0.8300, 0.8300, 0.8628, 0.8682, 0.8965, 0.8500, 0.8868, 0.8792]
j = 1
M = 8
while j< M:
    l = str(j)
    #data1 = numpy.loadtxt('../fit_PRad/fit_result/CF4_Model_' + l +'_1e4.txt', float)#load the fitted radius
    #data1 = numpy.loadtxt('../fit/loop_9models_roofit/roofit_result_PRad_range/R12_Model' + l +'_1e4.txt', float)#load the fitted radius
    #data1 = numpy.loadtxt('../fit/loop_9models/fit_result_fullRange_20210219/CF3_Model' + l +'_1e4.txt', float)#load the fitted radius
    #data1 = numpy.loadtxt('../PRadII_rebin/fit/fit_result/R11_PRadII_Model' + l +'_1e4.txt', float)
    data1 = numpy.loadtxt('../global_fit_new/result/poly10_xsfit_RE_GMfromRatio_check_Model' + l +'.txt', float)
    #data1 = numpy.loadtxt('../global_fit_new/result/R11_xsfit_RE_GMRandom_Model' + l +'.txt', float)
    #radius = data1[:]
    r = data1[:]
    i = 0
    radius = []
    while (i<len(r)):
        if (r[i]>-1.0 and r[i]<1.0):
        #if (r[i]>0.72 and r[i]<0.98):
            radius.append(r[i])
        i += 1
    R_mean = numpy.mean(radius)
    R_rms = numpy.std(radius)
    x.append(R_mean-R_input[j-1])
    rms.append(R_rms)
    RMSE.append(numpy.sqrt((R_mean-R_input[j-1])**2+R_rms**2)) 
    j += 1
print(x)
print(rms)
print(RMSE)
print(numpy.mean(numpy.absolute(x)))
print(numpy.mean(rms))
print(numpy.mean(RMSE))
