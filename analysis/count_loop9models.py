

import argparse
#import pickle
#import scipy
#from scipy.stats import norm
import numpy
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
#==========load txt file
x = []
rms = []
j = 1
M = 10
while j< M:
    l = str(j)
    data1 = numpy.loadtxt('../fit/loop_9models/fit_result/R11new_Model' + l +'_1e4.txt', float)#load the fitted radius
    #radius = data1[:,0]
    radius = data1[:]
    #print(len(radius))
    R_mean = numpy.mean(radius)
    R_rms = numpy.std(radius)
    x.append(R_mean)
    rms.append(R_rms)
    #B_max = numpy.amax(Bias)
    #print(R_mean)
    #print(R_rms)
    j += 1
print(x)
print(rms)
