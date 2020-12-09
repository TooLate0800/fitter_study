

import argparse
import pickle
import scipy
from scipy.stats import norm
import numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#==========load txt file
data1 = numpy.loadtxt('../fit/fit_result/test.txt', float)#load the fitted radius
#radius = data1[:,0]
radius = data1[:]
R_mean = numpy.mean(radius)
R_rms = numpy.std(radius)
#B_max = numpy.amax(Bias)
print(R_mean)
print(R_rms)
