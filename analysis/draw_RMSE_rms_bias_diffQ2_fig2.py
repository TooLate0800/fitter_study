#!/usr/bin/env python3

import argparse
import pickle

import numpy
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

from numpy import loadtxt, zeros, linspace
from pylab import plot, show, xlabel, ylabel, title
from math import sqrt
from matplotlib.pyplot import fill_between, xlim, ylim, xlabel, ylabel, subplots_adjust, xticks, yticks, legend, xscale, errorbar, subplots, locator_params, text, figure, savefig
from scipy.interpolate import interp1d
import matplotlib as mpl
import matplotlib.font_manager as font_manager
#import matplotlib.pyplot as plt
plt.rcParams["font.weight"] = "bold"

#fig = plt.Figure()
#ax = plt.gca()
#ax.minorticks_on()
#plt.grid(True)
#===========six Deuteron models + four theoretical =============
x = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
#bias = [0.00246, 0.00077, 0.00111, 0.00080, 0.00123, 0.00095, 0.00194, 0.00395, 0.00549, 0.01550, 0.01177]
#rms = [0.00984, 0.00510, 0.00594, 0.00441, 0.00275, 0.00251, 0.00199, 0.00170, 0.00151, 0.01461, 0.00507]
RMSE = [0.01019, 0.00524, 0.00610, 0.00456, 0.00316, 0.00281, 0.00310, 0.00463, 0.00597, 0.02151, 0.01323]
shift = [0.005]*11
data_1 = {#rms
'x' : [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
'y' : [0.01019, 0.00524, 0.00610, 0.00456, 0.00316, 0.00281, 0.00310, 0.00463, 0.00597, 0.02151, 0.01323],
'yerr' : [0.00984, 0.00510, 0.00594, 0.00441, 0.00275, 0.00251, 0.00199, 0.00170, 0.00151, 0.01461, 0.00507]}
data_2 = {#bias
'x' : numpy.add([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],shift),
'y' : [0.01019, 0.00524, 0.00610, 0.00456, 0.00316, 0.00281, 0.00310, 0.00463, 0.00597, 0.02151, 0.01323],
'yerr' :[0.00246, 0.00077, 0.00111, 0.00080, 0.00123, 0.00095, 0.00194, 0.00395, 0.00549, 0.01550, 0.01177]}

plt.figure(figsize=(10.0,6.5))
plt.xlim((0.0,1.01))
#plt.ylim((0.0001,0.01))
plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16, color='black',weight='bold')
plt.yticks([0.00, 0.01, 0.02, 0.03, 0.04], fontsize=16, color='black',weight='bold')

#plt.figure()
#for data in [data_1,data_2]:
plt.errorbar(**data_1, alpha=.75, fmt=':', capsize=3, capthick=1, lw = 2, label='Variance $\sigma_{\mathrm{total}}$')
data_1 = {
    'x': data_1['x'],
    'y1': [y - e for y, e in zip(data_1['y'], data_1['yerr'])],
    'y2': [y + e for y, e in zip(data_1['y'], data_1['yerr'])]}
plt.fill_between(**data_1, alpha=.25)
plt.errorbar(**data_2, alpha=.75, fmt=':', capsize=3, capthick=1, lw = 2, label='Bias $\delta r_\mathrm{p}$')
data_2 = {
    'x': data_2['x'],
    'y1': [y - e for y, e in zip(data_2['y'], data_2['yerr'])],
    'y2': [y + e for y, e in zip(data_2['y'], data_2['yerr'])]}
plt.fill_between(**data_2, alpha=.25)
plt.plot(x,RMSE,'o',markersize=8, color='black', label='RMSE')
plt.legend(fontsize=18)

plt.xlabel(r'Q$^2_\mathrm{max}$ [(GeV/c)$^2$]', fontsize=20, weight='bold')
plt.ylabel('RMSE [fm]', fontsize=20, weight='bold')

plt.yscale('log')
plt.savefig("RMSEgraph_fig2.png")
plt.savefig("Mainz_fitter_fig2.pdf")
#plt.show()

