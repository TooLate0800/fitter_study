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
#bias = [0.0006, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0008, 0.0009, 0.0009, 0.0010, 0.0011, 0.0011, 0.0012, 0.0013, 0.0021, 0.0039, 0.0063]
#rms = [0.0205, 0.0095, 0.0061, 0.0040, 0.0031, 0.0023, 0.0018, 0.0016, 0.0013, 0.0012, 0.0010, 0.0009, 0.0009, 0.0008, 0.0006, 0.0004, 0.0003]
x = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.3, 0.4]
RMSE = [0.0206, 0.0096, 0.0062, 0.0041, 0.0033, 0.0026, 0.0021, 0.0019, 0.0018, 0.0017, 0.0016, 0.0016, 0.0016, 0.0017, 0.0023, 0.0040, 0.0063]
shift = [0.002]*17
data_1 = {#rms
'x' : [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.3, 0.4],
'y' : [0.0206, 0.0096, 0.0062, 0.0041, 0.0033, 0.0026, 0.0021, 0.0019, 0.0018, 0.0017, 0.0016, 0.0016, 0.0016, 0.0017, 0.0023, 0.0040, 0.0063],
'yerr' : [0.0205, 0.0095, 0.0061, 0.0040, 0.0031, 0.0023, 0.0018, 0.0016, 0.0013, 0.0012, 0.0010, 0.0009, 0.0009, 0.0008, 0.0006, 0.0004, 0.0003]}
data_2 = {#bias
'x' : numpy.add([0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.3, 0.4],shift),
'y' : [0.0206, 0.0096, 0.0062, 0.0041, 0.0033, 0.0026, 0.0021, 0.0019, 0.0018, 0.0017, 0.0016, 0.0016, 0.0016, 0.0017, 0.0023, 0.0040, 0.0063],
'yerr' : [0.0006, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0008, 0.0009, 0.0009, 0.0010, 0.0011, 0.0011, 0.0012, 0.0013, 0.0021, 0.0039, 0.0063]}

plt.figure(figsize=(10.0,6.5))
plt.xlim((0.0,0.41))
#plt.ylim((0.0001,0.01))
plt.xticks([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40], ['0.00', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30', '0.35', '0.40'], fontsize=16, color='black',weight='bold')
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
plt.plot(x,RMSE,'o', markersize=8, color='black', label='RMSE')
plt.legend(fontsize=18)

plt.xlabel(r'Q$^{2}_\mathrm{max}$ [(GeV/c)$^2$]', fontsize=20, weight='bold')
plt.ylabel('RMSE [fm]', fontsize=20, weight='bold')

plt.yscale('log')
plt.savefig("RMSEgraph.png")
plt.savefig("Mainz_fitter_fig1.pdf")
#plt.show()

