#!/usr/bin/env python3

import argparse
import pickle

import numpy
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

fig = plt.Figure()
ax = plt.gca()
ax.minorticks_on()
plt.grid(True)
#===========six Deuteron models + four theoretical =============
#x = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
#bias = [0.00246, 0.00077, 0.00111, 0.00080, 0.00123, 0.00095, 0.00194, 0.00395, 0.00549, 0.01550, 0.01177]
#rms = [0.00984, 0.00510, 0.00594, 0.00441, 0.00275, 0.00251, 0.00199, 0.00170, 0.00151, 0.01461, 0.00507]
#RMSE = [0.01019, 0.00524, 0.00610, 0.00456, 0.00316, 0.00281, 0.00310, 0.00463, 0.00597, 0.02151, 0.01323]
#x = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.3, 0.4]
#bias = [0.0373, 0.0006, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0008, 0.0009, 0.0009, 0.0010, 0.0011, 0.0011, 0.0012, 0.0013, 0.0021, 0.0039, 0.0063]
#rms = [0.1000, 0.0205, 0.0095, 0.0061, 0.0040, 0.0031, 0.0023, 0.0018, 0.0016, 0.0013, 0.0012, 0.0010, 0.0009, 0.0009, 0.0008, 0.0006, 0.0004, 0.0003]
#RMSE = [0.1070, 0.0206, 0.0096, 0.0062, 0.0041, 0.0033, 0.0026, 0.0021, 0.0019, 0.0018, 0.0017, 0.0016, 0.0016, 0.0016, 0.0017, 0.0023, 0.0040, 0.0063]
data_1 = {
'x' : [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.3, 0.4],
'y' : [0.0206, 0.0096, 0.0062, 0.0041, 0.0033, 0.0026, 0.0021, 0.0019, 0.0018, 0.0017, 0.0016, 0.0016, 0.0016, 0.0017, 0.0023, 0.0040, 0.0063],
'yerr' : [0.0205, 0.0095, 0.0061, 0.0040, 0.0031, 0.0023, 0.0018, 0.0016, 0.0013, 0.0012, 0.0010, 0.0009, 0.0009, 0.0008, 0.0006, 0.0004, 0.0003]}
data_2 = {
'x' : [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.3, 0.4],
'y' : [0.0206, 0.0096, 0.0062, 0.0041, 0.0033, 0.0026, 0.0021, 0.0019, 0.0018, 0.0017, 0.0016, 0.0016, 0.0016, 0.0017, 0.0023, 0.0040, 0.0063],
'yerr' : [0.0006, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0008, 0.0009, 0.0009, 0.0010, 0.0011, 0.0011, 0.0012, 0.0013, 0.0021, 0.0039, 0.0063]}

#plt.title('Rational(1,1)',loc='center')
plt.figure(figsize=(10.0,5.5))
#plt.xlim((0.0,1.05))
#plt.ylim((0.0,0.04))
#plt.xlim((0.0,0.41))
#plt.ylim((0.0,0.22))
plt.xlim((0.0,0.41))
#plt.ylim((-0.001,0.042))
plt.xticks([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40], ['0.00', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30', '0.35', '0.40'], fontsize=16, color='black')
plt.yticks([0.00, 0.01, 0.02, 0.03, 0.04], fontsize=16, color='black')

#plt.errorbar(x[0], RMSE[0], yerr=rms[0], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[1], RMSE[1], yerr=rms[1], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[2], RMSE[2], yerr=rms[2], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[3], RMSE[3], yerr=rms[3], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[4], RMSE[4], yerr=rms[4], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[5], RMSE[5], yerr=rms[5], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[6], RMSE[6], yerr=rms[6], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[7], RMSE[7], yerr=rms[7], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[8], RMSE[8], yerr=rms[8], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[9], RMSE[9], yerr=rms[9], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[10], RMSE[10], yerr=rms[10], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[11], RMSE[11], yerr=rms[11], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[12], RMSE[12], yerr=rms[12], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[13], RMSE[13], yerr=rms[13], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[14], RMSE[14], yerr=rms[14], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[15], RMSE[15], yerr=rms[15], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[16], RMSE[16], yerr=rms[16], fmt ='o',color='black')#, color = 'r.')
##plt.errorbar(x[17], RMSE[17], yerr=rms[17], fmt ='o',color='black')#, color = 'r.')
#plt.errorbar(x[0], RMSE[0], yerr=bias[0], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[1], RMSE[1], yerr=bias[1], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[2], RMSE[2], yerr=bias[2], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[3], RMSE[3], yerr=bias[3], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[4], RMSE[4], yerr=bias[4], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[5], RMSE[5], yerr=bias[5], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[6], RMSE[6], yerr=bias[6], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[7], RMSE[7], yerr=bias[7], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[8], RMSE[8], yerr=bias[8], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[9], RMSE[9], yerr=bias[9], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[10], RMSE[10], yerr=bias[10], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[11], RMSE[11], yerr=bias[11], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[12], RMSE[12], yerr=bias[12], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[13], RMSE[13], yerr=bias[13], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[14], RMSE[14], yerr=bias[14], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[15], RMSE[15], yerr=bias[15], fmt ='o',color='blue')#, color = 'r.')
#plt.errorbar(x[16], RMSE[16], yerr=bias[16], fmt ='o',color='blue')#, color = 'r.')
##plt.errorbar(x[17]+0.005, RMSE[17], yerr=bias[17], fmt ='o',color='blue')#, color = 'r.')
#plt.figure()
#for data in [data_1,data_2]:
plt.errorbar(**data_1, alpha=.75, fmt=':', capsize=3, capthick=1, color='black', lw = 2)
data_1 = {
    'x': data_1['x'],
    'y1': [y - e for y, e in zip(data_1['y'], data_1['yerr'])],
    'y2': [y + e for y, e in zip(data_1['y'], data_1['yerr'])]}
plt.fill_between(**data_1, alpha=.25)
plt.errorbar(**data_2, alpha=.75, fmt=':', capsize=3, capthick=1, color='blue', lw = 2)
data_2 = {
    'x': data_2['x'],
    'y1': [y - e for y, e in zip(data_2['y'], data_2['yerr'])],
    'y2': [y + e for y, e in zip(data_2['y'], data_2['yerr'])]}
plt.fill_between(**data_2, alpha=.25)

plt.xlabel(r'Q$^2$ [(GeV/c)$^2$]', fontsize=20)
plt.ylabel('RMSE [fm]', fontsize=20)

plt.yscale('log')
plt.savefig("RMSEgraph.png")
#plt.show()

