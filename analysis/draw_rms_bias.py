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
plt.xlabel(r'$\delta R / fm$')
#===========six Deuteron models + four theoretical =============
y = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
#x = [-0.016322381,-0.01617413148,-0.00278999391,-0.002765298366,-0.002217697,-0.001562968,0.01272,0.00005,0.00743,0.0318,0.0276,-0.00001]#R11
#rms = [0.00239945936,0.002399977687,0.002418490395,0.002397701699,0.00236337406,0.002351917322,0.00264,0.00246,0.00232,0.0022,0.0023,0.0019]
#x = [-0.015873231,-0.01570853288,-0.00218836091,-0.002149178366,-0.002170536,-0.001621175,0.00849,0.00180,0.00641,0.0283,0.0273,-0.0110]#R13
#rms = [0.002461742247,0.002464264356,0.002487791592,0.002464297325,0.002443944568,0.002426219819, 0.002456761378,0.002519708452,0.002429524409,0.0025,0.0024,0.0030]
x = [-0.01322,-0.01301,0.00018,0.00025,-0.00006,0.00038,0.00406,0.00397,0.00754,0.12433,0.02648,0.02918]
rms = [0.00099,0.00099,0.00097,0.00099,0.00113,0.00111,0.00102,0.00114,0.00110,0.00629,0.00104,0.00111]
#plt.title('Polynomial(4)',loc='center')
#plt.title('Rational(2,1)',loc='center')
#plt.title('monopole',loc='center')
#plt.title('cf(3)',loc='center')
plt.figure(figsize=(2.8,5.5))
#plt.figure(figsize=(2.5,5.5))
#plt.figure(figsize=(2.5,5.5))
plt.xlim((-0.05,0.05))
plt.ylim((0.5,12.5))
#plt.xlim((-0.02,0.02))

#print(y[5])
#plt.errorbar(x[0], y[0], xerr=rms[0], fmt ='o',color='red')#, color = 'r.')
#plt.errorbar(x[1], y[1], xerr=rms[1], fmt ='o',color='blue')#, color = 'b.')
plt.errorbar(x[0], y[0], xerr=rms[0]*5/11, fmt ='o',color='black')#, color = 'r.')
plt.errorbar(x[1], y[1], xerr=rms[1]*5/11, fmt ='o',color='grey')#, color = 'b.')
plt.errorbar(x[2], y[2], xerr=rms[2]*5/11, fmt ='or')#, color = 'g.')
plt.errorbar(x[3], y[3], xerr=rms[3]*5/11, fmt ='o',color='salmon')#, color = 'k.')
plt.errorbar(x[4], y[4], xerr=rms[4]*5/11, fmt ='o',color='blue')#, color = 'm.')
plt.errorbar(x[5], y[5], xerr=rms[5]*5/11, fmt ='o',color='cornflowerblue')#, color = 'olive')
plt.errorbar(x[6], y[6], xerr=rms[6]*5/11, fmt ='o',color='magenta')#, color = 'brown')
plt.errorbar(x[7], y[7], xerr=rms[7]*5/11, fmt ='o',color='violet')#, color = 'teal')
plt.errorbar(x[8], y[8], xerr=rms[8]*5/11, fmt ='o',color='green')#, color = 'orange')
plt.errorbar(x[9], y[9], xerr=rms[9]*5/11, fmt ='o')#, color = 'orange')
plt.errorbar(x[10], y[10], xerr=rms[10]*5/11, fmt ='o')#, color = 'orange')
plt.errorbar(x[11], y[11], xerr=rms[11]*5/11, fmt ='o')#, color = 'orange')
#plt.errorbar(x[12], y[12], xerr=rms[12]*5/11, fmt ='o')#, color = 'orange')
#plt.errorbar(x[5], y[5], xerr=rms[5], fmt ='o',color='magenta')#, color = 'olive')
#plt.errorbar(x[6], y[6], xerr=rms[6], fmt ='o',color='cornflowerblue')#, color = 'brown')
plt.savefig("mygraph.png")
#plt.show()

