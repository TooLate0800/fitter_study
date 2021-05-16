

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
y = [1, 2, 3, 4, 5, 6, 7, 8, 9]
R_input=[0.8300, 0.8300, 0.8300, 0.8628, 0.8682, 0.8965, 0.8500, 0.8868, 0.8792]
j = 1
M = 10
while j< M:
    l = str(j)
    #data1 = numpy.loadtxt('../fit_PRad/fit_result/CF4_Model_' + l +'_1e4.txt', float)#load the fitted radius
    #data1 = numpy.loadtxt('../fit/loop_9models_roofit/roofit_result_PRad_range/R12_Model' + l +'_1e4.txt', float)#load the fitted radius
    #data1 = numpy.loadtxt('../fit/loop_9models/fit_result_fullRange_20210219/CF3_Model' + l +'_1e4.txt', float)#load the fitted radius
    #data1 = numpy.loadtxt('../PRadII_rebin/fit/fit_result/R11_PRadII_Model' + l +'_1e4.txt', float)
    #data1 = numpy.loadtxt('../global_fit_new/result/R11_xsfit_all_GMfromRatio_Model' + l +'.txt', float)
    data1 = numpy.loadtxt('../global_fit_new/result/R11_xsfit_all_GMfromRatio_Cut_01GeV_Model' + l +'.txt', float)
    #data1 = numpy.loadtxt('../global_fit/result/z2_model' + l +'_1e4.txt', float)
    #radius = data1[:]
    r = data1[:,0]#chi2
    i = 0
    radius = []
    while (i<len(r)):
        #if (r[i]>-1.0 and r[i]<1.0):
        #if (r[i]>0.72 and r[i]<0.98):
        radius.append(r[i])
        i += 1
    R_mean = numpy.mean(radius)
    R_rms = numpy.std(radius)
    x.append(R_mean-R_input[j-1])
    rms.append(R_rms)
    j += 1
print(x)
print(rms)
print(numpy.mean(numpy.absolute(x)))
#plt.figure(figsize=(3.5,5.5))
plt.figure(figsize=(3.0,5.5))
#plt.figure(figsize=(4.0,5.5))
#fig = plt.Figure()
ax = plt.gca()
ax.xaxis.label.set_size(15)
#plt.xlim((-0.1,0.1))
#plt.xlim((-0.07,0.07))
#plt.xlim((-0.05,0.05))
#plt.xlim((-0.16,0.16))
plt.ylim((0.5,9.5))
plt.grid(True)
plt.xlabel(r'$\delta R / fm$',fontsize=15)
#plt.xticks([-0.12, -0.08, -0.04, 0, 0.04, 0.08, 0.12])
#plt.xticks([-0.04, -0.02, 0, 0.02, 0.04])
#plt.xticks([-0.04, -0.02, 0, 0.02, 0.04])
#plt.xticks([-0.1, -0.05, 0, 0.05, 0.1])
#plt.xticks([-0.2, -0.1, 0, 0.1, 0.2])
#plt.xticks([-0.5, -0.25, 0, 0.25, 0.5])
plt.yticks([])

y1 = np.linspace(0.5,9.5)
x1 = y1-y1 
plt.plot(x1,y1,linestyle=(0, (5, 5)),color='black')




#plt.title('CF (3)',loc='center', fontsize=15)
#plt.title('PolynomialZ (4)',loc='center', fontsize=15)
plt.title('Rational (1,1)',loc='center', fontsize=15)
#plt.title('Gaussian',loc='center', fontsize=15)

#plt.title('Gaussian',loc='center',fontdict={'fontsize': 20, 'fontweight': 'medium'})
plt.errorbar(x[0], y[0], xerr=rms[0], fmt ='o',color='black')#, color = 'r.')
plt.errorbar(x[1], y[1], xerr=rms[1], fmt ='o',color='grey')#, color = 'b.')
plt.errorbar(x[2], y[2], xerr=rms[2], fmt ='or')#, color = 'g.')
plt.errorbar(x[3], y[3], xerr=rms[3], fmt ='o',color='salmon')#, color = 'k.')
plt.errorbar(x[4], y[4], xerr=rms[4], fmt ='o',color='blue')#, color = 'm.')
plt.errorbar(x[5], y[5], xerr=rms[5], fmt ='o',color='cornflowerblue')#, color = 'olive')
plt.errorbar(x[6], y[6], xerr=rms[6], fmt ='o',color='magenta')#, color = 'brown')
plt.errorbar(x[7], y[7], xerr=rms[7], fmt ='o',color='violet')#, color = 'teal')
plt.errorbar(x[8], y[8], xerr=rms[8], fmt ='o',color='green')#, color = 'orange')
#plt.errorbar(x[9], y[9], xerr=rms[9], fmt ='o')#, color = 'orange')


plt.savefig("mygraph.png")
