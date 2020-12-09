================
Fitter study 
Author : Jingyi Zhou
Last update: 20201209
================
Work Flow:

1. data/
 1) CrossSections.dat: Mainz cross section data
 2) Carl-norm.dat: Mainz G_E data from Keith Griffioen and Carl Carlson

2. generator/
 1) What it does: obtain Q^2 and stat uncertainty of GE from the input file in /data, use a selected model to generate pseudo-data, the central values of GE are smeared by the overall normalization and the statistical uncertainties.
 2) How to run: (on the mepg computer) make -j -> ./RobustGenerator
 3) Modify the code: Lines before 203 are fixed. They are the definitions of different form factor models. The main body of the code starts from line 203. The details are in the comments of the code. 

3. fit/
 1) What it does: fit the pseudo-data generated at generator/robust_table/ and write the fitted radius to the file at fit/fit_result/
 2) How to run: source Source_this_first -> make O=filename -> ./filename (for example: filename = fit_GE_1)
 3) In this package: these are the fitting program using ROOT Minimizer
    fit_GE_1: Use Rational(1,1) to fit the one set of data with one floating parameter
    fit_GE_2: Use Rational(1,1) to fit the two sets of data with two floating parameters
    fit_GE_3: Use Rational(1,1) to fit the three sets of data with three floating parameter
    xy_GE_Q2_expan_fit : second order polynomial fit
    xy_GE_z_fit: second order polynomial z-expansion fit
    xy_GE_z3_fit: third order polynomial z-expansion fit
    xy_GE_z4_fit: fourth order polynomial z-expansion fit
    fit_3para_loop: define any fitter with three fitting parameters
    fit_4para_loop: define any fitter with four fitting parameters 
 4) Modify the code: open filename.C, the details are in the comments of fit_GE_1.C
    pay attention to the following: 
       a. the definition of the fitter, radius should be one of the fitting parameter
       b. the input and output files
       c. repetition times
 5) roofit.C: Using TGraphErrors to fit, can not handle various sets of data with different floating parameters in one fit 

4. analysis/
 1) What it does: 
    count.py: calculate the rms and mean value of the fitted radius distribution   
    draw_rms_bias.py: draw the bias-variance figures as the ones in Xuefei's paper:https://arxiv.org/pdf/1803.01629.pdf
 2) How to run: source set_py3.8.sh -> python3.8 filename.py, the figure is saved as mygraph.png in this directory, to open it, type "display mygraph.png"
