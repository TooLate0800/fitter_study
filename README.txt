================
Fitter study for the proton radius extraction
Author : Jingyi Zhou
Last update: 20210516
================
Work Flow:

1. Raw data
Directory: data/
 1) CrossSections.dat: original Mainz cross section data from Bernauer's thesis(https://www.osti.gov/etdeweb/servlets/purl/21403504)
 2) Carl-norm.dat: Mainz G_E data from Keith Griffioen and Carl Carlson, the three columns are: Q^2(fm^-2), GE, dGE (arXiv:1509.06676)
 3) Carl-norm_#GeV.dat: Mainz normalized GE data at different Q^2_max cut off
 4) The cross section data at different Q^2_max cut off is in: cs_data/. They are generated based on CrossSections.dat
    Read.C: open CrossSections.dat and then produce 31 separated files(#.txt) corresponding to 31 datasets
    Read_cut.C: open CrossSections.dat and then produce various separated files(#.txt) at different Q^2_max cut off, also produce the normalization factor table(norm_table.dat) with four columns: dataset #, # of data points, two normalization factors(if there is only one factor then print -1)  

2. Generators
a. generator/RobustGenerator_loop_cut.cxx: 
   Obtain Q^2 and total uncertainty of GE from the input file in /data, use a selected model to generate pseudo-data, the central values of GE are smeared by the overall normalization and the statistical uncertainties.
b. generator_31norms/RobustGenerator_loop_CutQ2.cxx: 
   Obtain Q^2 and total uncertainty of xs from the input file in /data, use a selected GE+GM model to generate pseudo-data, the central values of xs are smeared by the overall normalization and the statistical uncertainties.

3. Fitting
a. fit/
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
    loop_9models/*.C: similar to those in fit/, but loop over all the 9 models automatically 
 4) roofit.C: Using TGraphErrors to fit, can not handle various sets of data with different floating parameters in one fit 
b. global_fit_new/R11_fit_xs_norms_cut_loop.C: Use Rational(1,1) to fit the xs pseudo-data at different Q^2_max cut off, loop over different models

4. analysis/
   count.py: calculate the rms and mean value of the fitted radius distribution   
   draw_rms_bias.py: draw the bias-variance figures as the ones in Xuefei's paper:https://arxiv.org/pdf/1803.01629.pdf
   count_draw_loop.py: combine the above two functions and make a plot
