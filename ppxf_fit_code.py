from __future__ import print_function
try:
    import pyfits
except:
    from astropy.io import fits as pyfits
from scipy import ndimage, signal
import numpy as np
from time import clock
import glob
import matplotlib.pyplot as plt
from ppxf.ppxf import ppxf
import ppxf.ppxf_util as util
from pysynphot import observation
from pysynphot import spectrum
np.set_printoptions(threshold=100000)
from astropy.table import Table
from astropy.stats import sigma_clip
import scipy as sp
from scipy import signal
from PyAstronomy import pyasl
from PyAstronomy import funcFit as fuf
import pandas as pd
from collections import defaultdict
import csv
import pdb
import os
import matplotlib
import scipy
from scipy import ndimage
from numpy.ma import masked_array
import math

##Path to directories:
#gcdir = '/uufs/astro.utah.edu/common/uuastro/astro_data/zasowski/research/lg_pops/data/m31_gcs/'
#gcdir = '/uufs/astro.utah.edu/common/uuastro/astro_data/zasowski/research/lg_pops/m31_gcs/'##Sakari's Gc's
sspdir = '/uufs/astro.utah.edu/common/home/u1085975/LG_POPS/ALIST_lib/'
stellib = '/uufs/astro.utah.edu/common/home/u6018717/OliviaLGPops/pPXF/'
Miles_lib = '/uufs/astro.utah.edu/common/home/u1085975/LG_POPS/MILES_library_v1/'
#------------------------------------------------------------------------------------------
##Getting the MILES spectrum ready to use for comparing:
a = pyfits.open(Miles_lib + 'age10mh0.0.fits')[0].data
b = np.arange(1680,50000,0.9)
ixs = np.where((b >= 15100.8015) & (b <= 16999.8074))[0]
specin = a[ixs]
wavein = b[ixs]
noisein = specin*0+10
#############################################################
##Getting the GC spectrum:
##Sakari GC's
#gc = Table.read(ogdir+'apred_unshifted_spectrum_AP00422107+4132142.dat',format='ascii', names=('wave','flux','error'))
##Other APOGEE GCs:
# gc = Table.read(ogdir+'AP00430957+4121321_unshifted.dat',format='ascii',names=('wave','flux','error', 'extra'))
# wavein = gc['wave']
# specin = gc['flux']
# noisein = gc['error']
#------------------------------------------------------------------------------------------
##Defining the new wavelength range:
wavnew = np.arange(np.min(wavein),np.max(wavein),0.17)
##Reading in the model spectra
coarse_grid = []
coarse_grid.extend(glob.glob(sspdir + 'age*mh0.0al0.0.fits'))
##Defining wavelength range to fit:
wavstart = 1.585e4
wavend = 1.644e4
##########################################
c = 299792.458 # speed of light (km/s)
FWHM_gal = 60*2.355 # FWHM in Angstroms
FWHM_gc = 1.6e4/22500 # FWHM in Angstromsfor GC's
#------------------------------------------------------------------------------------------
##Rebinning data according to new wavelength array
galaxy, noise, gal_wave, velscale, log_wave = rebin_data(wavein, specin, noisein, wavnew, wavstart, wavend)
print('Data velscale = %d'%velscale)
##Masking based on bad pixels and skylines: (Galen's code)
skylines = np.loadtxt('/uufs/astro.utah.edu/common/home/u1085975/LG_POPS/input_files/rousselot2000.dat') # table of sky lines
goodPixels = skyline_masks(skylines, wavstart, wavend, gal_wave)
#------------------------------------------------------------------------------------------
##Olivia's stellar template (Galen's code)
# Read in the Stellar Library
star = glob.glob(stellib + 'StellarLib/aspcapStar*.fits') # ASPCAP star files
star.sort()
#FWHM_tem = 1.6e4/22500 # FWHM in Angstroms
# Extract the wavelength range and logarithmically rebin one spectrum
# to the same velocity scale of the galaxy spectrum, to determine
# the size needed for the array which will contain the template spectra.
hdu1 = pyfits.open(star[0]) # open first star
starfit = hdu1[3].data # best fit flux array
starfit[starfit==0] = 1. # replace 0's in chip gaps with 1's
wavelength_start = hdu1[1].header['CRVAL1'] # header with min wavelength info
wavelength_logdisp = hdu1[1].header['CDELT1'] # header with log dispersion info
num_wavelength = hdu1[1].header['NAXIS1'] # header with wavelength number info
stel_wave = 10**(wavelength_start + np.arange(0, wavelength_logdisp*num_wavelength, wavelength_logdisp))
                                        # wavelength array in Angstroms
    
#------------------------------------------------------------------------------------------
chi2_fit,vel_fit,disp_fit,velerr_fit,disperr_fit,met_fit,alp_fit,age_fit = [], [], [], [], [], [], [], []
wt_used, age_used, met_used, alp_used, lum_used = [], [], [], [], []

alist = coarse_grid

#templates, temp_wave, velscale, sigma = rebin_template_gc(alist, velscale, stel_wave, wavnew, FWHM_gc)
templates, temp_wave, velscale = rebin_template_miles(alist, velscale, stel_wave, wavnew)



print('Model velscale = %d'%velscale)
velscaletemp = c*(temp_wave[1]-temp_wave[0])
print('Verifying velscale value = %d'%velscaletemp)

# Relation between velocity and redshift in pPXF; from NED z = -0.001001
dv = (temp_wave[0]-log_wave[0])*c # km/s

vel =  0 # Initial estimate of the galaxy velocity in km/s
zest = np.exp(vel/c) - 1

#------------------------------------------------------------------------------------------
for n in range(len(noise)):
    if math.isnan(float(noise[n])) == True:
        noise[n] = 0
    if noise[n] <= 0:
        noise[n] = abs(noise[n] + 0.0000001)
    if math.isnan(float(galaxy[n])) == True:
        galaxy[n] = 0

#------------------------------------------------------------------------------------------
if len(alist) > 1:
    for a in range(len(alist)):
    ###ppxf fit:          
        # Here the actual fit starts. The best fit is plotted on the screen.
        plt.clf()
        start = [vel, 0] # (km/s), starting guess for [V,sigma]
        bounds = [[0, 10],[0,10]] # (km/s) lower and upper limits for vel and disp
        t = clock()
        fixed = [1,0]
        fig1 = plt.figure(figsize=(10,6))
        ax = fig1.gca()
        pp = ppxf(templates[:,a], galaxy, noise*0+0.001, velscale, start, bounds,
           plot=True, moments=2, mdegree=5, lam = gal_wave, fixed=None, degree=0,
              vsyst=dv, goodpixels=goodPixels, clean=False)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Wavelength [nm]',fontsize=15)
        plt.ylabel('Relative Flux',fontsize=15)
        plt.show()

        print("Formal errors:")
        print("     dV    dsigma   dh3      dh4")
        print("".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2)))

        print('Elapsed time in PPXF: %.2f s' % (clock() - t))
        print('dv:',dv)
        print('Vel error:',pp.error[0])
        print('Sigma error:',pp.error[1])
        print('Vel:',pp.sol[0])
        print('Sigma:',pp.sol[1])
        print('Chi2:',pp.chi2)



        for i in range(len(pp.weights)):
            if pp.weights[i] != 0:
                print(pp.weights[i], alist[i])

        if pp.weights[i] != 0:        
            temp_used = alist[a]
            weight_used = [[],[],[],[]]
            lumin = []
            weight_used[0].append(pp.weights[i])
            weight_used[1].append(float(temp_used.split('ge')[1].split('mh')[0]))
            weight_used[2].append(float(temp_used.split('mh')[1].split('al')[0]))
            weight_used[3].append(float(temp_used.split('al')[1].split('.fits')[0]))
            lumin.append(pyfits.getheader(str(temp_used))[9])

        lum_used.append(lumin)
        chi2_fit.append(pp.chi2)
        vel_fit.append(pp.sol[0])
        disp_fit.append(pp.sol[1])
        met_fit.append(weight_used[2])
        alp_fit.append(weight_used[3])
        age_fit.append(weight_used[1])
else:
    plt.clf()
    start = [vel, 0] # (km/s), starting guess for [V,sigma]
    bounds = [[0, 10],[0,10]] # (km/s) lower and upper limits for vel and disp
    t = clock()
    fixed = [1,0]
    fig1 = plt.figure(figsize=(10,6))
    ax = fig1.gca()
    pp = ppxf(templates, galaxy, noise, velscale, start, bounds,
       plot=True, moments=2, mdegree=5, lam = gal_wave, fixed=None, degree=0,
          vsyst=dv, goodpixels=goodPixels, clean=False)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('Wavelength [nm]',fontsize=15)
    plt.ylabel('Relative Flux',fontsize=15)
    plt.show()

    print("Formal errors:")
    print("     dV    dsigma   dh3      dh4")
    print("".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2)))

    print('Elapsed time in PPXF: %.2f s' % (clock() - t))
    print('dv:',dv)
    print('Vel error:',pp.error[0])
    print('Sigma error:',pp.error[1])
    print('Vel:',pp.sol[0])
    print('Sigma:',pp.sol[1])
    print('Chi2:',pp.chi2)



    for i in range(len(pp.weights)):
        if pp.weights[i] != 0:
            print(pp.weights[i], alist[i])

    if pp.weights[i] != 0:        
        temp_used = alist[0]
        weight_used = [[],[],[],[]]
        lumin = []
        weight_used[0].append(pp.weights[i])
        weight_used[1].append(float(temp_used.split('ge')[1].split('mh')[0]))
        weight_used[2].append(float(temp_used.split('mh')[1].split('al')[0]))
        weight_used[3].append(float(temp_used.split('al')[1].split('.fits')[0]))
        lumin.append(pyfits.getheader(str(temp_used))[9])

    lum_used.append(lumin)
    chi2_fit.append(pp.chi2)
    vel_fit.append(pp.sol[0])
    disp_fit.append(pp.sol[1])
    met_fit.append(weight_used[2])
    alp_fit.append(weight_used[3])
    age_fit.append(weight_used[1])
    
    
 ##Getting the average of minimum chi2 based parameters:   
#----------------------------------------------------------------------------------
age = np.concatenate(age_fit).ravel()
met = np.concatenate(met_fit).ravel()
chi2_fit = np.asarray(chi2_fit)
#-------------------------------------------
likeli = np.exp(-np.asarray(chi2_fit)/2)
ind = np.argsort(1/likeli)
like_age = np.average(age[ind][0:3])
like_met = np.average(met[ind][0:3])
like_alpha = np.average(alp[ind][0:3])
like_vel = np.average(vel[ind][0:3])
like_disp = np.average(disp[ind][0:3])
