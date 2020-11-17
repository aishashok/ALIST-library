##Example ssp: pull out 10 Gyr, solar metallicity and solar alpha Padova-based spectral model

import numpy as np #(version: 1.18.1)
import matplotlib.pyplot as plt #(version: 2.1.0)
from astropy.table import Table, Column #(version: 4.0.1.post1)
from astropy.io import fits #(version: 4.0.1.post1)

#Reading the .fits file to access the age, [M/H] and [$\alpha$/M] as well as the spectral models.
#Here we are using astropy package in python.
sspgrid = Table.read('A-LIST_padova.fits') #Contains the table with age, [M/H] and [$\alpha$/M]
ssp_spec = fits.open('A-LIST_padova.fits')[2].data # contains the spectral models
ssp_spec_uncert = fits.open('A-LIST_padova.fits')[3].data #contains the uncertainty spectra for each model

##Selecting the models based on the quality cuts:
ixs_quality = np.where((sspgrid['lumfrac'] > 0.32) & (sspgrid['deltatemp'] > -200) & (sspgrid['deltatemp'] < 350))[0]

##Define the age, [M/H] and [$\alpha$/M] needed to read-in:
age, m_h, a_m = 10, 0.0, 0.0
##Pulling out a spectrum using the age, [M/H] and [$\alpha$/M] values from the table:
model_id = np.where((sspgrid[ixs_quality]['AGE'] == age) & (sspgrid[ixs_quality]['M_H'] == m_h) & \
                    (sspgrid[ixs_quality]['A_M'] == a_m))[0]
model_spec = ssp_spec[model_id][0]

##Defining the wavelength range using the information provided in the header.
wavstart = float(fits.getheader('A-LIST_padova.fits')['CRVAL1'].split()[0]) ##Starting wavelength in log scale
wavdelt = float(fits.getheader('A-LIST_padova.fits')['CDELT1'].split()[0]) ##Delta wavelength in log space
numpix = int(fits.getheader('A-LIST_padova.fits')['NWAVE'].split()[0]) ##Total number of pixels available in the spectra

wavend = wavstart + (wavdelt*numpix) ##The wavelength value for the last pixel
##Using the information avilable above, we create an array of wavelength in the linear space using numpy's linspace:
wavelength = 10**(np.linspace(wavstart, wavend, num = numpix))