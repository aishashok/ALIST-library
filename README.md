# ALIST-library of spectral templates.


Using the SDSS-IV/APOGEE stars from the latest (incremental) DR16 data, we create a new empirical stellar library to generate a new high-resolution, NIR SSP spectral model library.  
We provide these SSP spectral models based on 2 isochrones: Padova and MIST models.
More information on these models are provided in our paper: Ashok et al. 2020 (submitted)

# Access to A-LIST
The A-LIST spectral models are available for download as ".fits" files based on both the isochrones labeled: 'ALIST_Padova.fits' for the Padova-based and 'ALIST_MIST.fits' for the MIST-based models. Each .fits file contains the following:
HDU1: Table containing the parameters of the models. 
HDU2: The model Spectra (shape of this array is 4D matrix--Age, metallicity, alpha-abundance, no.of pixels)
HDU3: The Variance spectra (4D matrix)

# Table information
The columns available in the table are:
['AGE'] = Age of the spectral model in Gyr
['M_H'] = Metallicity([M/H]) of the spectral model (unit is dex)
['A_M'] = Alpha-abundance of the spectral model (unit is dex)
['lumfrac'] = Fractional luminosity of the spectral model (For more info, refer Section 4.1.1 in the paper).
['deltatemp'] = Difference in Effective temperature between the spectral model and the underlying SSP (For more info, refer Section 4.1.2 in the paper).
[calc_MH'] = The mean metallicity calculated based on the APOGEE stars in the SSP model (For more info, refer Section 4.1.3 in the paper).
['calc_AM'] = The mean alpha-abundance calculated based on the APOGEE stars in the SSP model (For more info, refer Section 4.1.3 in the paper)
    
# Getting Started

Here is a sample code to read in the models:

import numpy as np #(version: 1.18.1)
import matplotlib.pyplot as plt #(version: 2.1.0)
from astropy.table import Table, Column #(version: 4.0.1.post1)
from astropy.io import fits #(version: 4.0.1.post1)

#Reading the .fits file to access the age, [M/H] and [$\alpha$/M] as well as the spectral models.
#Here we are using astropy package in python.
sspgrid = Table.read('A-LIST_padova_may2020.fits') #Contains the table with age, [M/H] and [$\alpha$/M]
ssp_spec = fits.open('A-LIST_padova_may2020.fits')[2].data # contains the spectral models
ssp_spec_uncert = fits.open('A-LIST_padova_may2020.fits')[3].data #contains the uncertainty spectra for each model

##Selecting the models based on the quality cuts:
ixs_quality = np.where((sspgrid['lumfrac'] > 0.32) & (sspgrid['deltatemp'] > -200) & (sspgrid['deltatemp'] < 350))[0]

##Define the age, [M/H] and [$\alpha$/M] needed to read-in:
age, m_h, a_m = 10, 0.0, 0.0
##Pulling out a spectrum using the age, [M/H] and [$\alpha$/M] values from the table:
model_id = np.where((sspgrid[ixs_quality]['AGE'] == age) & (sspgrid[ixs_quality]['M_H'] == m_h) & \
                    (sspgrid[ixs_quality]['A_M'] == a_m))[0]
model_spec = ssp_spec[model_id][0]

##Defining the wavelength range using the information provided in the header.
wavstart = float(fits.getheader('A-LIST_padova_may2020.fits')['CRVAL1'].split()[0]) ##Starting wavelength in log scale
wavdelt = float(fits.getheader('A-LIST_padova_may2020.fits')['CDELT1'].split()[0]) ##Delta wavelength in log space
numpix = int(fits.getheader('A-LIST_padova_may2020.fits')['NWAVE'].split()[0]) ##Total number of pixels available in the spectra

wavend = wavstart + (wavdelt*numpix) ##The wavelength value for the last pixel
##Using the information avilable above, we create an array of wavelength in the linear space using numpy's linspace:
wavelength = 10**(np.linspace(wavstart, wavend, num = numpix))

Follow this sample code in 

