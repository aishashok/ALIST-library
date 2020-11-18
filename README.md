# ALIST-library of spectral templates.
Using the SDSS-IV/APOGEE stars from the latest (incremental) DR16 data, we create a new empirical stellar library to generate a new high-resolution, NIR SSP spectral model library.

We provide these SSP spectral models based on 2 isochrones: Padova and MIST models.

More information on these models are provided in our paper: Ashok et al. 2020 (submitted)

# A-LIST models information
The Spectral models are available in :

age (Gyr) = 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12

[M/H] (dex) = -2.2, -1.9, -1.6, -1.3, -1.0, 1-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4

[alpha/M] (dex) = -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4


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
The A-LIST spectral models are available for download as ".fits" files based on both the isochrones labeled: 'ALIST_Padova.fits' for the Padova-based and 'ALIST_MIST.fits' for the MIST-based models. 

Each .fits file contains the following:

HDU1: Table containing the parameters of the models. 

HDU2: The model Spectra (shape of this array is no. of models X no.of pixels)

HDU3: The Variance spectra (same shape as the model spectra)

**The Variance spectra provided is not the uncertainty/error spectra. It gives the scatter in flux per pixel of all the stellar spectra used to create the model (see Section 4 in the paper).

A sample code to read in a spectral model is provided here: [alist_call.py](alist_call.py)
