# ALIST-library of spectral templates.

Using the SDSS-IV/APOGEE stars from their latest DR16 data, we generate a new empirical stellar library to generate a new high-resolution, NIR SSP spectral model library. We genrate this library using 2 different isochrones: Padova-based and MIST-based. 

The files conatining the SSP templates are:

Filename: ALIST_padova_may2020.fits
Data Model:
HDU1 = Table containing the age, metallicity and alpha abundance and luminosity of the models
HDU2 = Array (NX8575) of all Padova-based models.

(To be uploaded):
Filename: ALIST_MIST_may2020.fits 
Data Model:
HDU1 = Table containing the age, metallicity and alpha abundance and luminosity of the models
HDU2 = Array (NX8575) of all Padova-based models.

An example code to read in the models using Python:

sspgrid = Table.read('ALIST_padova_may2020.fits')
ssp_spec = fits.open('ALIST_padova_may2020.fits')[2].data
age, feh, afe = 10, 0.0, 0.0
model_id = np.where((sspgrid['AGE'] == age) & (sspgrid['FE_H'] == feh) & (sspgrid['A_FE'] == afe))[0]
model_spec = ssp_spec[model_id][0]

Complete code to get the model is given in the file : 'alist_call.py'
