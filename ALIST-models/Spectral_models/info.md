In this folder, you can find the A-LIST spectral models based on Paodva isochrones and MIST isochrones as 'ALIST_Padova.fits' and 'ALIST_MIST.fits' respectively.

Each files is in the .fits format which contains the following data:

HDU1: A table containing the Age, [M/H], [alpha/M], fractional luminosity, delta temperature, calculated [M/H] and calculated [alpha/M]
(For more information on last 4 columns, refer Section 4.1 in the paper)

HDU2: A-LIST spectral models. This is a (nx8575) array where n is the length of the table above. All models that have zero flux imply those models do not exist in out library (corresponding NaN's in the table)

As mentioned in our paper (Section 5.4), we suggest using a quality cut based on the fractional luminsoity and delta temperature to select the models.
A sample code to read in these spectra is given in this repository ([alist_call.py](alist_call.py))


