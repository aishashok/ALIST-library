##Example ssp: pull out 10 Gyr, solar metallicity and solar alpha


from astropy.io import fits
from astropy.table import Table, columns
import numpy as np
import matplotlib.pyplot as plt


sspgrid = Table.read('ALIST_padova_may2020.fits')
ssp_spec = fits.open('ALIST_padova_may2020.fits')[2].data

age, feh, afe = 10, 0.0, 0.0 ##Models you want to get

model_id = np.where((sspgrid['AGE'] == age) & (sspgrid['FE_H'] == feh) & (sspgrid['A_FE'] == afe))[0]
model_spec = ssp_spec[model_id][0]

wave0 = 4.179
space = 5.880565355183108e-6
n = 8575
waven = wave0+(space*n)
wavelength = 10**(np.linspace(wave0, waven, num = n))

plt.plot(wavelength, model_spec)
plt.xlabel(r'$\lambda$($\AA$)')
plt.ylabel('flux')
plt.title(r'age = %d, [Fe/H] = %d, [$\alpha$/Fe] = %d'%(sspgrid['AGE'][model_id], sspgrid['FE_H'][model_id],sspgrid['A_FE'][model_id]))
plt.show()