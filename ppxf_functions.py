def rebin_data(wavein, specin, noisein, wavnew, wavstart, wavend):
    """"""
    
    spec = spectrum.ArraySourceSpectrum(wave=wavein, flux=specin)
    f = np.ones(len(wavein))
    filt = spectrum.ArraySpectralElement(wavein, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew, force='taper')

    # rebin uncertainty according to new wavelength array
    specn = spectrum.ArraySourceSpectrum(wave=wavein, flux=noisein)
    fn = np.ones(len(wavein))
    filtn = spectrum.ArraySpectralElement(wavein, fn, waveunits='angstrom')
    obsn = observation.Observation(specn, filtn, binset=wavnew, force='taper')

    # Read in a data spectrum and define the wavelength range
    lamgal = wavnew 
    dataflux = obs.binflux 
    noiseflux = obsn.binflux 
    ##Within required wavelength range:
    usepixg=np.where(np.logical_and(lamgal > wavstart,lamgal < wavend))
    data_lin = dataflux[usepixg] 
    noise_lin = noiseflux[usepixg] 
    lamRange1 = [np.min(lamgal[usepixg]),np.max(lamgal[usepixg])]

    ##Log rebinning the data and noise:
    galaxy, logLam1, velscale = util.log_rebin(lamRange1,data_lin) 
    noise, logLam1, velscale = util.log_rebin(lamRange1,noise_lin) 
    ##Normalize the rebinned spectrum
    noise = noise/np.median(galaxy) 
    galaxy = galaxy/np.median(galaxy) 
    gal_wave = np.exp(logLam1)
    
    return galaxy, noise, gal_wave, velscale, logLam1


def skyline_masks(skylines, wavstart, wavend, gal_wave):
    
    lamsky = skylines[:,0]   # wavelengths of sky lines
    fluxsky = skylines[:,1]  # Theoretical intensity of lines (absolute flux arbitrary)
    usepixs=np.where(np.logical_and(lamsky > wavstart,lamsky < wavend))# define wavelength range to use                          
    lamsky1 = lamsky[usepixs] # small wavelength region
    fluxsky1 = fluxsky[usepixs] # corresponding flux region
    dl = 3 # set variable for deviation from sky lines (Angstroms)

    minskyline = fluxsky1 > 1 # define minimum strength of sky lines to flag
    skymask = lamsky1[minskyline] # corresponding wavelengths for those sky lines

    zeros = np.zeros_like(gal_wave) # zeros array with same size as galaxy wavelength used in fit
    # add 1 to array elements with same indices as galaxy wavelengths within dl of sky line wavelengths
    for i in range(len(gal_wave)):
        for j in range(len(skymask)):
            if np.logical_and(gal_wave[i] < skymask[j]+dl,gal_wave[i] > skymask[j]-dl):
                zeros[i] = zeros[i]+1
                break
            else:
                zeros[i] = 0

    # zeros array with same size as galaxy wavelength used in fit
    # add 1 to array elements in gap 1
    gap1zeros = np.zeros_like(gal_wave)
    for i in range(len(gal_wave)):
        if np.logical_and(gal_wave[i] < 1.59e4,gal_wave[i] > 1.576e4): 
        #if np.logical_and(gal_wave[i] < 1.585e4,gal_wave[i] > 1.578e4):
            gap1zeros[i] = gap1zeros[i]+1


    # zeros array with same size as galaxy wavelength used in fit
    # add 1 to array elements in gap 1
    gap2zeros = np.zeros_like(gal_wave)
    for i in range(len(gal_wave)):
        if np.logical_and(gal_wave[i] < 1.6525e4,gal_wave[i] > 1.6385e4):
        #if np.logical_and(gal_wave[i] < 1.647e4,gal_wave[i] > 1.6405e4):
            gap2zeros[i] = gap2zeros[i]+1


    # combined gaps array
    gapzeros = gap1zeros + gap2zeros

    goodPixels = np.where(np.logical_and(zeros==0,gapzeros==0))[0] # extract indices of zeros from sky and gap masks
                                                                    # these are the good pixels
        
    return goodPixels


def rebin_template_gc(stars, velscale, stel_wave, wavnew, FWHM_gal):
    

    ##Rebin on spectrum to the same velscale of the data spectrum, to find the array size for the template spectra:
    hdu1 = pyfits.open(stars[0]) 
    starfit = hdu1[2].data 
    starfit[starfit==0] = 1. # replace 0's in chip gaps with 1's
    # rebin template according to new wavelength array
    specs = spectrum.ArraySourceSpectrum(wave=stel_wave, flux=starfit)
    fs = np.ones(len(stel_wave))
    filts = spectrum.ArraySpectralElement(stel_wave, fs, waveunits='angstrom')
    obss = observation.Observation(specs, filts, binset=wavnew, force='taper')
    lamRange2 = [np.min(wavnew),np.max(wavnew)] # wavelength range from array
    temp, logLam2,velscale = util.log_rebin(lamRange2,obss.binflux, velscale = velscale)

    FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)/1.6e4*c # FWHM difference in km/s
    sigma = FWHM_dif/2.355/velscale # Sigma difference in pixels    
    
    templates = np.empty((temp.size,len(stars))) # templates array
    for j in range(len(stars)):
        hdu1 = pyfits.open(stars[j])
        starfit = hdu1[2].data
        starfit[starfit==0] = 1.
        specs = spectrum.ArraySourceSpectrum(wave=stel_wave, flux=starfit)
        fs = np.ones(len(stel_wave))
        filts = spectrum.ArraySpectralElement(stel_wave, fs, waveunits='angstrom')
        obss = observation.Observation(specs, filts, binset=wavnew, force='taper')
        starfitNew, logLam2, velscale = util.log_rebin(lamRange2, obss.binflux, velscale = velscale)
        ##normalizing templates:
        templates[:,j] = np.nan_to_num(starfitNew/np.median(starfitNew))
        
    return templates, logLam2, velscale, sigma



def rebin_template_miles(stars, velscale, stel_wave, wavnew):
    
    # Sigma difference in pixels
    sigma = np.sqrt(60**2 - 4.144**2)/velscale
    
    ##Rebin on spectrum to the same velscale of the data spectrum, to find the array size for the template spectra:
    hdu1 = pyfits.open(stars[0]) 
    starfit = hdu1[2].data 
    starfit[starfit==0] = 1. # replace 0's in chip gaps with 1's
    ##Reducing the template resolution:
    starfit = scipy.ndimage.filters.gaussian_filter(starfit, sigma, order=0, output=None, mode='reflect', cval=0.0, truncate=2.0)
    # rebin template according to new wavelength array
    specs = spectrum.ArraySourceSpectrum(wave=stel_wave, flux=starfit)
    fs = np.ones(len(stel_wave))
    filts = spectrum.ArraySpectralElement(stel_wave, fs, waveunits='angstrom')
    obss = observation.Observation(specs, filts, binset=wavnew, force='taper')
    lamRange2 = [np.min(wavnew),np.max(wavnew)] # wavelength range from array
    temp, logLam2,velscale = util.log_rebin(lamRange2,obss.binflux, velscale = velscale)

    templates = np.empty((temp.size,len(stars))) # templates array
    for j in range(len(stars)):
        hdu1 = pyfits.open(stars[j])
        starfit = hdu1[2].data
        starfit[starfit==0] = 1.
        tempfit = ndimage.gaussian_filter1d(starfit,sigma)
        #tempfit = scipy.ndimage.filters.gaussian_filter(starfit, sigma, order = 0, output=None, mode='reflect', cval=0.0, truncate=2.0)        
        specs = spectrum.ArraySourceSpectrum(wave=stel_wave, flux=tempfit)
        fs = np.ones(len(stel_wave))
        filts = spectrum.ArraySpectralElement(stel_wave, fs, waveunits='angstrom')
        obss = observation.Observation(specs, filts, binset=wavnew, force='taper')
        starfitNew, logLam2, velscale = util.log_rebin(lamRange2, obss.binflux, velscale = velscale)
        ##normalizing templates:
        templates[:,j] = np.nan_to_num(starfitNew/np.median(starfitNew))
        
    return templates, logLam2, velscale


