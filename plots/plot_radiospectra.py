import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import healpy as hp
#import astropy.io.fits as pyfits
import fitsio
import numpy as np

import plot_lib as plib

print('Healpy version', hp.__version__)

GeV = 1.60218e-10

def get_header_info(fits_map_filename):
    h = fitsio.read_header(fits_map_filename, ext=1)
    NSIDE = h["NSIDE"]
    print ("NSIDE: %i" % NSIDE)
    print ("approx. resolution: %3.1f deg" % (hp.nside2resol(NSIDE, arcmin=True) / 60))
    print ("number of pixels: %i" % hp.nside2npix(NSIDE))
    print ("Process: %s" % h["PROCESS"])
    print ("Units: %s" % h["TUNIT1"])
    freq = h["FREQ"]
    print ("Frequency: %3.1f SI" % (freq))
    return NSIDE

def get_frequency_array(fits_map_filename):
    fits = fitsio.FITS(fits_map_filename)
    size = len(fits)
    print ("energy array size : %i" % (size - 1))
    freq = []
    for i in range(1,size):
        h = fitsio.read_header(fits_map_filename, ext=i)
        freq.append(h["FREQ"])
    print ("frequency range : %3.1e - %3.1e" % (np.min(freq), np.max(freq)))
    return np.array(freq) / 1e6

def get_fullsky_spectrum(fits_map_filename, NSIDE):
    NPIX = hp.nside2npix(NSIDE)
    fits = fitsio.FITS(fits_map_filename)
    size = len(fits)
    flux = []
    for imap in range(1, size):
        img = fits[imap].read()
        value = 0.
        for ipix in range(NPIX):
            value += img[ipix][0]
        flux.append(value / float(NPIX))
    fits.close()
    return np.array(flux)
    
def get_masked_spectrum(fits_map_filename, mask, NSIDE):
    NPIX = hp.nside2npix(NSIDE)
    assert (NPIX == len(mask))
    fits = fitsio.FITS(fits_map_filename)
    size = len(fits)
    flux = []
    visible_pixels = sum(mask)
    for imap in range(1, size):
        img = fits[imap].read()
        value = 0.
        for ipix in range(NPIX):
            value += img[ipix][0] * mask[ipix]
        flux.append(value / float(visible_pixels))
    fits.close()
    return np.array(flux)

def make_intermediate_mask(NSIDE):
    NPIX = hp.nside2npix(NSIDE)
    mask = np.zeros(NPIX)
    for ipix in range(NPIX):
        ang = hp.pix2ang(NSIDE, ipix, lonlat=True)
        b = abs(ang[1])
        l = ang[0]
        if (b > 10. and b < 45.) and (l > 40. and l < 340.):
            mask[ipix] = 1
    return mask

def make_inner_mask(NSIDE):
    NPIX = hp.nside2npix(NSIDE)
    mask = np.zeros(NPIX)
    for ipix in range(NPIX):
        ang = hp.pix2ang(NSIDE, ipix, lonlat=True)
        b = abs(ang[1])
        l = ang[0]
        if b < 5. and l < 10.:
            mask[ipix] = 1
    return mask
    
def make_spectra(fits_map_filename, spectrum_filename):
    NSIDE = get_header_info(fits_map_filename)
    freq = get_frequency_array(fits_map_filename)
    fullsky = get_fullsky_spectrum(fits_map_filename, NSIDE)
    mask = make_intermediate_mask(NSIDE)
    inner = get_masked_spectrum(fits_map_filename, mask, NSIDE)
    mask = make_inner_mask(NSIDE)
    gc = get_masked_spectrum(fits_map_filename, mask, NSIDE)

    f = open(spectrum_filename, "w")
    for i in range(len(freq)):
        f.write("%10.3e %10.3e %10.3e %10.3e\n" % (freq[i], fullsky[i], inner[i], gc[i]))
    f.close()

# MAIN

def read_spectrum_from_file(filename, icol):
    E, flux = np.loadtxt(filename, skiprows=0, usecols=(0,icol), unpack=True)
    return E, flux
    
def plot_radiospectra(icol, plot_title, plot_filename):
    fig, ax = plib.set_plot_style((4.5, 3.5))

    freq, flux = read_spectrum_from_file('spectrum-Synchro-64.txt', icol)
    ax.plot(freq, freq * freq * flux, label='Synchrotron (B regular + B turbulent)', color='tab:blue')

    freq, flux = read_spectrum_from_file('spectrum-Synchro-noturb-64.txt', icol)
    ax.plot(freq, freq * freq * flux, linestyle=':', label='Synchrotron (B regular only)', color='tab:blue')

    freq, flux = read_spectrum_from_file('spectrum-Synchro-absorption-64.txt', icol)
    ax.plot(freq, freq * freq * flux, linestyle='--', label='Synchrotron (including absorption)', color='tab:orange')

    freq, flux = read_spectrum_from_file('spectrum-FreeFree-64.txt', icol)
    ax.plot(freq, freq * freq * flux, color='tab:red', label='Free-free')

    ax.set_xscale('log')
    ax.set_xlabel(r'$\nu$ [MHz]', fontsize=12)
    ax.set_xlim([1, 1e5])

    ax.set_yscale('log')
    ax.set_ylabel(r'$\nu^2$ T$_{\rm b}$ [MHz$^2$ K]', fontsize=12)
    ax.set_ylim([1e2, 1e8])

    ax.set_title(plot_title, fontsize=11)

    ax.legend(fontsize=9, loc='lower left')
    plib.savefig(plt, plot_filename)
    
### MAiN ###
#make_spectra('output/spectrum-Synchro-noturb-64.fits.gz', 'spectrum-Synchro-noturb-64.txt')
#make_spectra('output/spectrum-Synchro-64.fits.gz', 'spectrum-Synchro-64.txt')
#make_spectra('output/spectrum-FreeFree-64.fits.gz', 'spectrum-FreeFree-64.txt')
#make_spectra('output/spectrum-Synchro-absorption-64.fits.gz', 'spectrum-Synchro-absorption-64.txt')

plot_radiospectra(1, 'Full Sky', 'radio_spectrum_fullysky')
plot_radiospectra(2, 'Intermediate Latitudes', 'radio_spectrum_intermediate')
