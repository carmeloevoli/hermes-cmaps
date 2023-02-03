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
    E_gamma = h["ENERGY"] / GeV
    print ("Energy: %3.1f GeV" % (E_gamma))
    return NSIDE

def get_energy_array(fits_map_filename):
    fits = fitsio.FITS(fits_map_filename)
    size = len(fits)
    print ("energy array size : %i" % (size - 1))
    E = []
    for i in range(1,size):
        h = fitsio.read_header(fits_map_filename, ext=i)
        E.append(h["ENERGY"] / GeV)
    print ("energy range : %3.1e - %3.1e" % (np.min(E), np.max(E)))
    return np.array(E)

def get_fullsky_spectrum(fits_map_filename, NSIDE):
    NPIX = hp.nside2npix(NSIDE)
    fits = fitsio.FITS(fits_map_filename)
    size = len(fits)
    flux = []
    for imap in range(1, size):
        img = fits[imap].read()
        value = 0.
        for ipix in range(NPIX):
            if img[ipix][0] > 0:
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
            if img[ipix][0] > 0:
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
        if b > 10. and b < 40.:
            mask[ipix] = 1
    return mask
    
def make_inner_mask(NSIDE):
    NPIX = hp.nside2npix(NSIDE)
    mask = np.zeros(NPIX)
    for ipix in range(NPIX):
        ang = hp.pix2ang(NSIDE, ipix, lonlat=True)
        b = abs(ang[1])
        l = ang[0]
        if b < 8. and (l < 80. or l > 280.):
            mask[ipix] = 1
    return mask
    
def make_outer_mask(NSIDE):
    NPIX = hp.nside2npix(NSIDE)
    mask = np.zeros(NPIX)
    for ipix in range(NPIX):
        ang = hp.pix2ang(NSIDE, ipix, lonlat=True)
        b = abs(ang[1])
        l = ang[0]
        if b < 8. and (l > 80. and l < 280.):
            mask[ipix] = 1
    return mask

def make_spectra(fits_map_filename, spectrum_filename):
    NSIDE = get_header_info(fits_map_filename)
    E = get_energy_array(fits_map_filename)
    fullsky = get_fullsky_spectrum(fits_map_filename, NSIDE)
    mask = make_intermediate_mask(NSIDE)
    intermediate = get_masked_spectrum(fits_map_filename, mask, NSIDE)
    mask = make_inner_mask(NSIDE)
    inner = get_masked_spectrum(fits_map_filename, mask, NSIDE)
    mask = make_outer_mask(NSIDE)
    outer = get_masked_spectrum(fits_map_filename, mask, NSIDE)
    
    f = open(spectrum_filename, "w")
    for i in range(len(E)):
        f.write("%10.3e %10.3e %10.3e %10.3e %10.3e\n" % (E[i], fullsky[i], intermediate[i], inner[i], outer[i]))
    f.close()

def read_flux_from_file(filename, icol):
    E, flux = np.loadtxt(filename, skiprows=0, usecols=(0,icol), unpack=True)
    return E, flux
    
def plot_flux(icol, plot_title, plot_filename):
    fig, ax = plib.set_plot_style((4.5, 3.5))

    E, flux_HI = read_flux_from_file('spectrum-PionDecay-HI-100MeV-1TeV-64.txt', icol)
    ax.plot(E, E * E * flux_HI, color='tab:green', linestyle='--', label='Pion Decay (HI)')

    E, flux_H2 = read_flux_from_file('spectrum-PionDecay-H2-100MeV-1TeV-64.txt', icol)
    ax.plot(E, E * E * flux_H2, color='tab:green', linestyle=':', label='Pion Decay (H2)')

    ax.plot(E, E * E * (flux_HI + flux_H2), color='tab:green', linestyle='-', label='Pion Decay (total)')
    
    E, flux_IC = read_flux_from_file('spectrum-InverseCompton-100MeV-1TeV-64.txt', icol)
    ax.plot(E, E * E * flux_IC, color='tab:red', linestyle='-', label='Inverse Compton')

    E, flux_Bremss_HI = read_flux_from_file('spectrum-Bremsstrahlung-HI-100MeV-1TeV-64.txt', icol)
    E, flux_Bremss_H2 = read_flux_from_file('spectrum-Bremsstrahlung-HI-100MeV-1TeV-64.txt', icol)
    ax.plot(E, E * E * (flux_Bremss_HI + flux_Bremss_H2), color='tab:orange', linestyle='-', label='Bremsstrahlung')

    flux = flux_Bremss_HI + flux_Bremss_H2 + flux_HI + flux_H2 + flux_IC
    ax.plot(E, E * E * flux, color='tab:blue', linestyle='-', label='Total')

    ax.set_xscale('log')
    ax.set_xlabel(r'E [GeV]', fontsize=12)
    ax.set_xlim([0.1, 1e3])

    ax.set_yscale('log')
    ax.set_ylabel(r'E$^2$ Flux [GeV m$^{-2}$ s$^{-1}$ sr$^{-1}$]', fontsize=12)
    ax.set_ylim([1e-5, 1])

    ax.set_title(plot_title, fontsize=11)

    ax.legend(fontsize=7)
    plib.savefig(plt, plot_filename)
    
def plot_flux_HE(icol, plot_title, plot_filename):
    fig, ax = plib.set_plot_style((4.5, 3.5))

    E, flux_HI = read_flux_from_file('spectrum-PionDecay-HI-1TeV-1PeV-64.txt', icol)
    E, flux_H2 = read_flux_from_file('spectrum-PionDecay-H2-1TeV-1PeV-64.txt', icol)
    ax.plot(E, E * E * (flux_HI + flux_H2), color='tab:green', linestyle='-', label=r'Pion Decay ($\gamma$)')

    E, flux_HI = read_flux_from_file('spectrum-PionDecayNu-HI-1TeV-1PeV-64.txt', icol)
    E, flux_H2 = read_flux_from_file('spectrum-PionDecayNu-H2-1TeV-1PeV-64.txt', icol)
    ax.plot(E, E * E * (flux_HI + flux_H2), color='tab:orange', linestyle='--', label=r'Pion Decay ($\nu$ all-flavours)')

    E, flux_DM = read_flux_from_file('spectrum-DarkMatter-1TeV-30TeV-256.txt', icol)
    ax.plot(E, E * E * 10. * flux_DM, color='tab:red', linestyle='-', label='Dark Matter (prompt)')

    E, flux_DM = read_flux_from_file('spectrum-DarkMatterSecondary-1TeV-30TeV-64.txt', icol)
    ax.plot(E, E * E * 10. * flux_DM, color='tab:red', linestyle='--', label='Dark Matter (secondary)')

    ax.set_xscale('log')
    ax.set_xlabel(r'E [GeV]', fontsize=12)
    ax.set_xlim([1e3, 1e6])

    ax.set_yscale('log')
    ax.set_ylabel(r'E$^2$ Flux [GeV m$^{-2}$ s$^{-1}$ sr$^{-1}$]', fontsize=12)
    ax.set_ylim([1e-8, 1e-2])

    ax.set_title(plot_title, fontsize=11)

    ax.legend(fontsize=7)
    plib.savefig(plt, plot_filename)

def plot_flux_absorption(plot_title, plot_filename):
    fig, ax = plib.set_plot_style((4.5, 3.5))

    icol = 3
    
    E, flux_HI = read_flux_from_file('spectrum-PionDecayAbsorption-HI-1TeV-1PeV-16.txt', icol)
    E, flux_H2 = read_flux_from_file('spectrum-PionDecayAbsorption-H2-1TeV-1PeV-16.txt', icol)
    ax.plot(E, E * E * (flux_HI + flux_H2), color='tab:blue', linestyle='-', label=r'Pion Decay')

    E, flux_HI = read_flux_from_file('spectrum-PionDecay-HI-1TeV-1PeV-16.txt', icol)
    E, flux_H2 = read_flux_from_file('spectrum-PionDecay-H2-1TeV-1PeV-16.txt', icol)
    ax.plot(E, E * E * (flux_HI + flux_H2), color='tab:blue', linestyle='--', label=r'Pion Decay (no absorption)')

    ax.set_xscale('log')
    ax.set_xlabel(r'E [GeV]', fontsize=12)
    ax.set_xlim([1e3, 1e6])

    ax.set_yscale('log')
    ax.set_ylabel(r'E$^2$ Flux [GeV m$^{-2}$ s$^{-1}$ sr$^{-1}$]', fontsize=12)
    ax.set_ylim([1e-7, 1e-2])

    ax.set_title(plot_title, fontsize=11)

    ax.legend(fontsize=7)
    plib.savefig(plt, plot_filename)
    
def plot_flux_absorption_ratio(plot_title, plot_filename):
    fig, ax = plib.set_plot_style((4.5, 3.5))

    icol = 3
    
    E, flux_HI = read_flux_from_file('spectrum-PionDecayAbsorption-HI-1TeV-1PeV-16.txt', icol)
    E, flux_H2 = read_flux_from_file('spectrum-PionDecayAbsorption-H2-1TeV-1PeV-16.txt', icol)
    y_w = flux_HI + flux_H2
    
    E, flux_HI = read_flux_from_file('spectrum-PionDecay-HI-1TeV-1PeV-16.txt', icol)
    E, flux_H2 = read_flux_from_file('spectrum-PionDecay-H2-1TeV-1PeV-16.txt', icol)
    y_no = flux_HI + flux_H2
    
    r = (y_w - y_no) / y_no
    ax.plot(E, r - max(r))
    
    ax.set_xscale('log')
    ax.set_xlabel(r'E [GeV]', fontsize=12)
    ax.set_xlim([1e3, 1e6])

#    ax.set_yscale('log')
    ax.set_ylabel(r'Flux relative ratio', fontsize=12)
    ax.set_ylim([-1, 0])

    ax.set_title(plot_title, fontsize=11)

    ax.legend(fontsize=7)
    plib.savefig(plt, plot_filename)
    
### LE Energy ###

#make_spectra('output/spectrum-PionDecay-HI-100MeV-1TeV-64.fits.gz', 'spectrum-PionDecay-HI-100MeV-1TeV-64.txt')
#make_spectra('output/spectrum-PionDecay-H2-100MeV-1TeV-64.fits.gz', 'spectrum-PionDecay-H2-100MeV-1TeV-64.txt')
#make_spectra('output/spectrum-InverseCompton-100MeV-1TeV-64.fits.gz', 'spectrum-InverseCompton-100MeV-1TeV-64.txt')
#make_spectra('output/spectrum-Bremsstrahlung-HI-100MeV-1TeV-64.fits.gz', 'spectrum-Bremsstrahlung-HI-100MeV-1TeV-64.txt')
#make_spectra('output/spectrum-Bremsstrahlung-H2-100MeV-1TeV-64.fits.gz', 'spectrum-Bremsstrahlung-H2-100MeV-1TeV-64.txt')

#plot_flux(1, 'Full Sky', 'gamma_spectrum_belowtev_fullysky')
#plot_flux(2, 'Intermediate Latitudes', 'gamma_spectrum_belowtev_intermediate')
#plot_flux(3, 'Inner Galaxy', 'gamma_spectrum_belowtev_inner')
#plot_flux(4, 'Outer Galaxy', 'gamma_spectrum_belowtev_outer')

### HE Energy ###

#make_spectra('fits/spectrum-PionDecay-HI-1TeV-1PeV-64.fits.gz', 'spectrum-PionDecay-HI-1TeV-1PeV-64.txt')
#make_spectra('fits/spectrum-PionDecay-H2-1TeV-1PeV-64.fits.gz', 'spectrum-PionDecay-H2-1TeV-1PeV-64.txt')
#make_spectra('fits/spectrum-PionDecayNu-HI-1TeV-1PeV-64.fits.gz', 'spectrum-PionDecayNu-HI-1TeV-1PeV-64.txt')
#make_spectra('fits/spectrum-PionDecayNu-H2-1TeV-1PeV-64.fits.gz', 'spectrum-PionDecayNu-H2-1TeV-1PeV-64.txt')
#make_spectra('fits/spectrum-DarkMatter-1TeV-30TeV-64.fits.gz', 'spectrum-DarkMatter-1TeV-30TeV-64.txt')
#make_spectra('fits/spectrum-DarkMatter-1TeV-30TeV-256.fits.gz', 'spectrum-DarkMatter-1TeV-30TeV-256.txt')
#make_spectra('fits/spectrum-DarkMatterSecondary-1TeV-30TeV-64.fits.gz', 'spectrum-DarkMatterSecondary-1TeV-30TeV-64.txt')
#make_spectra('output/spectrum-PionDecayAbsorption-HI-1TeV-1PeV-16.fits.gz', 'spectrum-PionDecayAbsorption-HI-1TeV-1PeV-16.txt')

#plot_flux_HE(1, 'Full Sky', 'gamma_spectrum_abovetev_fullysky')
#plot_flux_HE(2, 'Intermediate Latitudes', 'gamma_spectrum_abovetev_intermediate')
#plot_flux_HE(3, 'Inner Galaxy', 'gamma_spectrum_abovetev_inner')
#plot_flux_HE(4, 'Outer Galaxy', 'gamma_spectrum_abovetev_outer')

#make_spectra('fits/spectrum-PionDecayAbsorption-HI-1TeV-1PeV-16.fits.gz', 'spectrum-PionDecayAbsorption-HI-1TeV-1PeV-16.txt')
#make_spectra('fits/spectrum-PionDecay-HI-1TeV-1PeV-16.fits.gz', 'spectrum-PionDecay-HI-1TeV-1PeV-16.txt')

#make_spectra('fits/spectrum-PionDecayAbsorption-H2-1TeV-1PeV-16.fits.gz', 'spectrum-PionDecayAbsorption-H2-1TeV-1PeV-16.txt')
#make_spectra('fits/spectrum-PionDecay-H2-1TeV-1PeV-16.fits.gz', 'spectrum-PionDecay-H2-1TeV-1PeV-16.txt')

#plot_flux_absorption('Inner Galaxy', 'gamma_spectrum_abovetev_absorption')
#plot_flux_absorption_ratio('Inner Galaxy', 'gamma_ratio_abovetev_absorption')

#make_spectra('fitsfiles/spectrum-PionDecay-HI-1TeV-1PeV-64.fits.gz', 'spectrum-PionDecay-HI-1TeV-1PeV-64.txt')
#make_spectra('fitsfiles/spectrum-PionDecayWa-HI-1TeV-1PeV-64.fits.gz', 'spectrum-PionDecayWa-H2-1TeV-1PeV-64.txt')

fig, ax = plib.set_plot_style((4.5, 3.5))

icol = 3

E, flux_HI = read_flux_from_file('spectrum-PionDecay-HI-1TeV-1PeV-64.txt', icol)
E2 = np.power(E, 2.5)
ax.plot(E, E2 * (flux_HI), color='tab:blue', linestyle='-', label=r'Pion Decay')

E, flux_HI = read_flux_from_file('spectrum-PionDecayWa-H2-1TeV-1PeV-64.txt', icol)
E2 = np.power(E, 2.5)
ax.plot(E, E2 * (flux_HI), color='tab:blue', linestyle='--', label=r'Pion Decay (with absorption)')

ax.set_xscale('log')
ax.set_xlabel(r'E [GeV]', fontsize=12)
ax.set_xlim([1e3, 1e6])

ax.set_yscale('log')
ax.set_ylabel(r'E$^{2.5}$ Flux [GeV m$^{-2}$ s$^{-1}$ sr$^{-1}$]', fontsize=12)
ax.set_ylim([1e-3, 1e-1])

ax.legend(fontsize=7)
plib.savefig(plt, 'test_absorption')

