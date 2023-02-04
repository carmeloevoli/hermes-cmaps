import healpy as hp
import fitsio
import numpy as np

import my_masks as masks

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
    GeV = 1.60218e-10
    fits = fitsio.FITS(fits_map_filename)
    size = len(fits)
    print ("energy array size : %i" % (size - 1))
    E = []
    for i in range(1,size):
        h = fitsio.read_header(fits_map_filename, ext=i)
        E.append(h["ENERGY"] / GeV)
    print ("energy range : %3.1e - %3.1e" % (np.min(E), np.max(E)))
    return np.array(E)

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
    
def compute_spectrum(fits_map_filename, mask, NSIDE):
    E = get_energy_array(fits_map_filename)
    spectrum = get_masked_spectrum(fits_map_filename, mask, NSIDE)
    return E, spectrum
    
def write_spectrum(E, spectrum, filename):
    f = open(filename, "w")
    f.write("#\n")
    for i in range(len(E)):
        f.write("%10.5e %10.5e\n" % (E[i], spectrum[i]))
    f.close()

def convert_fit_to_spectrum(filename):
    fits_map_filename = '../build/' + filename + '.fits.gz'
    NSIDE = get_header_info(fits_map_filename)
    full = masks.make_fullsky_mask(NSIDE)
    inner = masks.make_inner_mask(NSIDE)
    E, spectrum = compute_spectrum(fits_map_filename, full, NSIDE)
    write_spectrum(E, spectrum, filename + '_full.txt')
    E, spectrum = compute_spectrum(fits_map_filename, inner, NSIDE)
    write_spectrum(E, spectrum, filename + '_inner.txt')

if __name__== "__main__":
    filename = 'spectrum-KA-Pi0-HI-128'
    convert_fit_to_spectrum(filename)

    filename = 'spectrum-AAfrag-Pi0-HI-128'
    convert_fit_to_spectrum(filename)
