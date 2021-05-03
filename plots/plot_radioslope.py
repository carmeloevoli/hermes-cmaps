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
    return NSIDE

def get_map(fits_map_filename):
    h = fitsio.read_header(fits_map_filename, ext=1)
    n_entries = h["NAXIS2"]
    NSIDE = h["NSIDE"]
    assert(n_entries == hp.nside2npix(NSIDE))
    fits = fitsio.FITS(fits_map_filename, iter_row_buffer=10000)
    hmap = []
    for i in range(n_entries):
        flux = fits[1][i][0]
        hmap.append(flux)
    print ("Read map from %3.1e to %3.1e with %i pixels" % (min(hmap), max(hmap), n_entries))
    print ("Mean flux: ", np.mean(hmap))
    print ("Total flux: ", sum(hmap))
    return np.array(hmap)

def compute_map_slope(map_nu1, map_nu2, nu1, nu2, b, l):
    b_size, l_size = len(b), len(l)
    slopes = np.zeros((b_size, l_size))
    for i in range(b_size):
        for j in range(l_size):
            ipix = hp.ang2pix(nside, l[j], b[i], lonlat=True)
            slopes[i][j] = (np.log(map_nu1[ipix]) - np.log(map_nu2[ipix])) / (np.log(nu1) - np.log(nu2))
    return slopes
    
def make_cbar(image, units):
    cb = fig.colorbar(image, orientation='horizontal', pad=0.15)
    cb.ax.set_xlabel(units, fontsize=10)
    cb.ax.labelpad = 2
    cb.ax.tick_params('both', length=0, width=1., which='major', pad=4, bottom=True, top=True, left=True, right=True)
    cb.outline.set_linewidth(0.8)
    
# MAIN

output_filename = 'SynchrotronSlope-408MHz-cartesian-128'
title = r'Synchrotron Slope'
units = r'$\beta$'
min_map, max_map = -3.1, -2.9

fig, ax = plib.set_plot_style((4.5, 4))

nside = get_header_info('fits/map-Synchro-408MHz-128.fits.gz')

map_nu1 = get_map('fits/map-Synchro-noturb-408MHz-128.fits.gz')
map_nu2 = get_map('fits/map-Synchro-noturb-412MHz-128.fits.gz')
nu1, nu2 = 408., 412.

b = np.linspace(-80., +80., 160 * 2)
l = np.linspace(-150., +150., 300 * 2)

map_2d = compute_map_slope(map_nu1, map_nu2, nu1, nu2, b, l)

image = ax.pcolor(l, b, map_2d,
                  cmap='jet',
                  vmin=min_map,
                  vmax=max_map,
                  shading='auto',
                  edgecolors='face')

make_cbar(image, units)

ax.invert_xaxis()

ax.set_title(title, pad=5, fontsize=11)
#ax.grid(True)
ax.set_xlabel(r'l [deg]')
ax.set_ylabel(r'b [deg]')

plib.savefig(plt, output_filename)

