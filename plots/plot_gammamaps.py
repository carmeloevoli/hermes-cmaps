import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import healpy as hp
#import astropy.io.fits as pyfits
import fitsio
import numpy as np
import plot_lib as plib

import matplotlib.cm as cm

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
    return NSIDE, E_gamma

def get_map(fits_map_filename):
    h = fitsio.read_header(fits_map_filename, ext=1)
    n_entries = h["NAXIS2"]
    NSIDE = h["NSIDE"]
    assert(n_entries == hp.nside2npix(NSIDE))
    fits = fitsio.FITS(fits_map_filename, iter_row_buffer=10000)
    hmap = []
    for i in range(n_entries):
        flux = fits[1][i][0]
        hmap.append(max(flux, 1e-20))
    print ("Read map from %3.1e to %3.1e with %i pixels" % (min(hmap), max(hmap), n_entries))
    print ("Mean flux: ", np.mean(hmap))
    print ("Total flux: ", sum(hmap))
    return np.array(hmap)

def convert_map_2D(map, nside, E_gamma, b, l):
    b_size, l_size = len(b), len(l)
    map_2d = np.zeros((b_size, l_size))
    for i in range(b_size):
        for j in range(l_size):
            ipix = hp.ang2pix(nside, l[j], b[i], lonlat=True)
            map_2d[i][j] = map[ipix]
    return np.log10(E_gamma * E_gamma * map_2d)

def make_cbar(fig, image):
    cb = fig.colorbar(image, orientation='horizontal', pad=0.15)
    cb.ax.set_xlabel(r'log E$^2$ Flux [GeV cm$^{-2}$ sr$^{-1}$ s$^{-1}$]', fontsize=10)
    cb.ax.labelpad = 2
    cb.ax.tick_params('both', length=0, width=1., which='major', pad=4, bottom=True, top=True, left=True, right=True)
    cb.outline.set_linewidth(0.8)

def plot_map(fits_map_filename, output_filename, title, min_map, max_map):
    fig, ax = plib.set_plot_style((4.5, 4))
    nside, E_gamma = get_header_info(fits_map_filename)
    map = get_map(fits_map_filename)
    
    b = np.linspace(-80., +80., 160 * 2)
    l = np.linspace(-150., +150., 300 * 2)

    map_2d = convert_map_2D(map, nside, E_gamma, b, l)

    image = ax.pcolor(l, b, map_2d,
                      cmap='jet',
                      vmin=min_map,
                      vmax=max_map,
                      shading='auto',
                      edgecolors='face')

    make_cbar(fig, image)

    ax.invert_xaxis()

    ax.set_title(title, pad=5, fontsize=11)
    #ax.grid(True)
    ax.set_xlabel(r'l [deg]')
    ax.set_ylabel(r'b [deg]')
    plib.savefig(plt, output_filename)

def plot_ring(fits_map_filename, output_filename, title, cmap):
    fig, ax = plib.set_plot_style((4.5, 4))
    nside, E_gamma = get_header_info(fits_map_filename)
    map = get_map(fits_map_filename)

    b = np.linspace(-20., +20., 40 * 4)
    l = np.linspace(-150., +150., 300 * 2)

    map_2d = convert_map_2D(map, nside, E_gamma, b, l)

    max_map = np.max(map_2d)
    min_map = max_map - 2.5
    
    image = ax.pcolor(l, b, map_2d,
                      cmap=cmap,
                      vmin=min_map,
                      vmax=max_map,
                      shading='auto',
                      edgecolors='face')

    make_cbar(fig, image)

    ax.invert_xaxis()

    ax.set_title(title, pad=5, fontsize=11)
    #ax.grid(True)
    ax.set_xlabel(r'l [deg]')
    ax.set_ylabel(r'b [deg]')
    ax.set_xlim([-150, 150])
    ax.set_ylim([-20, 20])

    plib.savefig(plt, output_filename)
    #plt.savefig(output_filename + ".png")

# MAIN
#plot_map('fitsfiles/map-Pi0-HI-10GeV-256.fits.gz ', 'PionDecay-HI-10GeV-cartesian-256', r'Pion Decay HI - E$_\gamma$ = 10 GeV', -4., -1.)
#plot_map('fitsfiles/map-Pi0-H2-10GeV-256.fits.gz', 'PionDecay-H2-10GeV-cartesian-256', r'Pion Decay H2 - E$_\gamma$ = 10 GeV', -3., -0.5)

for i in range(0,11):
#    plot_ring('fitsfiles/map-ring-Pi0-HI-5GeV-256-' + str(i) + '.fits.gz', 'PionDecay-HI-5GeV-256-' + str(i), r'Pion Decay HI - E$_\gamma$ = 5 GeV', 'OrRd')
    plot_ring('fitsfiles/map-ring-Pi0-H2-5GeV-256-' + str(i) + '.fits.gz', 'PionDecay-H2-5GeV-256-' + str(i), r'Pion Decay H2 - E$_\gamma$ = 5 GeV', 'PuBuGn')

#plot_map('output/map-IC-10GeV-128.fits.gz', 'IC-10GeV-cartesian-128', r'Inverse Compton - E$_\gamma$ = 10 GeV')
