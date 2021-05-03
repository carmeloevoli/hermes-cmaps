import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import healpy as hp
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

def convert_map_2D(map, nside, b, l):
    b_size, l_size = len(b), len(l)
    map_2d = np.zeros((b_size, l_size))
    for i in range(b_size):
        for j in range(l_size):
            ipix = hp.ang2pix(nside, l[j], b[i], lonlat=True)
            map_2d[i][j] = map[ipix]
    return map_2d
    
def make_cbar(fig, image, units):
    cb = fig.colorbar(image, orientation='horizontal', pad=0.15)
    cb.ax.set_xlabel(units, fontsize=10)
    cb.ax.labelpad = 2
    cb.ax.tick_params('both', length=0, width=1., which='major', pad=4, bottom=True, top=True, left=True, right=True)
    cb.outline.set_linewidth(0.8)

def plot_map(fits_map_filename, output_filename, title, units, min_map, max_map, doLog=True):
    fig, ax = plib.set_plot_style((4.5, 4))

    nside = get_header_info(fits_map_filename)
    map = get_map(fits_map_filename)

    b = np.linspace(-80., +80., 160 * 2)
    l = np.linspace(-150., +150., 300 * 2)

    map_2d = convert_map_2D(map, nside, b, l)
    if doLog:
        map_2d = np.log10(map_2d)

    image = ax.pcolor(l, b, map_2d,
                      cmap='jet',
                      vmin=min_map,
                      vmax=max_map,
                      shading='auto',
                      edgecolors='face')

    make_cbar(fig, image, units)

    ax.invert_xaxis()

    ax.set_title(title, pad=5, fontsize=11)
    #ax.grid(True)
    ax.set_xlabel(r'l [deg]')
    ax.set_ylabel(r'b [deg]')
    plib.savefig(plt, output_filename)

# MAiN
plot_map('output/map-DM-128.fits.gz', 'DispersionMeasure-cartesian-128', r'Dispersion Measure', r'DM [parsec/cm${^3}$]', 0., 2500., False)
plot_map('output/map-RM-128.fits.gz', 'RotationMeasure-cartesian-128', r'Rotation Measure', r'RM [rad/m$^2$]', -100., 100., False)
plot_map('output/map-RM-noturb-128.fits.gz', 'RotationMeasure-noturb-cartesian-128', r'Rotation Measure - B regular field', r'RM [rad/m$^2$]', -100., 100., False)
plot_map('output/map-Synchro-408MHz-128.fits.gz', 'Synchrotron-408MHz-cartesian-128', r'Synchrotron Intensity', r'log T$_{\rm b}$ [K]', -1, 2)



