import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import healpy as hp
#import astropy.io.fits as pyfits
import fitsio
import numpy as np
import plot_lib as plib

print('Healpy version', hp.__version__)

def get_header(fits_map_filename):
    h = fitsio.read_header(fits_map_filename, ext=0)
    n_lon = h["NAXIS1"]
    n_lat = h["NAXIS2"]
    n_rings = h["NAXIS3"]

    min_lon = h["CRVAL1"]
    delta_lon = h["CDELT1"]
    min_lat = h["CRVAL2"]
    delta_lat = h["CDELT2"]

    print (n_lon, n_lat, n_rings)
    print (min_lon, delta_lon, min_lat, delta_lat)
    
    return h

def get_img(fits_map_filename, imap):
    h = fitsio.read_header(fits_map_filename, ext=0)
    n_entries = h["NAXIS3"]
    assert(n_entries == 12)
    fits = fitsio.FITS(fits_map_filename, iter_row_buffer=100000)
    img = fits[0].read()
    img = img[imap]
    return img

def make_cbar(fig, image):
    cb = fig.colorbar(image, orientation='horizontal', pad=0.15)
    cb.ax.set_xlabel(r'log N$_{\rm CO}$ [cm$^{-2}$]', fontsize=10)
    cb.ax.labelpad = 2
    cb.ax.tick_params('both', length=0, width=1., which='major', pad=4, bottom=True, top=True, left=True, right=True)
    cb.outline.set_linewidth(0.8)

def plot_map(fits_map_filename, imap, output_filename, title):
    fig, ax = plib.set_plot_style((4.5, 4))
    h = get_header(fits_map_filename)
    img = get_img(fits_map_filename, imap)

    n_lon = h["NAXIS1"]
    n_lat = h["NAXIS2"]
    b = np.linspace(-90., +90., n_lat)
    l = np.linspace(-180., +190., n_lon)
    
    logImg = np.log10(img + 1e-20)

    print (np.max(logImg), np.max(logImg) - 2.5)

    image = ax.pcolor(l, b, logImg,
                      cmap='Blues',
                      vmax=np.max(logImg),
                      vmin=np.max(logImg) - 2.5,
                      shading='auto',
                      edgecolors='face')

    make_cbar(fig, image)

    #ax.invert_xaxis()

    ax.set_title(title, pad=5, fontsize=11)
    #ax.grid(True)
    ax.set_xlabel(r'l [deg]')
    ax.set_ylabel(r'b [deg]')
    plt.savefig(output_filename + ".png")

# MAIN

fits_map_filename = 'WCOrings_COGAL.fits.gz'

print (get_header(fits_map_filename))

for i in range(12):
    plot_map(fits_map_filename, i, 'CO_ring_' + str(i), 'CO')



