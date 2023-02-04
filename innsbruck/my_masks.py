import healpy as hp
import numpy as np

def make_fullsky_mask(NSIDE):
    NPIX = hp.nside2npix(NSIDE)
    mask = np.ones(NPIX)
    return mask

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
