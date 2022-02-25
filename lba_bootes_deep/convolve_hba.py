import os
import matplotlib.pyplot as plt
import numpy as np
import aplpy as ap
import utils.plotting as pp
import utils.make_subim as mm
import utils.cutouts as cc
import astropy.io.fits as pf
import astropy.units as u
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky
from astropy.convolution import Gaussian2DKernel, convolve
from reproject import reproject_interp


clobber = False
#clobber=True

dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'
imfits = 'image_full_phase_m.NS_shift.int.facetRestored.blanked.fits'

cat = 'bootes_deep_lba_hbashift.cat.fits'
cat ='bootes_deep_lba.cat.matched.v0.1.fits'
#tcat = Table.read(dpath+cat)


gcat = 'bootes_deep_lba_hbashift.gcat.fits'
#tgcat = Table.read(dpath+gcat)

deepcat= '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.blanked.scaled.cat.fits'
deepidcat= '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/bootes_final_cross_match_catalogue-v1.0.fits'
deepim = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits'

#tdeepcat = Table.read(deepcat)
#tdeepidcat = Table.read(deepidcat)


def std_sigclip(x, nit=10, nclip=5.):
    
    x = x[np.isfinite(x)]
    
    for i in range(nit):
        std = np.std(x)
        mn = np.mean(x)
        
        x = x[ (np.abs (x-mn) - nclip*std < 0)]
    
    mn = np.mean(x)
    std = np.std(x)
    
    
    return mn, std


def get_scale(t):
    try:
        m, s = std_sigclip(t[0].data)
    except:
        return 0, 10, 1
    vmin = -1.*s +m
    vmax = 25.*s + m
    vmid = 10.*s + m
    return vmin, vmax, vmid

if 1:
    
    t = pf.getdata(dpath+imfits)
    tdeep = pf.getdata(deepim)
    
    h = pf.getheader(dpath+imfits)
    h['WCSAXES'] = 2
    h['NAXIS'] = 2
    for tkey in ['CRPIX3','CDELT3','CTYPE3','CRVAL3','CUNIT3','NAXIS3','CRPIX4','CDELT4','CTYPE4','CRVAL4','CUNIT4','NAXIS4']:
        if tkey in h.keys():
            h.remove(tkey)
    hdeep = pf.getheader(deepim)
    
    #pf.writeto(name+'hba_deep.fits', tdeep[0].data, header=hdeep, overwrite=True)
    
    '''
    convolve from
    0.001666666666666667 deg 
    0.00041666666666667 deg / pixel
    to
    0.004166666666666667 
    '''
    
    beamdeep = hdeep['BMAJ'] #/ abs(hdeep['CDELT2'])
    beam = h['BMAJ'] #/ abs(h['CDELT2'])
    conv = np.sqrt(beam**2 - beamdeep**2.)
    conv_pix_maj = conv / abs(hdeep['CDELT2'])
    
    beamdeepm = hdeep['BMIN'] #/ abs(hdeep['CDELT2'])
    beamm = h['BMIN'] #/ abs(h['CDELT2'])
    convm = np.sqrt(beamm**2 - beamdeepm**2.)
    conv_pix_min = convm / abs(hdeep['CDELT2'])
    
    conv_pix_pa = hdeep['BPA']
    kernel = Gaussian2DKernel(x_stddev=conv_pix_maj, y_stddev=conv_pix_min, theta=conv_pix_pa)
    
    
    
    #tdeep[0].data = tdeep[0].data * 1e3
    deepdat = tdeep[0].data
    
    deepdatC = convolve(deepdat, kernel)
    pf.writeto('hba_deep_convolved.fits', deepdatC, header=hdeep, overwrite=True)
    # concolvce deep to 15" resolution
    #deepdatC = 
    
    deepdatCr, footprint = reproject_interp(tdeep, h)
    hdeep = h
    
    hdeep['BMAJ'] = h['BMAJ']
    hdeep['BMIN'] = h['BMIN']
    hdeep['BPA'] = h['BPA']
    hdeep['CDELT1'] = h['CDELT1']
    hdeep['CRPIX1'] = h['CRPIX1']
    hdeep['CDELT2'] = h['CDELT2']
    hdeep['CRPIX2'] = h['CRPIX2']

    
    tdeep[0].data = deepdatCr
    
    pf.writeto('hba_deep_convolved_reproj.fits', deepdatCr, header=h, overwrite=True)
    
