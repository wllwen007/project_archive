import os, sys
import utils.plotting as pp
from utils.make_subim import extract_subim
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import leastsq 
import astropy.coordinates as ac
from astropy.table import Table, Column, MaskedColumn
import astropy.units as u
import astropy.io.fits as fits
from astropy.wcs import WCS
from scipy import ndimage

from facet_offsets import RegPoly

def med_in_bins(x,y,bins=10):
    xt = x.copy()
    xt.sort()
    xt = xt[10:-10]
    xbins = np.linspace(np.min(xt), np.max(xt),bins+1)
    mx = np.zeros(bins)
    my = np.zeros(bins)
    mys = np.zeros(bins)
    mxs = np.zeros((2,bins))
    for i in range(bins):
        mx[i] = (xbins[i]+xbins[i+1]) /2.
        s = (x > xbins[i]) & (x <= xbins[i+1])
        my[i] = np.median(y[s])
        mys[i] = 1.2533 * np.std(y[s]) / np.sqrt(np.sum(s))
    for i in range(bins):
        mxs[0,i] = xbins[i+1]-mx[i]
        mxs[1,i] = mx[i]-xbins[i]
    
    return mx,my,mys,mxs

#clobber = True
clobber = 0
stilts = '/net/beerze/data2/wwilliams/software/bin/stilts '


import configparser
config = configparser.ConfigParser()
config.read('image.config')
plotpath = config['SETUP']['plotpath']
dpath = config['SETUP']['dpath']
regfile = config['fluxes']['regfile']
image = config['fluxes']['image']

maskfits = 'image_full_phase_m.NS.mask01.fits'

#cattable = config['fluxes']['cattable']
for cattable in [config['fluxes']['cattable'] , config['fluxes']['gcattable'] ]:
    print (cattable)
    cat = Table.read(cattable)
    #for c in cat:
        #c['RA'],c['DEC']


    hdu = fits.open(dpath+maskfits)[0]
    t = extract_subim(dpath+maskfits,217.5,34.5,10)
    w = WCS(t[0].header)

    maskdat = t[0].data
    maskdat = ndimage.binary_opening(maskdat)  #structure=np.ones((3,3)  - dilate the mask a bit - solves offset by 1 and rounding errors
    #fp=fits.open(dpath+maskfits)
    #f=fp[0].data[0,0]
    #prhd=fp[0].header

    #w = WCS(prhd)
    x,y =w.all_world2pix(cat['RA'],cat['DEC'],1)
    x = np.array(x,dtype=int)
    y = np.array(y,dtype=int)

    inmask = maskdat[y,x] 
    inmask = np.array(inmask,dtype=bool)

    pp.paper_single()
    ax1 = plt.subplot(1,1,1)
    ax1.imshow(maskdat)
    #ax1.scatter(cat['RA'],cat['DEC'])
    ax1.scatter(x[inmask],y[inmask])
    ax1.scatter(x[~inmask],y[~inmask])
    #w=wcs.WCS(prhd)



    sel = (cat['Peak_flux']/cat['Isl_rms'] > 5. ) & (cat['Maj'] < 30./3600.)
    ax1 = plt.subplot(1,1,1)
    plt.scatter(cat['Peak_flux'][inmask&sel]/cat['Isl_rms'][inmask&sel], cat['Total_flux'][inmask&sel]/cat['Peak_flux'][inmask&sel])
    plt.scatter(cat['Peak_flux'][~inmask&sel]/cat['Isl_rms'][~inmask&sel], cat['Total_flux'][~inmask&sel]/cat['Peak_flux'][~inmask&sel])


    if 'FLAG_MASKED' not in cat.columns:
        cat.add_column(Column(name='FLAG_MASKED', data=inmask))
    else:
        cat['FLAG_MASKED'] = inmask
        
    cat.write(cattable,overwrite=True)
