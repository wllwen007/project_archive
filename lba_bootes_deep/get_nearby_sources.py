import os, sys
import utils.plotting as pp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import leastsq 
import astropy.coordinates as ac
from astropy.table import Table, Column, MaskedColumn
import astropy.units as u
import astropy.io.fits as fits

from facet_offsets import RegPoly


import configparser
config = configparser.ConfigParser()
config.read('image.config')
plotpath = config['SETUP']['plotpath']
dpath = config['SETUP']['dpath']
cattable = config['fluxes']['cattable']
cattable = config['fluxes']['cattable']
regfile = config['fluxes']['regfile']
image = config['fluxes']['image']

name = sys.argv[1]

cat = Table.read(cattable)
catsel = cat[cat['Source_Name']==name]
Ccat = ac.SkyCoord(cat['RA'], cat['DEC'])
Ccatsel = ac.SkyCoord(catsel['RA'], catsel['DEC'])

catdeep = Table.read('/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/data_release/bootes/final_cross_match_catalogue-v1.0.fits')
Ccatdeep = ac.SkyCoord(catdeep['RA'], catdeep['DEC'], unit='deg')

sep = Ccatsel.separation(Ccat)
isort = np.argsort(sep)

cat = cat[isort]
sep = sep[isort]

print('closest LBA sources to',name)
for i in range(5):
    print('{0:s} {1:.2f} mJy {2:.2f} asec {3:.2f} asec'.format(cat['Source_Name'][i], cat['Total_flux'][i]*1000.,cat['Maj'][i]*3600, sep[i].to(u.arcsec).value))
    
    
    
sep = Ccatsel.separation(Ccatdeep)
isort = np.argsort(sep)

cat = catdeep[isort]
sep = sep[isort]

print('closest Deep HBA-opt sources to',name)
for i in range(5):
    print('{0:s} {1:.2f} mJy {2:.2f} asec {3:.2f} asec'.format(cat['Source_Name'][i], cat['Total_flux'][i]*1000.,cat['Maj'][i]*3600, sep[i].to(u.arcsec).value))
