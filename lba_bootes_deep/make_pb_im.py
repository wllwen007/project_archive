import astropy.io.fits as pf
import numpy as np


import configparser
config = configparser.ConfigParser()
config.read('image.config')
dpath = config['SETUP']['dpath']

appimfits = config['pbcor']['appimfits']
intimfits = config['pbcor']['intimfits']
pbimfits = config['pbcor']['pbimfits']
bappimfits = config['pbcor']['bappimfits']
bintimfits = config['pbcor']['bintimfits']


if 1:
    apim = pf.getdata(dpath+appimfits)
    aphead = pf.getheader(dpath+appimfits)
    corim = pf.getdata(dpath+intimfits)
    corhead = pf.getheader(dpath+intimfits)

    pbim = apim/corim


    pf.writeto(dpath+pbimfits, data=pbim, header=aphead, clobber=True)


    mask = (pbim < 0.3)

    apim[mask] = np.nan
    corim[mask] = np.nan

    pf.writeto(dpath+bappimfits, data=apim, header=aphead, clobber=True)
    pf.writeto(dpath+bintimfits, data=corim, header=corhead, clobber=True)

