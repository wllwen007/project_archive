import matplotlib.pyplot as plt
import numpy as np
import aplpy as ap
import utils.plotting as pp
import utils.make_subim as mm
import astropy.io.fits as pf
from astropy.wcs import WCS
from astropy.table import Table, Column
import pyregion 

import configparser
config = configparser.ConfigParser()
config.read('image.config')
dpath = config['SETUP']['dpath']
plotpath = config['SETUP']['plotpath']
imfits = config['image']['imfits']
rmsfits = config['image']['rmsfits']
starmask = config['deeptable']['starmask']



dat = pf.getdata(imfits)
head = pf.getheader(imfits)
dat = dat[np.isfinite(dat)]

ra = head['CRVAL1']
dec = head['CRVAL2']
ps = head['CDELT2']
t = mm.extract_subim(imfits, ra,dec, 2.85)
trms = mm.extract_subim(rmsfits, ra,dec, 2.85)

tstarmask = mm.extract_subim(starmask, ra,dec, 2.85)

wrms = WCS(trms[0].header)
wmask = WCS(tstarmask[0].header)


reg = pyregion.open('bootes_2degcirc.reg')
mask2 = reg.get_mask(trms[0])
reg = pyregion.open('bootes_1degcirc.reg')
mask1 = reg.get_mask(trms[0])

trms[0].data = trms[0].data*1000.

#sys.exit()
f = plt.figure(figsize=(6.64,6.))
ax1 = ap.FITSFigure(trms[0],figure=f)
ax1.ticks.set_color('k') 
ax1.show_colorscale(stretch='linear', cmap=plt.cm.cubehelix_r)
ax1.show_contour(levels=[0.9,1.8],colors='k')
ax1.add_colorbar()
ax1.colorbar.set_axis_label_text('RMS (mJy/bm)') 
ax1.recenter(ra,dec,width=5.6,height=5.6)

#f.tight_layout()
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_rms')


trmsdat = trms[0].data


xin,yin = np.where(np.isnan(tstarmask[0].data))
ra,dec = wmask.all_pix2world(xin,yin,1)
xin1,yin1 = wrms.all_world2pix(ra,dec,1)
xin1 = np.array(xin1,dtype=int)
yin1 = np.array(yin1,dtype=int)
trmsdatmask = trmsdat.copy()
trmsdatmask[xin1,yin1] *= np.nan
trmsdatmask = trmsdatmask.flatten()


#f1,ax1 = pp.paper_single_ax(projection=wrms)
#f,ax = pp.paper_single_ax(projection=wmask)  
#ax1.imshow(trmsdatmask)   
#ax.imshow(tstarmask[0].data)       
#sys.exit()




trmsdat_rin = trmsdat.copy()
trmsdat_rin[~mask2] *=np.nan
trmsdat_rin = trmsdat_rin.flatten()
trmsdat_rout = trmsdat.copy()
trmsdat_rout[mask2] *=np.nan
trmsdat_rout = trmsdat_rout.flatten()


trmsdat_rin1 = trmsdat.copy()
trmsdat_rin1[~mask1] *=np.nan
trmsdat_rin1 = trmsdat_rin1.flatten()
trmsdat_rout1 = trmsdat.copy()
trmsdat_rout1[mask1] *=np.nan
trmsdat_rout1 = trmsdat_rout1.flatten()


trsmdat = trmsdat.flatten()
#xrms = np.arange(0.6, 3.0, 0.01)
#xrms1 = np.arange(3., 20., 1.)
#xrms = np.hstack((xrms,xrms1)) 


A = np.sum(np.isfinite(trmsdat))*ps*ps
xrms = np.linspace(0.4, 3.5, 100)
Arms = np.array([np.sum(trmsdat < xrmsi) for xrmsi in xrms]) * ps*ps
Afrac = Arms/A




### RMS is in mJy

f,ax = pp.paper_single_ax(AR=0.9) 
ax.plot(xrms, Arms)
ax2 = plt.twinx()
ax2.plot(xrms, Afrac)
ax2.set_ylabel('Fractional Area $<\sigma$')
plt.subplots_adjust(right=0.85)
ax2.minorticks_on()     
pp.set_attrib(ax, xlabel='RMS (mJy/bm)', ylabel='Area $<\sigma$ (deg$^2$)', xlim=(0.6, 3.0))
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_rms_area')


# save table rms versus area covered (used for completeness, etc...)
trmstab = Table([Column(data=Arms, name='Area', unit='deg2'),
                 Column(data=Afrac, name='Area_frac'),
                 Column(data=xrms, name='rms',unit='mJy/beam') ])

trmstab.write(rmsfits.replace('.rms.fits', '.Area_rms.fits'), overwrite=True)


## starmask

A = np.sum(np.isfinite(trmsdatmask))*ps*ps
xrms = np.linspace(0.4, 3.5, 100)
Arms = np.array([np.sum(trmsdatmask < xrmsi) for xrmsi in xrms]) * ps*ps
Afrac = Arms/A

f,ax = pp.paper_single_ax(AR=0.9) 
ax.plot(xrms, Arms)
ax2 = plt.twinx()
ax2.plot(xrms, Afrac)
ax2.set_ylabel('Fractional Area $<\sigma$')
plt.subplots_adjust(right=0.85)
ax2.minorticks_on()     
pp.set_attrib(ax, xlabel='RMS (mJy/bm)', ylabel='Area $<\sigma$ (deg$^2$)', xlim=(0.6, 3.0))
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_rms_area_starmask')


# save table rms versus area covered (used for completeness, etc...)
trmstab = Table([Column(data=Arms, name='Area', unit='deg2'),
                 Column(data=Afrac, name='Area_frac'),
                 Column(data=xrms, name='rms',unit='mJy/beam') ])

trmstab.write(rmsfits.replace('.rms.fits', '.Area_rms_starmask.fits'), overwrite=True)



for rad in [1,2]:
    ## inner

    if rad == 1:
        trmsdati = trmsdat_rin
        trmsdato = trmsdat_rout
    A = np.sum(np.isfinite(trmsdati))*ps*ps
    xrms = np.linspace(0.4, 3.5, 100)
    Arms = np.array([np.sum(trmsdati < xrmsi) for xrmsi in xrms]) * ps*ps
    Afrac = Arms/A




    ### RMS is in mJy

    f,ax = pp.paper_single_ax(AR=0.9) 
    ax.plot(xrms, Arms)
    ax2 = plt.twinx()
    ax2.plot(xrms, Afrac)
    ax2.set_ylabel('Fractional Area $<\sigma$')
    plt.subplots_adjust(right=0.85)
    ax2.minorticks_on()     
    pp.set_attrib(ax, xlabel='RMS (mJy/bm)', ylabel='Area $<\sigma$ (deg$^2$)', xlim=(0.6, 3.0))
    pp.fig_save_many(f, plotpath+f'bootes_deep_lba_image_rms_area_rin{rad:d}')


    # save table rms versus area covered (used for completeness, etc...)
    trmstab = Table([Column(data=Arms, name='Area', unit='deg2'),
                    Column(data=Afrac, name='Area_frac'),
                    Column(data=xrms, name='rms',unit='mJy/beam') ])

    trmstab.write(rmsfits.replace('.rms.fits', f'.Area_rms_rin{rad:d}.fits'), overwrite=True)


    ## outer


    A = np.sum(np.isfinite(trmsdato))*ps*ps
    xrms = np.linspace(0.4, 3.5, 100)
    Arms = np.array([np.sum(trmsdato < xrmsi) for xrmsi in xrms]) * ps*ps
    Afrac = Arms/A




    ### RMS is in mJy

    f,ax = pp.paper_single_ax(AR=0.9) 
    ax.plot(xrms, Arms)
    ax2 = plt.twinx()
    ax2.plot(xrms, Afrac)
    ax2.set_ylabel('Fractional Area $<\sigma$')
    plt.subplots_adjust(right=0.85)
    ax2.minorticks_on()     
    pp.set_attrib(ax, xlabel='RMS (mJy/bm)', ylabel='Area $<\sigma$ (deg$^2$)', xlim=(0.6, 3.0))
    pp.fig_save_many(f, plotpath+f'bootes_deep_lba_image_rms_area_rout{rad:d}')


    # save table rms versus area covered (used for completeness, etc...)
    trmstab = Table([Column(data=Arms, name='Area', unit='deg2'),
                    Column(data=Afrac, name='Area_frac'),
                    Column(data=xrms, name='rms',unit='mJy/beam') ])

    trmstab.write(rmsfits.replace('.rms.fits', f'.Area_rms_rout{rad:d}.fits'), overwrite=True)
