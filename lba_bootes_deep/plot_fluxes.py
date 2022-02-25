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
cattable = config['fluxes']['cattable']
cattable = config['fluxes']['cattable']
regfile = config['fluxes']['regfile']
image = config['fluxes']['image']


cat = Table.read(cattable)
cat.sort('Source_id')




fp=fits.open(dpath+image)
#f=fp[0].data[0,0]
prhd=fp[0].header
#w=wcs.WCS(prhd)
r = RegPoly(regfile, prhd['CRVAL1'], prhd['CRVAL2'])

r.add_facet_labels(cat)

#for ci in ['Total_flux','E_Total_flux','Peak_flux','E_Peak_flux','Isl_rms' ]:
    #cat[ci] = cat[ci]*1e3   #to mJy

#trms['Area']
#trms['rms']

# select only the 'proper' detected sources.... what to do about completeness of higher order wavelet scales??
print('Catalogue contains {0} sources'.format(len(cat)))
sel = (cat['Peak_flux'] / cat['Isl_rms']) > 5.
#cat = cat[sel]
print('Catalogue contains {0} sources after Peak/noise>5'.format(sum(sel)))

tflux = cat['Total_flux'] 
pflux = cat['Peak_flux']
rms = cat['Isl_rms']

C0 = ac.SkyCoord(218.0125, 34.2758333333,unit='deg')
sep = C0.separation(ac.SkyCoord(cat['RA'],cat['DEC']))

maskS = (cat['S_Code'] == 'S') & (sel) & (cat['Maj'] < 30./3600)
print('Catalogue contains {0} sources S, Peak/noise>5, Maj<30'.format(sum(maskS)))
maskS = (cat['S_Code'] == 'S') & (sel) & (cat['Maj'] < 30./3600) & (cat['FLAG_MASKED'])
print('Catalogue contains {0} sources S, Peak/noise>5, Maj<30, deconv'.format(sum(maskS)))

f,ax = pp.paper_single_ax()
ax.scatter(sep, tflux/pflux, marker='o',s=1,c='gray')
for i in range(5):
    mm = cat['Facet']==i
    ax.scatter(sep[maskS&mm], tflux[maskS&mm]/pflux[maskS&mm], marker='o',s=7,label=str(i))
plt.legend()
mx,my,mys,mxs = med_in_bins(sep[maskS].value, tflux[maskS]/pflux[maskS], bins=6)
ax.errorbar(mx, my, mys, mxs, c='C1',linestyle='none')
y1,y2 = ax.get_ylim()
#ax.vlines(5, y1,y2, color='k')
ax.hlines(1., 0, 5, color='k')
pp.set_attrib(ax, xlabel='distance from pointing centre (deg)', ylabel='$S_{\mathrm{total}}/S_{\mathrm{peak}}$',xlim=(0,3), ylim=(0,5))
pp.fig_save_many(f, plotpath+'flux_ratio_sep_facet')


f,ax = pp.paper_single_ax()
ax.scatter(sep, tflux/pflux, marker='o',s=1,c='gray')
ax.scatter(sep[maskS], tflux[maskS]/pflux[maskS], marker='o',s=7,c='C0')
mx,my,mys,mxs = med_in_bins(sep[maskS].value, tflux[maskS]/pflux[maskS], bins=6)
ax.errorbar(mx, my, mys, mxs, c='C1',linestyle='none')
with open('smearing_radius.txt','w') as fout:
    fout.write('sep smear\n')
    for ii in range(len(mx)):
        fout.write('{:f} {:f}\n'.format(mx[ii],my[ii]))
y1,y2 = ax.get_ylim()
#ax.vlines(5, y1,y2, color='k')
ax.hlines(1., 0, 5, color='k')
pp.set_attrib(ax, xlabel='distance from pointing centre (deg)', ylabel='$S_{\mathrm{total}}/S_{\mathrm{peak}}$',xlim=(0,3), ylim=(0,5))
pp.fig_save_many(f, plotpath+'flux_ratio_sep')


f,ax = pp.paper_single_ax()
ax.scatter(pflux/rms, tflux/pflux, marker='.')
ax.loglog()
y1,y2 = ax.get_ylim()
ax.vlines(5, y1,y2, color='k')
ax.hlines(1.4, 0., 1e3, color='k')
pp.set_attrib(ax, xlabel='$S_{\mathrm{peak}}/\sigma_{\mathrm{local}}$', ylabel='$S_{\mathrm{total}}/S_{\mathrm{peak}}$',xlim=(4e-1,2e3))
pp.fig_save_many(f, plotpath+'flux_ratio_snr')

#sys.exit()

fmin = np.min(tflux)
fmax = np.max(tflux)

'''
# S&H: scale to RCB below 325MHz (avoids issues with CasA dec at low flux densities in usual Baars77)
# ratio of Baars 1977 to KPW @ 1400 MHz is 1.029

nu > 325: RCB and KPW consistent
 Roger, Costain & Bridle (1973)  - RCB
 Kellerman, Pauliny-Toth & Williams (1969) - KPW 
nu < 325: need to scale from KPW to RCB 

 B77 data needs to be scaled to RCB: using polynomial fit to ratios of Baars to KPW
ff = [38,  
81.5,
178,
400,
750,
1400,
2700,
5000,
10600,
15000]
frat = [0.981, 1.020, 1.051, 1.065, 1.046, 1.029, 1.011, 0.993, 0.974, 0.966]
# VLSS - B77
# VLSSr - RCB (S&H)
# wenss : average scaling to B77, then to KPW
# 6C already on RCB
# VLA-P on SH
# GMRT 610 SH
# NVSS is on B77
# WSRT1400 scaled to NVSS
# FIRST is on B77

'''

data_rows = [
             ('LOFAR60', '/net/beerze/data2/wwilliams/projects/bootes/Bootes_62MHz_vanWeeren.fits', 62., 'mJy', 60., 30., 30., 'RAJ2000','DEJ2000', 'S','e_S', 1e2, 1, 0.2),
             ('VLSSr', '/net/beerze/data2/wwilliams/projects/bootes/Bootes_VLSSr.fits', 74., 'Jy', 30., 75., 75., 'RA','DEC', 'Sp','e_Sp', 530., 1, 0.2),
             ('LOFAR150-DEEP', '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.blanked.scaled.cat.fits', 144., 'Jy', 30., 6., 6., 'RA','DEC', 'Total_flux','E_Total_flux', 10., 0.859, 0.1),
             ('LOFAR150', '/net/beerze/data2/wwilliams/projects/bootes/Bootes_150MHz_Williams2016.fits', 150., 'mJy', 30., 7.5, 15., 'RAJ2000','DEJ2000', 'Sint','e_Sint', 1e0, 1, 0.15),
             #('TGSS', '/net/beerze/data2/wwilliams/projects/ddf/TGSSADR1_7sigma_catalog.fits', 150., 'mJy', 30., 25., 25., 'RA','DEC', 'Total_flux','E_Total_flux', 65., 1, 0.2),
             ('6C', '/net/beerze/data2/wwilliams/projects/bootes/Bootes_6C_151MHz.fits', 151., 'Jy', 30., 60., 60., '_RAJ2000','_DEJ2000', 'F151p','F151p', 250., 1, 0.2),
             ('T-RaMiSu', '/net/beerze/data2/wwilliams/projects/bootes/Bootes_153MHz_Williams2013.fits', 153., 'mJy', 30., 25., 15., 'RAJ2000','DEJ2000', 'Si','e_Si', 10., 1, 0.15),
             ('VLA-P', '/net/beerze/data2/wwilliams/projects/bootes/Bootes_325MHz_Coppejans.fits', 325., 'mJy', 30., 5.6, 15., 'RAJ2000','DEJ2000', 'Si','e_Si', 2., 0.982/1.05, 0.15),
             ('WENSS', '/net/beerze/data2/wwilliams/projects/bootes/WENSS.cat.fits', 325., 'mJy', 30., 54., 30.,'_RAJ2000','_DEJ2000', 'Sint','Sint', 15., 0.982, 0.15),
             ('GMRT610', '/net/beerze/data2/wwilliams/projects/bootes/BOOTES610.MOSAIC.pybdsm.fits', 610., 'Jy', 30., 5.5, 5.5, 'RA','DEC', 'Total_flux','E_Total_flux', 2., 1, 0.15),
             ('RACS', '/net/beerze/data2/wwilliams/projects/radio_spectra/data/bootes_racs_new_clean.fits', 887.5, 'mJy', 30., 25., 25., 'ra_deg_cont','dec_deg_cont', 'flux_int','flux_int_err', 2., 1.041, 0.1),
             ('WSRT1400','/net/beerze/data2/wwilliams/projects/bootes/Bootes_1.4GHz_deVries.fits' , 1400., 'mJy', 60., 54., 54., '_RAJ2000','_DEJ2000', 'Si1.4GHz','e_Si1.4GHz', 1., 1.029, 0.1),
             ('NVSS', '/net/beerze/data2/wwilliams/projects/ddf/NVSS.fits', 1400., 'Jy', 60., 45., 45.,'RA','DEC', 'Total_flux','E_Total_flux', 2., 1.029, 0.1),
             ('FIRST', '/net/beerze/data2/wwilliams/projects/bootes/FIRST.cat.fits', 1400., 'mJy', 30., 5., 15.,'_RAJ2000','_DEJ2000', 'Fint','Fint', 2., 1.029, 0.1)
             
              ]
cat_info = Table(rows=data_rows, names=('SID', 'Cat_Name', 'Freq', 'FluxUnit', 'SizeLimit', 'resolution', 'mradius', 'racol','deccol', 'fluxcol','fluxerrcol', 'FluxLimit','scale_flux','scale_flux_err'))



Ccat = ac.SkyCoord(cat['RA'], cat['DEC'])
idx, sep,_ = Ccat.match_to_catalog_sky(Ccat,nthneighbor=2)
cat.add_column(Column(name='NNSep',data=sep.to('arcsec').value))


cat_info.add_column(Column(data=np.ones(len(cat_info)),name='FScale'))
#cat_info.add_column(Column(data=np.ones(len(cat_info)),name='FScale_Error'))
cat_info.add_column(Column(data=np.ones(len(cat_info)),name='E_u_FScale'))
cat_info.add_column(Column(data=np.ones(len(cat_info)),name='E_l_FScale'))
cat_info.add_column(Column(data=np.ones(len(cat_info)),name='E_FScale'))

fS = 1.


racol1='RA'
deccol1='DEC'
fluxcol1='Total_flux'
fluxerrcol1='E_Total_flux'

def bootstrap(data,function,iters):
    result=np.zeros(iters)
    for i in range(iters):
        x=data[np.random.randint(low=0,high=len(data),
                                 size=len(data))]
        result[i]=function(x)
    return result

for icat in range(len(cat_info)):
    

    sid = cat_info['SID'][icat]
    freq = cat_info['Freq'][icat]
    
    print(sid)
    
    mradius = max(cat_info['mradius'][icat], 15.)
    mradius_isol = max(cat_info['resolution'][icat], 15.)


    cat_name = cat_info['Cat_Name'][icat]
    mcat = Table.read(cat_name)
    
    fscale = cat_info['scale_flux'][icat]
    racol2 = cat_info['racol'][icat]
    deccol2 = cat_info['deccol'][icat]
    fluxcol2 = cat_info['fluxcol'][icat]
    fluxerrcol2 = cat_info['fluxerrcol'][icat]
    
    if mcat[racol2].unit =='DEGREE':
        mcat[racol2].unit='deg'
    if mcat[deccol2].unit =='DEGREE':
        mcat[deccol2].unit='deg'
    
    # calc nearest neighbour for mcat
    mCcat = ac.SkyCoord(mcat[racol2], mcat[deccol2], unit='deg')
    idx, sep,_ = mCcat.match_to_catalog_sky(mCcat,nthneighbor=2)
    mcat.add_column(Column(name='NNSep',data=sep.to('arcsec').value))
    
    
    # find nearest match for each cat source
    idx, sep,_ = Ccat.match_to_catalog_sky(mCcat)

    matchedcat = mcat[idx]
    match = (sep  < mradius*u.arcsec/2)
    
    

    print('{m:d} matches within {mr:.1f} arcsec'.format(m=np.sum(match),mr=mradius))
    
    #print (cols)
    
    # scale to common flux scale
    
    # hack the lofar60 column - its and int and masked
    flux2 = np.array(matchedcat[fluxcol2], dtype='float')*fscale
    fluxerr2 = np.array(matchedcat[fluxerrcol2], dtype='float')*fscale
    # WENSS/6C has 0 flux errors, assume 10%
    if sid in ['WENSS','6C']:
        fluxerr2 = 0.1*flux2
    flux2[flux2<-1e6] = np.nan
    fluxerr2[np.isnan(flux2)] = np.nan
    
    flux2[~match] = np.nan
    fluxerr2[~match] = np.nan
    
    
    if cat_info['FluxUnit'][icat] == 'mJy':
        flux2 = flux2 *1e-3
        fluxerr2 = fluxerr2 *1e-3
        
    print ('flux range', np.nanmin(flux2), np.nanmax(flux2))
    
    fscale54 = (54./cat_info['Freq'][icat]) **-0.78
    
    #print(flux2)
    
    isol = (cat['NNSep'] > mradius_isol/2)
    misol = (matchedcat['NNSep'] > mradius)
    
    print('{m:d} cat isolated within {mr:.1f} arcsec'.format(m=np.sum(isol&match),mr=mradius_isol))
    print('{m:d} match isolated within {mr:.1f} arcsec'.format(m=np.sum(misol&match),mr=mradius))
    
    sizelimt = (cat['Maj']*3600 < cat_info['SizeLimit'][icat])
    
    
    
    #'''
    #The threshold was set to
    #be the flux density that a source at the completeness limit of the
    #comparison survey multiplied by the LOFAR flux density that
    #such a source would have for a spectral index of 1.5 (except for
    #the lower frequency survey VLSSr, where a flat spectral index
    #was used). '''
    
    fxlim1 = cat_info['FluxLimit'][icat] *1e-3  # in Jy
    fylim1 =  fxlim1 * (54.393/freq)**-1.5
    X0 = np.log10(fxlim1)
    Y0 = np.log10(fylim1)
    M = -1.
    C = Y0 - M*X0
    
    
    
    thresh = (cat[fluxcol1] >= 10**(M*np.log10(flux2) + C  ))
    
    threshdeep = (flux2 > cat_info['FluxLimit'][icat])
    
    #goodmatch = match & isol & misol & thresh & sizelimt
    goodmatch = match & isol & misol & thresh 
    
    print(np.sum(goodmatch), ' good matches')
    
    ra2 = matchedcat[racol2]
    dec2 = matchedcat[deccol2]
    ra2[~match] = np.nan
    dec2[~match] = np.nan
    
    cat.add_column(Column(name=sid+'_RA', data=ra2))
    cat.add_column(Column(name=sid+'_DEC', data=dec2))
    cat.add_column(Column(name=sid+'_Flux', data=flux2))
    cat.add_column(Column(name=sid+'_Flux_Err', data=fluxerr2))
    cat.add_column(Column(name=sid+'_Flux_Scaled', data=flux2*fscale54))
    cat.add_column(Column(name=sid+'_Flux_Scaled_Err', data=fluxerr2*fscale54))
    cat.add_column(Column(name=sid+'_Matched', data=match))
    cat.add_column(Column(name=sid+'_Flag', data=goodmatch))
    
    
    
    #### plots per cat ###

    f,ax = pp.paper_single_ax()
    ax.scatter(3600*(cat[racol1]-matchedcat[racol2]), 3600*(cat[deccol1]-matchedcat[deccol2]))
    ax.scatter(3600*(cat[racol1]-matchedcat[racol2])[goodmatch], 3600*(cat[deccol1]-matchedcat[deccol2])[goodmatch])
    pp.invert_xlim(ax)
    pp.set_attrib(ax,xlabel='$\Delta$ RA (arcsec)', ylabel='$\Delta$ Dec (arcsec)', title=sid)
    pp.fig_save_many(f, plotpath+'{sid}_pos_offsets'.format(sid=sid))


    fmmin = np.nanmin(flux2)
    fmmax = np.nanmax(flux2)
        
    fx = np.linspace(fmmin,fmmax, 10)
    alp = np.array((0, -1.5))
    fysp = fx * ((54.393/freq)**alp[:,np.newaxis])
    
    
    fxthresh_r = np.linspace(fmmin,fmmax, 10)
    Xr = np.log10(fxthresh_r)
    Yr = M*Xr + C    
    fythresh_r = 10**Yr
    
    f,ax = pp.paper_single_ax()
    #ax.scatter(fS*tmlofar[fluxcol1],fS*flux2, marker='.', c='gray')
    ax.scatter(flux2, fS*cat[fluxcol1], marker='.', c='C0')
    ax.scatter(flux2[goodmatch], fS*cat[fluxcol1][goodmatch], marker='+', c='C1')
    y1,y2 = ax.get_ylim()
    #ax.scatter(newflux2[~good],fS*newflux1[~good], marker='+', c='C1')
    ax.plot(fx, fysp.T, 'k', alpha=0.5)
    ax.plot(fxthresh_r, fythresh_r, 'C4', alpha=0.5)
    ax.loglog()
    ax.vlines(cat_info['FluxLimit'][icat],y1,y2)
    pp.set_attrib(ax, ylabel='$S_{54\mathrm{\,MHz}}$ (Jy)',xlabel='$S_{{{f:.0f}\mathrm{{\,MHz}}}}$ (Jy)'.format(f=freq), title=sid)
    pp.fig_save_many(f, plotpath+'{sid}_fluxes1'.format(sid=sid))
    
    frac = flux2/cat[fluxcol1]
    
    print(np.median(frac[goodmatch]), np.std(frac[goodmatch]))
    
    med = np.median(frac[goodmatch])
    err1, err2 = np.percentile(bootstrap(frac[goodmatch],np.median,1000),(16,84))

    EFScale = cat_info['scale_flux_err'][icat]
    
    fscale = fscale54
    cat_info['FScale'][icat] = med
    cat_info['E_u_FScale'][icat] =np.sqrt( (err2-med)**2. +  (EFScale*med)**2.)
    cat_info['E_l_FScale'][icat] =np.sqrt( (med-err1)**2. +  (EFScale*med)**2.)
    cat_info['E_FScale'][icat] =np.sqrt(np.max((med-err1,err2-med ))**2.+ (EFScale*med)**2.)
    
    f,ax = pp.paper_single_ax()
    #ax.scatter(fS*tmlofar[fluxcol1],fS*flux2, marker='.', c='gray')
    ax.scatter(fS*cat[fluxcol1],fscale*flux2, marker='.', c='C0')
    ax.scatter(fS*cat[fluxcol1][goodmatch],fscale*flux2[goodmatch], marker='+', c='C1')
    #ax.plot(fx, fysp.T, 'k', alpha=0.5)
    ax.plot(fx, fx, 'k', alpha=0.5)
    ax.loglog()
    pp.set_attrib(ax, xlabel='$S_{54\mathrm{\,MHz}}$ (Jy)',ylabel='Scaled $S_{{{f:.0f}\mathrm{{\,MHz}}}}$ (Jy)'.format(f=freq), title=sid)
    pp.fig_save_many(f, plotpath+'{sid}_fluxes'.format(sid=sid))


    f,ax = pp.paper_single_ax()
    #ax.scatter(fS*tmlofar[fluxcol1],flux2/tmlofar[fluxcol1], marker='.', c='gray')
    ax.scatter(fS*cat[fluxcol1],fscale*flux2/cat[fluxcol1], marker='.', c='C0')
    ax.scatter(fS*cat[fluxcol1][goodmatch],fscale*flux2[goodmatch]/cat[fluxcol1][goodmatch], marker='+', c='C1')
    ax.plot(fx, np.ones_like(fx), 'k', alpha=0.5)
    ax.plot(fx, np.ones_like(fx)*np.median(fscale*flux2[goodmatch]/cat[fluxcol1][goodmatch]), 'k', alpha=0.5,linestyle='dashed')
    ax.loglog()
    pp.set_attrib(ax, xlabel='$S_{54\mathrm{\,MHz}}$ (Jy)',ylabel='scaled $S_{{{f:.0f}\mathrm{{\,MHz}}}}/S_{{{f2:.0f}\mathrm{{\,MHz}}}}$'.format(f=freq,f2=54.), title=sid)
    pp.fig_save_many(f, plotpath+'{sid}_flux_ratios'.format(sid=sid))
    
    

    #goodmatch = match & isol & misol & thresh 
    mm = match & isol & misol & thresh 
    f,ax = pp.paper_single_ax()
    xx = fscale*flux2[goodmatch]/cat[fluxcol1][goodmatch]
    v1,v2 = np.percentile(xx[np.isfinite(xx)],[12,84])
    #ax.scatter(fS*tmlofar[fluxcol1],flux2/tmlofar[fluxcol1], marker='.', c='gray')
    #c=ax.scatter(cat[racol1][goodmatch], cat[deccol1][goodmatch],c=fscale*flux2[goodmatch]/cat[fluxcol1][goodmatch], vmin=v1,vmax=v2, marker='o')
    c=ax.scatter(cat[racol1][mm], cat[deccol1][mm],c=fscale*flux2[mm]/cat[fluxcol1][mm], vmin=v1,vmax=v2, marker='o')
    cbar=plt.colorbar(c)
    cbar.set_label('Fratio')
    pp.set_attrib(ax, xlabel='RA',ylabel='DEC', title=sid)
    pp.fig_save_many(f, plotpath+'{sid}_pos_flux_ratios'.format(sid=sid))
    
    #mm = match   & isol  &thresh
    #f,ax = pp.paper_single_ax()
    #xx = fscale*flux2[goodmatch]/cat[fluxcol1][goodmatch]
    #v1,v2 = np.percentile(xx[np.isfinite(xx)],[12,84])
    ##ax.scatter(fS*tmlofar[fluxcol1],flux2/tmlofar[fluxcol1], marker='.', c='gray')
    ##c=ax.scatter(cat[racol1][goodmatch], cat[deccol1][goodmatch],c=fscale*flux2[goodmatch]/cat[fluxcol1][goodmatch], vmin=v1,vmax=v2, marker='o')
    #c=ax.scatter(mcat[racol2], mcat[deccol2],c=mcat['NNSep'], marker='o')
    #cbar=plt.colorbar(c)
    #cbar.set_label('nnsep')
    #pp.set_attrib(ax, xlabel='RA',ylabel='DEC', title=sid)
    #pp.fig_save_many(f, plotpath+'{sid}_pos_nnsep'.format(sid=sid))
    
    #sys.exit()
    plt.close('all')
    
    #/S_{54\mathrm{\,MHz}}
#cat['WSRT1400_Flux'] = cat['WSRT1400_Flux']/1000  ### to Jy from mJy
    


def fit_leastsq(p0, datax, datay, datayerr, function):

    errfunc = lambda p, x, y, yerr: (function(p,x) - y)/(yerr**2.)

    pfit, pcov, infodict, errmsg, success = \
        leastsq(errfunc, p0, args=(datax, datay, datayerr), \
                          full_output=1, epsfcn=0.0001)

    if (len(datay) > len(p0)) and pcov is not None:
        s_sq = (errfunc(pfit, datax, datay, datayerr)**2).sum()/(len(datay)-len(p0))
        pcov = pcov * s_sq
    else:
        pcov = np.inf
    if (len(datay) <= len(p0)):
        pfit = np.nan*p0

    error = [] 
    for i in range(len(pfit)):
        try:
          error.append(np.absolute(pcov[i][i])**0.5)
        except:
          error.append( 0.00 )
    pfit_leastsq = pfit
    perr_leastsq = np.array(error) 
    return pfit_leastsq, perr_leastsq 


#print("\n# Fit parameters and parameter errors from lestsq method :")
#print("pfit = ", pfit)
#print("perr = ", perr)

funcLine = lambda tpl,x : tpl[0]*x + tpl[1]
ErrorFunc = lambda tpl,x,y,ye: (funcLine(tpl,x)-y)/(ye**2.)

funcLine2 = lambda tpl,x : tpl[0]*x**2. + tpl[1]*x + tpl[2]
ErrorFunc2 = lambda tpl,x,y,ye: (funcLine2(tpl,x)-y)/(ye**2.)

def booterrs(func, p, perr, x, nit=100):
    
    #y = func(p,x)
    ys = np.nan*np.ones(nit)
    for i in range(nit):
        pp = [ np.random.normal(pi,pei) for pi,pei in zip(p, perr)]
        ys[i] = func(pp,x)
    
    return np.nanmean(ys), np.nanstd(ys)

def fit_alp(freqs, fluxes,fluxerrs,freq0=53., freqrange=None):
    freqs = np.array(freqs)
    fluxes = np.array(fluxes)
    fluxerrs = np.array(fluxerrs)
    
    p0 = np.array((-0.7, np.log10(np.nanmean(fluxes))))
    
    # leastsq finds the set of parameters in the tuple tpl that minimizes
    # ErrorFunc=yfit-yExperimental
    yerr = np.abs(np.log(10.) * fluxerrs/fluxes)
    x = np.log10(freqs)
    y = np.log10(fluxes)
    
    # fit_leastsq takes into account the yerrs
    p , p_err = fit_leastsq(p0, x, y, yerr, funcLine)
    
    #S53, S53e =  p[1], p_err[1]
    
    S53, S53e  = booterrs(funcLine, p, p_err, np.log10(freq0))
    S53e = np.abs( (10**S53)*np.log(10) * S53e)
    S53 = 10**S53
    
    if freqrange is None:
        return S53, S53e
    else:
        
        xr = np.log10(freqrange)
        lSrange = funcLine(p, xr)
        Srange = 10**lSrange
        return S53, S53e, Srange




def fit_alp2(freqs, fluxes,fluxerrs,freq0=53., freqrange=None):
    freqs = np.array(freqs)
    fluxes = np.array(fluxes)
    fluxerrs = np.array(fluxerrs)
    
    p20 = np.array((0,-0.7, np.log10(np.nanmean(fluxes))))
    
    yerr = np.abs(np.log(10.) * fluxerrs/fluxes)
    x = np.log10(freqs)
    y = np.log10(fluxes)
    
    p2, p2_err = fit_leastsq(p20, x, y, yerr, funcLine2)
    
    S53_2, S53e_2  = booterrs(funcLine2, p2, p2_err, np.log10(freq0))
    #S53_2, S53e_2 =  p2[2], p2_err[2]
    
    S53e_2= np.abs( (10**S53_2)*np.log(10) * S53e_2)
    S53_2 = 10**S53_2
    if freqrange is None:
        return S53_2, S53e_2
    else:
        xr = np.log10(freqrange)
        lSrange = funcLine2(p2, xr)
        Srange = 10**lSrange
        return S53_2, S53e_2, Srange


ffrange = np.logspace(np.log10(50), np.log10(1500),100)
S53, S53e, Srange = fit_alp(cat_info['Freq'], cat_info['FScale'],cat_info['E_FScale'],freq0=53., freqrange=ffrange)
print ('Flux scaling {:5.3f} +/- {:5.2f}'.format(S53, S53e))
norm = mpl.colors.Normalize(vmin=0, vmax=len(cat_info))
cmap = mpl.cm.get_cmap('viridis')
f,ax = pp.paper_single_ax(TW=8, AR=0.6)
ax.plot(54., 1., c='C0',marker='o')
for icat,cc in enumerate(cat_info):
    c = cmap(norm(icat))
    ax.errorbar(np.array([cc['Freq']]),np.array([cc['FScale']]),np.array([[cc['E_u_FScale'],cc['E_l_FScale']]]).T, marker='x',color=c, label=cc['SID'],capsize=5)
    #ax.errorbar(cc['Freq'],cc['FScale'],cc['E_FScale'], marker='+',color=c, label=cc['SID'])
pp.set_attrib(ax, xlabel='Frequency (MHz)', ylabel='Flux density ratio')
ax.loglog()
#ax.legend()
ax.legend(bbox_to_anchor=(1.0,1.05), loc="upper left")
plt.subplots_adjust(left=0.1,right=0.725)
ax.plot(ffrange, Srange, color='r', ls='dashed')
ax.set_ylim(0.05, 2)
ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
pp.fig_save_many(f, plotpath+'flux_density_scale')

cat_info['FScale'].format = '%.3f' 
cat_info['E_FScale'].format = '%.3f'  
cat_info['E_l_FScale'].format = '%.3f'  
cat_info['E_u_FScale'].format = '%.3f'  

with open('flux_density_scale.log','w') as ff:
    ff.write('Flux scaling {:5.3f} +/- {:5.2f}'.format(S53, S53e))
    ff.write(str(cat_info['SID','Freq',  'SizeLimit', 'resolution', 'mradius', 'FluxLimit', 'scale_flux', 'FScale', 'E_u_FScale', 'E_l_FScale', 'E_FScale']) )
print(cat_info['SID','Freq',  'SizeLimit', 'resolution', 'mradius', 'FluxLimit', 'scale_flux', 'FScale', 'E_u_FScale', 'E_l_FScale', 'E_FScale']) 

#sids = ['LOFAR60', 'VLSS', 'LOFAR150', 'TGSS', 'VLAP','WENSS', 'WSRT1400','NVSS']

#tl = [s+'_Flag' for s in sids]

def plot_spec(i):
    f,ax = pp.paper_single_ax()
    ax.loglog()
    #ax.scatter(53., cat['Total_flux'][i], c='k')
    ax.errorbar(53., cat['Total_flux'][i], cat['E_Total_flux'][i], c='k')
    for f,s in zip(freqs,sids):
        #ax.scatter(f, cat[s+'_Flux'][i])
        ax.errorbar(f, cat[s+'_Flux'][i], cat[s+'_Flux_Err'][i])
    return
#mask3 = np.sum(cat[tl],axis=1) > 3

# for using only VLSS, WENSS and NVSS
mask3 = cat['VLSSr_Flag'] & cat['WENSS_Flag'] & cat['NVSS_Flag']
mask60 = cat['LOFAR60_Flag']


#cat = cat[mask3]
alp = np.nan*np.ones(len(cat))
s53 = np.nan*np.ones(len(cat))
s53_2 = np.nan*np.ones(len(cat))
s53_full = np.nan*np.ones(len(cat))
s53_2_full = np.nan*np.ones(len(cat))
s53e = np.nan*np.ones(len(cat))
s53e_2 = np.nan*np.ones(len(cat))
s53e_full = np.nan*np.ones(len(cat))
s53e_2_full = np.nan*np.ones(len(cat))
nb = np.nan*np.ones(len(cat))


freqs = cat_info['Freq']
sids = cat_info['SID']

ffrange = np.logspace(np.log10(50), np.log10(1500),100)
for i in range(len(cat)):
    plot = False
    
    ff = [ f for f, s in zip(freqs, sids) if cat[s+'_Flag'][i]]
    fF = np.array([ cat[s+'_Flux'][i] for f, s in zip(freqs, sids) if cat[s+'_Flag'][i]])
    fFe = np.array([ cat[s+'_Flux_Err'][i] for f, s in zip(freqs, sids) if cat[s+'_Flag'][i]])
    #ff = [ f for f, s in zip(freqs, sids) ]#if cat[s+'_Flag'][i]]
    #fF = np.array([ cat[s+'_Flux'][i] for f, s in zip(freqs, sids) ])#if cat[s+'_Flag'][i]])
    #fFe = np.array([ cat[s+'_Flux_Err'][i] for f, s in zip(freqs, sids) ])#if cat[s+'_Flag'][i]])
    
    fFe = np.sqrt(fFe**2+(0.1*fF)**2.)
    
    #if len(ff) <= 2:
        #continue
    
    if not os.path.isfile(plotpath+'seds/sed{i:04d}.png'.format(i=cat['Source_id'][i])):
        plot = True
    
    if plot:
        f = plt.figure(figsize=(12,6))
        ax2 = plt.subplot(121)
        ax = plt.subplot(122, sharex=ax2, sharey=ax2)
        #f,ax2 = pp.paper_single_ax()
        ax.loglog()
        ax2.loglog()
        ax2.set_title('all')
    
    
    
    nb[i] = len(ff)
    if len(ff) >= 3:
        s53_full[i], s53e_full[i], Sr = fit_alp(ff, fF,fFe, freqrange=ffrange)
        if plot:
            ax2.errorbar(ff, fF, fFe, c='C0', marker='o',ms=5, ls='')
            ax2.errorbar([75,324, 1400], [cat['VLSSr_Flux'][i], cat['WENSS_Flux'][i],cat['NVSS_Flux'][i]], [cat['VLSSr_Flux_Err'][i], cat['WENSS_Flux_Err'][i],cat['NVSS_Flux_Err'][i]], marker='o',c='C0',ms=10, ls='')
        
            ax2.errorbar(53., cat['Total_flux'][i], cat['E_Total_flux'][i], c='C1', marker='o', label='observed')
        
            ax2.errorbar(53., s53_full[i], s53e_full[i], c='C2', marker='o', label='linear')
            ax2.plot(ffrange, Sr, color='C2')
        
    if len(ff) >= 4:
        s53_2_full[i], s53e_2_full[i], Sr2 = fit_alp2(ff, fF,fFe, freqrange=ffrange)
        if plot:
            ax2.errorbar(53., s53_2_full[i], s53e_2_full[i], c='C3', marker='o', label='2nd order')
            ax2.plot(ffrange, Sr2, color='C3')
        #print (s53_full[i], s53_2_full[i], ff,fF)
    if plot:
        ax2.legend()
        
        
    if mask3[i]:
        s53[i], s53e[i], Sr = fit_alp([75,324, 1400], [cat['VLSSr_Flux'][i], cat['WENSS_Flux'][i],cat['NVSS_Flux'][i]], [cat['VLSSr_Flux_Err'][i], cat['WENSS_Flux_Err'][i],cat['NVSS_Flux_Err'][i]], freqrange=ffrange)
        s53_2[i], s53e_2[i], Sr2 = fit_alp2([75,324, 1400], [cat['VLSSr_Flux'][i], cat['WENSS_Flux'][i],cat['NVSS_Flux'][i]], [cat['VLSSr_Flux_Err'][i], cat['WENSS_Flux_Err'][i],cat['NVSS_Flux_Err'][i]], freqrange=ffrange)
        
        if plot:
            ax.errorbar([75,324, 1400], [cat['VLSSr_Flux'][i], cat['WENSS_Flux'][i],cat['NVSS_Flux'][i]], [cat['VLSSr_Flux_Err'][i], cat['WENSS_Flux_Err'][i],cat['NVSS_Flux_Err'][i]], marker='o',c='C0',ms=5, ls='')
            ax.errorbar(53., cat['Total_flux'][i], cat['E_Total_flux'][i], c='C1', marker='o', label='observed')
            ax.errorbar(53., s53[i], s53e[i], c='C2', marker='o', label='linear')
            ax.errorbar(53., s53_2[i], s53e_2[i], c='C3', marker='o', label='2nd order')
            ax.plot(ffrange, Sr, color='C2')
            ax.plot(ffrange, Sr2, color='C3')
            ax.set_title('VLSS+WENSS+NVSS')
            ax.legend()
    
    
        
    if plot:
        pp.set_attrib(ax2, xlabel=r'$\nu$ (MHz)', ylabel=r'$S_{\nu}$ (Jy)')
        pp.set_attrib(ax, xlabel=r'$\nu$ (MHz)')
        pp.fig_save_many(f, plotpath+'seds/sed{i:04d}'.format(i=cat['Source_id'][i]))
    
        plt.close(f)
    
cat.add_column(MaskedColumn(data=cat['Total_flux']/s53_2_full,name='S53_pred_2_f'))
cat.add_column(MaskedColumn(data=cat['Total_flux']/s53_full,name='S53_pred_1_f'))
    
f,ax = pp.paper_single_ax()
ax.errorbar(cat['Total_flux'][mask3], s53[mask3], s53e[mask3], cat['E_Total_flux'][mask3],fmt='.')
c = ax.scatter(cat['Total_flux'][mask3], s53[mask3])
#cbar = plt.colorbar(c)
#cbar.set_label('$\\alpha$')
ax.plot([0.5,20],[0.5,20],'k')
ax.loglog()
pp.set_attrib(ax, xlabel='$S_{53,obs}$ (Jy)', ylabel='$S_{53,pred}$ (Jy)')
pp.fig_save_many(f, plotpath+'pred_vs_meas_flux_ratios_nvss_wenss_vlss_1st')



s53_2_full[s53_2_full>20] = np.nan ## this is a bad fit??
s53_2_full[s53_2_full<1e-3] = np.nan ## this is a bad fit??
s53_full[s53_full>20] = np.nan ## this is a bad fit??
s53_full[s53_full<1e-3] = np.nan ## this is a bad fit??
s53_full[np.abs(s53e_full/s53_full)>0.1] = np.nan ## this is a bad fit??
s53_2_full[np.abs(s53e_2_full/s53_2_full)>0.1] = np.nan ## this is a bad fit??


f,ax = pp.paper_single_ax()
ax.errorbar(cat['Total_flux'][nb>3], s53_2_full[nb>3], s53e_2_full[nb>3], cat['E_Total_flux'][nb>3],fmt='.',label='$n=3$ (linear)')
ax.errorbar(cat['Total_flux'][nb==3], s53_full[nb==3], s53e_full[nb==3], cat['E_Total_flux'][nb==3],fmt='.',label='$n>3$ (2nd order)')
ax.errorbar(cat['Total_flux'][mask60], s53_full[mask60], s53e_full[mask60], cat['E_Total_flux'][mask60],fmt='.',label='60MHz')
#c = ax.scatter(cat['Total_flux'][nb>3], s53_full[nb>3],c=nb[nb>3],vmin=0,vmax=6)
#c = ax.scatter(cat['Total_flux'][nb==3], s53_full[nb==3],c=nb[nb==3],vmin=0,vmax=6)
#cbar = plt.colorbar(c)
#cbar.set_label('$\\alpha$')
ax.plot([5e-3,20],[5e-3,20],'k')
ax.loglog()
ax.legend()
pp.set_attrib(ax, xlabel='$S_{53,obs}$ (Jy)', ylabel='$S_{53,pred}$ (Jy)')
pp.fig_save_many(f, plotpath+'pred_vs_meas_flux_ratios')

s53_pred = s53_2_full.copy()
s53_pred[nb==3] = s53_full[nb==3]
s53_pred = cat['Total_flux']/s53_pred
f,ax = pp.paper_single_ax()
ax.hist(s53_pred,bins=20,histtype='step')
m_s53_pred = np.nanmean(s53_pred)
s_s53_pred = np.nanstd(s53_pred)
y1,y2 = ax.get_ylim()
ax.vlines(m_s53_pred, y1,y2,color='k')
ax.vlines(m_s53_pred+s_s53_pred, y1,y2,color='k')
ax.vlines(m_s53_pred-s_s53_pred, y1,y2,color='k')
pp.fig_save_many(f, plotpath+'pred_vs_meas_flux_ratios_hist')

mask_full = nb >=3

print('meas/pred(3) full {:.2f} \pm {:.2f} '.format(m_s53_pred, s_s53_pred))
print('meas/pred(2) full {:.2f} \pm {:.2f} '.format(np.mean((cat['Total_flux']/s53_full)[mask_full]), np.std((cat['Total_flux']/s53_full)[mask_full])))
print('meas/pred(2) full {:.2f} \pm {:.2f} '.format(np.mean((cat['Total_flux']/s53_2_full)[mask_full]), np.std((cat['Total_flux']/s53_2_full)[mask_full])))
print('VLSS/WENSS/NVSS only')
print('meas/pred(2) {:.2f} \pm {:.2f} '.format(np.mean((cat['Total_flux']/s53)[mask3]), np.std((cat['Total_flux']/s53)[mask3])))
print('meas/pred(3) {:.2f} \pm {:.2f} '.format(np.mean((cat['Total_flux']/s53_2)[mask3]), np.std((cat['Total_flux']/s53_2)[mask3])))



catvlap = Table.read(cat_info['Cat_Name'][cat_info['SID'] == 'VLA-P'][0])
catwenss = Table.read(cat_info['Cat_Name'][cat_info['SID'] == 'WENSS'][0])
Cvlap = ac.SkyCoord(catvlap[cat_info['racol'][cat_info['SID'] == 'VLA-P'][0]], catvlap[cat_info['deccol'][cat_info['SID'] == 'VLA-P'][0]])
Cwenss = ac.SkyCoord(catwenss[cat_info['racol'][cat_info['SID'] == 'WENSS'][0]], catwenss[cat_info['deccol'][cat_info['SID'] == 'WENSS'][0]])
idx, sep,_ = Cvlap.match_to_catalog_sky(Cvlap,nthneighbor=2)
misol = (sep > 54*u.arcsec) & (catvlap[cat_info['fluxcol'][cat_info['SID'] == 'VLA-P'][0]]>0.2)
idx, sep,_ = Cvlap.match_to_catalog_sky(Cwenss)
mm = (sep <10*u.arcsec)
catvlap = catvlap[mm&misol]
catwenss = catwenss[idx[mm&misol]]


plt.figure()
plt.scatter(cat['WENSS_Flux'],cat['VLA-P_Flux'])
plt.scatter(1e-3*catwenss[cat_info['fluxcol'][cat_info['SID'] == 'WENSS'][0]],1e-3*catvlap[cat_info['fluxcol'][cat_info['SID'] == 'VLA-P'][0]])
plt.plot([0,0.5],[0,0.5],c='k')
plt.plot([0,0.5],[0,0.5],c='k')
plt.loglog()
plt.xlabel('WENSS')
plt.ylabel('VLAP')
plt.xlim(0.01,1)
plt.ylim(0.001,1)
plt.savefig(plotpath+'vlap_vs_wenss.png')
print('median WENSS/VLAP = {f:.2f}'.format(f=np.nanmedian(cat['WENSS_Flux']/cat['VLA-P_Flux'])))


plt.figure()
plt.scatter(catvlap[cat_info['fluxcol'][cat_info['SID'] == 'VLA-P'][0]],catwenss[cat_info['fluxcol'][cat_info['SID'] == 'WENSS'][0]]/catvlap[cat_info['fluxcol'][cat_info['SID'] == 'VLA-P'][0]])
plt.semilogx()

plt.figure()
#plt.scatter(cat['WENSS_Flux'],cat['VLA-P_Flux'])
#plt.scatter(Cvlap.ra.value,Cvlap.dec.value,marker='+')
#plt.scatter(Cwenss.ra.value,Cwenss.dec.value,marker='x')
c=plt.scatter(catvlap[cat_info['racol'][cat_info['SID'] == 'VLA-P'][0]],catvlap[cat_info['deccol'][cat_info['SID'] == 'VLA-P'][0]],c=catwenss[cat_info['fluxcol'][cat_info['SID'] == 'WENSS'][0]]/catvlap[cat_info['fluxcol'][cat_info['SID'] == 'VLA-P'][0]], vmin=0.5,vmax=2)
plt.colorbar(c)




plt.figure()
plt.scatter(cat['NVSS_Flux'],cat['WSRT1400_Flux'])
plt.plot([0,0.5],[0,0.5],c='k')
plt.plot([0,0.5],[0,0.5],c='k')
plt.loglog()
plt.xlabel('NVSS')
plt.ylabel('WSRT')
plt.xlim(0.001,0.1)
plt.ylim(0.0001,0.1)
plt.savefig(plotpath+'wsrt_vs_nvss.png')
print('median NVSS/WSRT = {f:.2f}'.format(f=np.nanmedian(cat['NVSS_Flux']/cat['WSRT1400_Flux'])))




plt.figure()
plt.scatter(cat['Total_flux'],cat['VLA-P_Flux_Scaled'],alpha=0.5)
plt.plot([0,10],[0,10],c='k')
plt.plot([0,10],[0,10],c='k')
plt.loglog()
plt.xlabel('54 MHZ flux')
plt.ylabel('VLAP scaled')
plt.xlim(0.001,10)
plt.ylim(0.001,10)
plt.savefig(plotpath+'vlap_vs_lba.png')
plt.figure()
plt.scatter(cat['Total_flux'],cat['WENSS_Flux_Scaled'],alpha=0.5)
plt.plot([0,10],[0,10],c='k')
plt.plot([0,10],[0,10],c='k')
plt.loglog()
plt.xlabel('54 MHZ flux')
plt.ylabel('WENSS scaled')
plt.xlim(0.001,10)
plt.ylim(0.001,10)
plt.savefig(plotpath+'wenss_vs_lba.png')
