import os
import utils.plotting as pp
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table, Column, MaskedColumn, join
import astropy.units as u
import pyregion

#clobber = True
clobber = 1
stilts = '/net/beerze/data2/wwilliams/software/bin/stilts '

#lofar_cat_name = 'bootes_deep_lba.cat.fits'
lofar_cat_name = 'bootes_deep_lba_hbashift.cat.fits'

facet_regfile= 'facets_fin.reg'

## match to deep opt

#lofar_code_name = 'lba_match_deepwsrt_codes.fits'
             #('WSRT1400','/net/beerze/data2/wwilliams/projects/bootes/Bootes_1.4GHz_deVries.fits' , 1400., 'mJy', 60., 54., 54., '_RAJ2000','_DEJ2000', 'Si1.4GHz','e_Si1.4GHz', 1., 0.1),


m_lofar_cat_name = 'bootes_deep_lba.cat.match_{n}.fits'.format(n='deepwsrt')
m_lofar_cat_name = 'bootes_deep_lba_hbashift.cat.match_{n}.fits'.format(n='WSRT1400')
cat_name = '/net/beerze/data2/wwilliams/projects/bootes/Bootes_1.4GHz_deVries.fits'
mradius = 10.
if (not os.path.exists(m_lofar_cat_name)) or clobber:
    cmd = stilts+'tskymatch2 in1={in1} ra1=RA dec1=DEC ra2=_RAJ2000 dec2=_DEJ2000 in2={in2} error={rad:f} out={out} find=all'.format(in1=lofar_cat_name, in2=cat_name, out=m_lofar_cat_name, rad=mradius)
    os.system(cmd)
    
    
lcat = Table.read(lofar_cat_name)
#lcat_codes = Table.read(lofar_code_name)
hcat = Table.read(cat_name)
    
cat = Table.read(m_lofar_cat_name)
cat.sort('Source_id')

#cat = join(cat , lcat_codes)

#trms['Area']
#trms['rms']

# select only the 'proper' detected sources.... what to do about completeness of higher order wavelet scales??
print('Catalogue contains {0} sources'.format(len(cat)))
#sel = (cat['match_code']==1)
#cat = cat[sel]
#print('Catalogue contains {0} sources after clean selection'.format(len(cat)))
sel = (cat['Peak_flux'] / cat['Isl_rms']) > 5.
cat = cat[sel]
print('Catalogue contains {0} sources after S/N selection'.format(len(cat)))
sel = cat['GroupSize'].mask
cat = cat[sel]
print('Catalogue contains {0} sources after group==1 selection'.format(len(cat)))

tflux = cat['Total_flux'] 
pflux = cat['Peak_flux']
rms = cat['Isl_rms']


#maskS = (cat['S_Code'] == 'S')  & (cat['Maj'] < 30.)
#cat = cat[maskS]
#print('Catalogue contains {0} sources after S, size selection'.format(len(cat)))

#maskS = ~(cat['ID_OPTICAL'].mask)
#cat = cat[maskS]
#print('Catalogue contains {0} sources after optical selection'.format(len(cat)))

sel = ((cat['Peak_flux'] / cat['Isl_rms']) > 10.) & (cat['Maj']*3600 < 30.)
sel = ((cat['Peak_flux'] / cat['Isl_rms']) > 7.5) & (cat['Maj']*3600 < 25.)& (cat['S_Code']=='S') 
print('Catalogue contains {0} sources with high snr and small size selection'.format(len(cat[sel])))

N = len(cat)

dRA = 3600.*(cat['RA'] - cat['_RAJ2000']) *np.cos(cat['DEC']*np.pi/180.)
dDEC = 3600.*(cat['DEC'] - cat['_DEJ2000'])
d = np.sqrt(dRA**2. + dDEC**2.)

def stats(x):
    N = len(x)
    m = np.mean(x)
    s = np.std(x)
    me = s/np.sqrt(N)
    return m,s,me

dRAm, dRAs, dRAms = stats(dRA)
dDECm, dDECs, dDECms = stats(dDEC)
dm, ds, dms = stats(d)


print('dRA = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dRAm, me=dRAms, s=dRAs))
print('dDEC = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dDECm, me=dDECms, s=dDECs))
print('d = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dm, me=dms, s=ds))

dRAm, dRAs, dRAms = stats(dRA[sel])
dDECm, dDECs, dDECms = stats(dDEC[sel])
dm, ds, dms = stats(d[sel])

print('after clean selection')
print('dRA = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dRAm, me=dRAms, s=dRAs))
print('dDEC = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dDECm, me=dDECms, s=dDECs))
print('d = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dm, me=dms, s=ds))

f,ax = pp.paper_single_ax()
ax.hist(cat['Separation'], histtype='step')

f,ax = pp.paper_single_ax()
ax.hist(dRA, bins=25, histtype='step', label='dRA')
ax.hist(dDEC, bins=25, histtype='step', label='dDEC')
ax.hist(d, bins=25, histtype='step', label='d')
pp.set_attrib(ax,xlabel='separation (arcsec)', ylabel='N')

#f,ax = pp.paper_single_ax()
#c = ax.scatter(hcat['RA'], hcat['DEC'], c='gray', marker='.')
#c = ax.scatter(lcat['RA'], lcat['DEC'], c='k', marker='.')
#c = ax.scatter(cat['RA'], cat['DEC'], c='C0', marker='.')
##cbar = plt.colorbar(c)
##cbar.set_label('dRA (arcsec)')
##ax.hlines(1.4, 0.1, 1e3, color='k')
#pp.set_attrib(ax, xlabel='RA (J2000)', ylabel='DEC (J2000)')
#pp.fig_save_many(f, 'beerze_plots/deepwsrt_hbashift_astro_qual_offset_dRA_pos')

def ell(x0,y0,xs,ys):
    th = np.arange(0, 2*np.pi+0.001, 0.001)
    x = x0 + xs*np.cos(th)
    y = y0 + ys*np.sin(th)
    return x,y


f,ax = pp.paper_single_ax()
ax.scatter(dRA[sel], dDEC[sel], marker='.', c='C0', alpha=0.5)
#ax.scatter(dRA[sel], dDEC[sel], marker='.', c='C1', alpha=0.5)
#ax.hlines(1.4, 0.1, 1e3, color='k')
x1,x2 = ax.get_xlim()
y1,y2 = ax.get_ylim()
x,y=ell(dRAm,dDECm,dRAs,dDECs)
ax.plot(x,y,'k',alpha=0.8,linestyle='dashed')
ax.hlines(dDECm,-11,11,'k',alpha=0.8,linestyle='dashed')
ax.vlines(dRAm,-10,10,'k',alpha=0.8,linestyle='dashed')
pp.set_attrib(ax, xlabel='dRA (arcsec)', ylabel='dDEC (arcsec)')
ax.axis('square') 
ax.set_xlim(-10,10) 
ax.set_ylim(-9,9)  
#ax.set_ylim(-8,8) 
pp.fig_save_many(f, 'beerze_plots/deepwsrt_hbashift_astro_qual_offset')

r2 = pyregion.open(facet_regfile)
patch_list, artist_list = r2.get_mpl_patches_texts(origin=0)

f,ax = pp.paper_single_ax()
for p in patch_list: ax.add_patch(p)
#c = ax.scatter(lcat['RA'], lcat['DEC'], c='gray', marker='.', s=5)
c = ax.scatter(cat['RA'][sel], cat['DEC'][sel], c=dRA[sel], marker='.', vmin=-2*dRAs, vmax=2*dRAs)
cbar = plt.colorbar(c, extend='both')
cbar.set_label('dRA (arcsec)')
#ax.hlines(1.4, 0.1, 1e3, color='k')
pp.invert_xlim(ax)
pp.set_attrib(ax, xlabel='RA (J2000)', ylabel='DEC (J2000)',xlim=[221.3,214.8], ylim=[31.6,37.06])
pp.fig_save_many(f, 'beerze_plots/deepwsrt_hbashift_astro_qual_offset_dRA_pos')

patch_list, artist_list = r2.get_mpl_patches_texts(origin=0)
f,ax = pp.paper_single_ax()
for p in patch_list: ax.add_patch(p)
c = ax.scatter(cat['RA'][sel], cat['DEC'][sel], c=dDEC[sel], marker='.', vmin=-2*dDECs, vmax=2*dDECs)
cbar = plt.colorbar(c, extend='both')
cbar.set_label('dDEC (arcsec)')
#ax.hlines(1.4, 0.1, 1e3, color='k')
pp.invert_xlim(ax)
pp.set_attrib(ax, xlabel='RA (J2000)', ylabel='DEC (J2000)',xlim=[221.3,214.8], ylim=[31.6,37.06])
pp.fig_save_many(f, 'beerze_plots/deepwsrt_hbashift_astro_qual_offset_dDEC_pos')

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

#cat = cat[sel]
cat['Si1.4GHz'] = cat['Si1.4GHz']*1e-3

nu1 = 53.
nu2 = 1400.
f1min = 0.000646 * 5  # min rms
f2min = 0.001116 * 5  # median rms
alpha = np.log10(cat['Total_flux']/cat['Si1.4GHz']) / np.log10(nu1/nu2)
f,ax = pp.paper_single_ax()
ax.semilogx()
c = ax.scatter(cat['Total_flux'], alpha, marker='.',alpha=0.7)
#c = ax.scatter(cat['Total_flux'][sel], alpha[sel], marker='.')
mx,my,mys,_ = med_in_bins(np.log10(cat['Total_flux']), alpha)
mx = 10**mx
ax.errorbar(mx,my,mys, color='k')
xr = np.linspace(3e-2,20,10)
yr = np.log(xr/0.00015836431647618383 ) / np.log10(nu1/nu2)
ax.plot(xr,yr)
pp.set_attrib(ax, xlabel='$S_{53}$ (Jy)', ylabel='$\\alpha_{53}^{1400}$', ylim=[-2.2,1.8], xlim=[3e-3,20])
pp.format_log_axis10(ax, axis='x')
pp.fig_save_many(f, 'beerze_plots/deepwsrt_hbashift_flux_ratio')


f,ax = pp.paper_single_ax()
c = ax.scatter(np.log10(cat['Si1.4GHz']), alpha, marker='.',alpha=0.7)
#c = ax.scatter(cat['Total_flux'][sel], alpha[sel], marker='.')
mx,my,mys,mxs = med_in_bins(np.log10(cat['Si1.4GHz']), alpha)
#mx = 10**mx
#mxs = 10**mxs
ax.errorbar(mx,my,mys,mxs, color='k',ls='none')
#ax.semilogx()
xr = np.linspace(1e-4,10,10)
alphamin = np.log10(f1min/xr) / np.log10(nu1/nu2)
alphamin2 = np.log10(f2min/xr) / np.log10(nu1/nu2)
ax.fill_between(np.log10(xr),alphamin, 10*np.ones(len(xr)), color='gray', alpha=0.2)
ax.fill_between(np.log10(xr),alphamin2, 10*np.ones(len(xr)) , color='gray', alpha=0.2)
pp.set_attrib(ax, xlabel='$\log_{10} S_{1400}$ (Jy)', ylabel='$\\alpha_{53}^{1400}$', ylim=[-2.2,1.8], xlim=[-4.2,0.])
pp.format_log_axis10(ax, axis='x')
pp.fig_save_many(f, 'beerze_plots/deepwsrt_hbashift_flux_ratio2')
