import os
import utils.plotting as pp
import utils.cutouts as uc
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table, Column, MaskedColumn
import astropy.units as u
import astropy.io.fits as fits
import pyregion

#clobber = True
clobber = 0
stilts = '/net/beerze/data2/wwilliams/software/bin/stilts '

#lofar_cat_name = 'bootes_deep_lba.cat.fits'
lofar_cat_name = 'bootes_deep_lba_hbashift.cat.fits'

facet_regfile= 'facets_fin.reg'



## match to deep opt

m_lofar_cat_name = 'bootes_deep_lba.cat.match_{n}.fits'.format(n='deephba')
m_lofar_cat_name = 'bootes_deep_lba_hbashift.cat.match_{n}.fits'.format(n='deephba_to_lba')
deep_hba_cat_name = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.blanked.scaled.cat.fits'
mradius = 15.
if (not os.path.exists(m_lofar_cat_name)) or clobber:
    cmd = stilts+'tskymatch2 in1={in2} ra1=RA dec1=DEC ra2=RA dec2=DEC in2={in1} error={rad:f} out={out} find=best join=all1'.format(in1=lofar_cat_name, in2=deep_hba_cat_name, out=m_lofar_cat_name, rad=mradius)
    os.system(cmd)
    
m_lofar_cat_name_not = 'bootes_deep_lba_hbashift.cat.match_{n}.fits'.format(n='deephba_not_lba')
if (not os.path.exists(m_lofar_cat_name_not)) or clobber:
    cmd = stilts+'tskymatch2 in1={in2} ra1=RA dec1=DEC ra2=RA dec2=DEC in2={in1} error={rad:f} out={out} find=best join=1not2'.format(in1=lofar_cat_name, in2=deep_hba_cat_name, out=m_lofar_cat_name_not, rad=mradius)
    os.system(cmd)
    
deep_hba_cat_selfmatch_name = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.blanked.scaled.cat.selfmatch.fits'
mradius = 2*15.
if (not os.path.exists(deep_hba_cat_selfmatch_name)) or clobber:
    cmd = stilts+'tskymatch2 in1={in1} ra1=RA dec1=DEC ra2=RA dec2=DEC in2={in1} error={rad:f} out={out} find=all'.format(in1=deep_hba_cat_name, out=deep_hba_cat_selfmatch_name, rad=mradius)
    os.system(cmd)
    
# exclude the deep hba sources that are not isolated - will give spurious/wrong matches
hcat_dup = Table.read(deep_hba_cat_selfmatch_name)
hisol = hcat_dup['GroupSize'].mask ==True
hisol_sources = hcat_dup['Source_id_1'][hisol]

lcat = Table.read(lofar_cat_name)
hcat = Table.read(deep_hba_cat_name)
cat = Table.read(m_lofar_cat_name)
cat = cat[np.isfinite(cat['RA_2'])]
catnot = Table.read(m_lofar_cat_name_not)

# select only sources within the LBA coverage
incoverage = np.isfinite(uc.get_pixel_values(cat['RA_1'], cat['DEC_1'],'image_full_phase_m.NS_shift.app.facetRestored.blanked.fits'))
isol = np.array([s in hisol_sources for s in cat['Source_id_1']])
cat = cat[incoverage&isol]
incoverage = np.isfinite(uc.get_pixel_values(hcat['RA'], hcat['DEC'],'image_full_phase_m.NS_shift.app.facetRestored.blanked.fits'))
isol = np.array([s in hisol_sources for s in hcat['Source_id']])
hcat = hcat[incoverage&isol]
incoverage = np.isfinite(uc.get_pixel_values(catnot['RA'], catnot['DEC'],'image_full_phase_m.NS_shift.app.facetRestored.blanked.fits'))
isol = np.array([s in hisol_sources for s in catnot['Source_id']])
catnot = catnot[incoverage&isol]
    
print('LBA catalogue contains {0} sources'.format(len(lcat)))
print('Deep HBA catalogue contains {0} sources'.format(len(hcat)))
print('Matched Deep HBA to LBA catalogue contains {0} sources'.format(len(cat)))

f,ax = pp.paper_single_ax() 
ax.scatter(hcat['RA'], hcat['DEC'],marker='+',c='gray') 
ax.scatter(lcat['RA'], lcat['DEC'],marker='+',c='r') 
ax.scatter(cat['RA_2'], cat['DEC_2'],marker='x',c='b') 

hcat['Total_flux'] *= 1000  # to mJy
hcat['Maj'] *= 3600  # to arcsec

cat['Total_flux_1'] *= 1000  # to mJy
cat['Maj_1'] *= 3600  # to arcsec

catnot['Total_flux'] *= 1000  # to mJy
catnot['Maj'] *= 3600  # to arcsec

alpha = -0.6
Sfactor = (60./150.)**alpha

fbins = np.log10(np.logspace(-1, 2, 30))
abins = np.log10(np.logspace(0.5, 1.7, 29))
X, Y = np.meshgrid(fbins, abins)
Xl, Yl = np.meshgrid(10**fbins, 10**abins)
Hall, fedges, aedges = np.histogram2d(np.log10(hcat['Total_flux']*Sfactor), np.log10(hcat['Maj']),bins=(fbins,abins))
Hnot, fedges, aedges = np.histogram2d(np.log10(catnot['Total_flux']*Sfactor), np.log10(catnot['Maj']),bins=(fbins,abins))
Hmatched, fedges, aedges = np.histogram2d(np.log10(cat['Total_flux_1']*Sfactor), np.log10(cat['Maj_1']),bins=(fbins,abins))

f,ax = pp.paper_single_ax() 
#ax.set_aspect('equal')
c=ax.pcolormesh(X, Y, np.log10(Hall.T), cmap=plt.cm.Blues)
pp.set_attrib(ax,xlabel='$\log S_{60}$ (mJy)', ylabel='log Maj (arcsec)')
pp.add_colorbar(f,ax,c,label='$\log N$')
plt.subplots_adjust(right=0.85)


#f,ax = pp.paper_single_ax() 
##ax.set_aspect('equal')
#c=ax.pcolormesh(X, Y, (1.*Hnot/Hall).T, cmap=plt.cm.Blues)
#pp.set_attrib(ax,xlabel='$\log S_{60}$ (mJy)', ylabel='log Maj (arcsec)')
#pp.add_colorbar(f,ax,c,label='fraction missed')
#plt.subplots_adjust(right=0.85)


f,ax = pp.paper_single_ax() 
#ax.set_aspect('equal')
c=ax.pcolormesh(Xl, Yl, (1.*Hnot/Hall).T, cmap=plt.cm.Blues)
ax.set_xscale('log')
ax.set_yscale('log')
ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax.set_yticks([ 2, 5, 10,20])
ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
pp.set_attrib(ax,xlabel='$\log S_{60}$ (mJy)', ylabel='log Maj (arcsec)')
pp.add_colorbar(f,ax,c,label='fraction missed')
plt.subplots_adjust(right=0.85)



f,ax = pp.paper_single_ax() 
#ax.set_aspect('equal')
c=ax.pcolormesh(X, Y, (1.*Hmatched/Hall).T, cmap=plt.cm.Blues)
pp.set_attrib(ax,xlabel='$\log S_{60}$ (mJy)', ylabel='log Maj (arcsec)')
pp.add_colorbar(f,ax,c,label='fraction detected')
plt.subplots_adjust(right=0.85)



trms = Table.read('bootes_deep_lba.Area_rms.fits')
rms_area_interp = interp1d(trms['rms'],trms['Area_frac'],kind='cubic', bounds_error=False, fill_value=(0.,trms['Area_frac'][-1]))


sigfactor = 5.
smear_factor = 2

fbins = np.log10(np.logspace(-1, 4, 50))
fbinsl = 10**fbins
Nall_f,_ = np.histogram(np.log10(hcat['Total_flux']*Sfactor),bins=fbins)
Nmissed_f,_ = np.histogram(np.log10(catnot['Total_flux']*Sfactor),bins=fbins)

f,ax = pp.paper_single_ax() 
pp.plot_np_hist(ax, 10**fbins, 1-1.*Nmissed_f/Nall_f)

rms_compl = rms_area_interp(fbinsl/ sigfactor)   
rms_compl_smeared = rms_area_interp(fbinsl/ sigfactor / smear_factor)   
ax.plot(fbinsl, rms_compl)
ax.plot(fbinsl, rms_compl_smeared,marker='o')
ax.set_xscale('log')
ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
pp.set_attrib(ax,xlabel='$\log S_{60}$ (mJy)', ylabel='fraction detected')


#cbar=plt.colorbar(c)



cat.sort('Source_id_1')

sys.exit()
#trms['Area']
#trms['rms']

# select only the 'proper' detected sources.... what to do about completeness of higher order wavelet scales??
print('Catalogue contains {0} sources'.format(len(cat)))
sel = (cat['Peak_flux_1'] / cat['Isl_rms_1']) > 5.
cat = cat[sel]
print('Catalogue contains {0} sources after S/N selection'.format(len(cat)))
sel = cat['GroupSize'].mask
cat = cat[sel]
print('Catalogue contains {0} sources after group==1 selection'.format(len(cat)))

tflux = cat['Total_flux_1'] 
pflux = cat['Peak_flux_1']
rms = cat['Isl_rms_1']


#maskS = (cat['S_Code'] == 'S')  & (cat['Maj'] < 30.)
#cat = cat[maskS]
#print('Catalogue contains {0} sources after S, size selection'.format(len(cat)))

#maskS = ~(cat['ID_OPTICAL'].mask)
#cat = cat[maskS]
#print('Catalogue contains {0} sources after optical selection'.format(len(cat)))

sel = ((cat['Peak_flux_1'] / cat['Isl_rms_1']) > 10.) & (cat['Maj_1']*3600 < 30.)
sel = ((cat['Peak_flux_1'] / cat['Isl_rms_1']) > 7.5) & (cat['Maj_1']*3600 < 25.)& (cat['S_Code_1']=='S') 
print('Catalogue contains {0} sources with high snr and small size selection'.format(len(cat[sel])))

N = len(cat)

dRA = 3600.*(cat['RA_1'] - cat['RA_2']) *np.cos(cat['DEC_1']*np.pi/180.)
dDEC = 3600.*(cat['DEC_1'] - cat['DEC_2'])
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

f,ax = pp.paper_single_ax()
c = ax.scatter(hcat['RA'], hcat['DEC'], c='gray', marker='.')
c = ax.scatter(lcat['RA'], lcat['DEC'], c='k', marker='.')
c = ax.scatter(cat['RA_1'], cat['DEC_1'], c='C0', marker='.')
#cbar = plt.colorbar(c)
#cbar.set_label('dRA (arcsec)')
#ax.hlines(1.4, 0.1, 1e3, color='k')
pp.set_attrib(ax, xlabel='RA (J2000)', ylabel='DEC (J2000)')
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_astro_qual_offset_dRA_pos')

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
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_astro_qual_offset')

r2 = pyregion.open(facet_regfile)
patch_list, artist_list = r2.get_mpl_patches_texts(origin=0)

f,ax = pp.paper_single_ax()
for p in patch_list: ax.add_patch(p)
#c = ax.scatter(lcat['RA'], lcat['DEC'], c='gray', marker='.', s=5)
c = ax.scatter(cat['RA_1'][sel], cat['DEC_1'][sel], c=dRA[sel], marker='.', vmin=-2*dRAs, vmax=2*dRAs)
cbar = plt.colorbar(c, extend='both')
cbar.set_label('dRA (arcsec)')
#ax.hlines(1.4, 0.1, 1e3, color='k')
pp.invert_xlim(ax)
pp.set_attrib(ax, xlabel='RA (J2000)', ylabel='DEC (J2000)',xlim=[221.3,214.8], ylim=[31.6,37.06])
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_astro_qual_offset_dRA_pos')

patch_list, artist_list = r2.get_mpl_patches_texts(origin=0)
f,ax = pp.paper_single_ax()
for p in patch_list: ax.add_patch(p)
c = ax.scatter(cat['RA_1'][sel], cat['DEC_1'][sel], c=dDEC[sel], marker='.', vmin=-2*dDECs, vmax=2*dDECs)
cbar = plt.colorbar(c, extend='both')
cbar.set_label('dDEC (arcsec)')
#ax.hlines(1.4, 0.1, 1e3, color='k')
pp.invert_xlim(ax)
pp.set_attrib(ax, xlabel='RA (J2000)', ylabel='DEC (J2000)',xlim=[221.3,214.8], ylim=[31.6,37.06])
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_astro_qual_offset_dDEC_pos')

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

nu1 = 53.
nu2 = 146.
f1min = 0.000646 * 5  # min rms
f2min = 0.001116 * 5  # median rms
alpha = np.log10(cat['Total_flux_1']/cat['Total_flux_2']) / np.log10(nu1/nu2)
f,ax = pp.paper_single_ax()
ax.semilogx()
c = ax.scatter(cat['Total_flux_1'], alpha, marker='.',alpha=0.7)
#c = ax.scatter(cat['Total_flux_1'][sel], alpha[sel], marker='.')
mx,my,mys,_ = med_in_bins(np.log10(cat['Total_flux_1']), alpha)
mx = 10**mx
ax.errorbar(mx,my,mys, color='k')
xr = np.linspace(3e-2,20,10)
yr = np.log(xr/0.00015836431647618383 ) / np.log10(nu1/nu2)
ax.plot(xr,yr)
pp.set_attrib(ax, xlabel='$S_{53}$ (Jy)', ylabel='$\\alpha_{53}^{146}$', ylim=[-2.2,1.8], xlim=[3e-3,20])
pp.format_log_axis10(ax, axis='x')
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_flux_ratio')


f,ax = pp.paper_single_ax()
c = ax.scatter(np.log10(cat['Total_flux_2']), alpha, marker='.',alpha=0.7)
#c = ax.scatter(cat['Total_flux_1'][sel], alpha[sel], marker='.')
mx,my,mys,mxs = med_in_bins(np.log10(cat['Total_flux_1']), alpha)
#mx = 10**mx
#mxs = 10**mxs
ax.errorbar(mx,my,mys,mxs, color='k',ls=None)
#ax.semilogx()
xr = np.linspace(1e-3,10,10)
alphamin = np.log10(f1min/xr) / np.log10(nu1/nu2)
alphamin2 = np.log10(f2min/xr) / np.log10(nu1/nu2)
ax.fill_between(np.log10(xr),alphamin, 10*np.ones(len(xr)), color='gray', alpha=0.2)
ax.fill_between(np.log10(xr),alphamin2, 10*np.ones(len(xr)) , color='gray', alpha=0.2)
pp.set_attrib(ax, xlabel='$\log_{10} S_{146}$ (Jy)', ylabel='$\\alpha_{53}^{146}$', ylim=[-2.2,1.8], xlim=[-3.2,0.9])
pp.format_log_axis10(ax, axis='x')
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_flux_ratio2')
