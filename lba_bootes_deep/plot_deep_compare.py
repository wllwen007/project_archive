import os
import utils.plotting as pp
from utils.make_subim import extract_subim
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table, Column, MaskedColumn, join
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyregion
from astropy.wcs import WCS
from astropy.io import fits

#clobber = True
clobber = 1
stilts = '/net/beerze/data2/wwilliams/software/bin/stilts '

#lofar_cat_name = 'bootes_deep_lba.cat.fits'
lofar_cat_name = 'bootes_deep_lba_hbashift.cat.fits'

# TODO
#nudeep = 148.72 ??  - in header
nudeep = 144.

facet_regfile= 'facets_fin_finer.reg'

## match to deep opt

lofar_code_name = 'lba_match_deephba_codes.fits'

m_lofar_cat_name = 'bootes_deep_lba.cat.match_{n}.fits'.format(n='deephba')
m_lofar_cat_name = 'bootes_deep_lba_hbashift.cat.match_{n}.fits'.format(n='deephba')
cat_name = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.blanked.scaled.cat.fits'
fmatch_cat = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/Bootes/bootes_final_cross_match_catalogue-v1.0-agnclass.fits'
# this is these 2 matched in topcat
#fmatch_cat = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/Bootes/bootes_final_cross_match_catalogue-v1.0.fits'
#fagn_class_cat = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/SEDfits_v1.0/AGNclasses_Bootes_v1.fits'
m_lofar_deephba_cat_name = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/Bootes/radio_bootes_final_cross_match_catalogue-v1.0.fits'
'''
AGN_final=0  &  RadioAGN_final=0     -> star-forming galaxy
AGN_final=1  &  RadioAGN_final=0     -> 'radio-quiet' AGN
AGN_final=0  &  RadioAGN_final=1     -> 'jet-mode' radio AGN / LERG
AGN_final=1  &  RadioAGN_final=1     -> quasar-like radio AGN / HERG
AGN_final=-1 or RadioAGN_final=-1    -> no secure classification
'''


dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'
imfits = 'image_full_phase_m.NS_shift.int.facetRestored.blanked.fits'
hdu = fits.open(dpath+imfits)[0]
t = extract_subim(dpath+imfits,217.5,34.5,10)
wcs = WCS(t[0].header)

mradius = 10.
print ('match lba to deep hba')
if (not os.path.exists(m_lofar_cat_name)) or clobber:
    cmd = stilts+'tskymatch2 in1={in1} ra1=RA dec1=DEC ra2=RA dec2=DEC in2={in2} error={rad:f} out={out} find=all'.format(in1=lofar_cat_name, in2=cat_name, out=m_lofar_cat_name, rad=mradius)
    os.system(cmd)
    
mradius = 10.
print ('match deep hba to deep hba/opt')
if (not os.path.exists(m_lofar_deephba_cat_name)) or clobber:
    cmd = stilts+'tskymatch2 in1={in1} ra1=RA dec1=DEC ra2=RA dec2=DEC in2={in2} error={rad:f} out={out} find=all'.format(in1=cat_name, in2=fmatch_cat, out=m_lofar_deephba_cat_name, rad=mradius)
    os.system(cmd)
    
    
lcat = Table.read(lofar_cat_name)
lcat_codes = Table.read(lofar_code_name)
hcat = Table.read(cat_name)
hcatopt = Table.read(fmatch_cat)
#hcatopt = hcatopt[hcatopt['GroupSize'] == 1]

    
cat = Table.read(m_lofar_cat_name)
cat.sort('Source_id_1')

cat = join(cat , lcat_codes)


cat.add_column(Column(data=-1*np.ones(len(cat),dtype=int),name='has_hba_opt'))
cat.add_column(Column(data=-99*np.ones(len(cat),dtype=int),name='AGN_final'))
cat.add_column(Column(data=-99*np.ones(len(cat),dtype=int),name='RadioAGN_final'))
cat.add_column(Column(data=-1*np.ones(len(cat),dtype=float),name='Z_BEST'))
cat.add_column(Column(data=-1*np.ones(len(cat),dtype=float),name='SFR_conc'))
cat.add_column(Column(data=-1*np.ones(len(cat),dtype=float),name='Mass_conc'))
cat.add_column(Column(data=-1*np.ones(len(cat),dtype=float),name='Size'))
chcatopt = SkyCoord(hcatopt['RA'], hcatopt['DEC'], unit='deg')
for i in range(len(cat)):
    ci = SkyCoord(cat['RA_1'][i], cat['DEC_1'][i], unit='deg')
    sep = ci.separation(chcatopt)
    hi = np.argmin(sep)
    if (sep[hi] < (12.*u.arcsec)):
        if np.sum(sep < (12*u.arcsec)) == 1:
            cat['has_hba_opt'][i] = 1
            cat['AGN_final'][i] = hcatopt['AGN_final'][hi]
            cat['RadioAGN_final'][i] = hcatopt['RadioAGN_final'][hi]
            cat['Z_BEST'][i] = hcatopt['Z_BEST'][hi]
            cat['SFR_conc'][i] = hcatopt['SFR_conc'][hi]
            cat['Mass_conc'][i] = hcatopt['Mass_conc'][hi]
            if np.isfinite(hcatopt['LGZ_Size'][hi]):
                cat['Size'][i] = hcatopt['LGZ_Size'][hi]
            else:
                cat['Size'][i] = 2*hcatopt['DC_Maj'][hi]
            
        else:
            cat['has_hba_opt'][i] = 2
    #if cat['Source_id_2'][i] in hcatopt['Source_id']:   # source_id_2 is hba source_id
        #ii = np.where(cat['Source_id_2'][i]==hcatopt['Source_id'])[0]
        #if len(ii) == 1:
            #cat['has_hba_opt'] = 1
            #cat['AGN_final'][i] = hcatopt['AGN_final'][ii]
            #cat['RadioAGN_final'][i] = hcatopt['RadioAGN_final'][ii]
        #else:
            #cat['has_hba_opt'] = 2   # confused = multiple matches
            
           
cat['flag_SFG'] = (cat['AGN_final']==0)  &  (cat['RadioAGN_final']==0)     #-> star-forming galaxy
cat['flag_RQ'] = (cat['AGN_final']==1)  &  (cat['RadioAGN_final']==0 )    #-> 'radio-quiet' AGN
cat['flag_LERG'] = (cat['AGN_final']==0)  &  (cat['RadioAGN_final']==1)     #-> 'jet-mode' radio AGN / LERG
cat['flag_HERG'] = (cat['AGN_final']==1)  &  (cat['RadioAGN_final']==1 )    #-> quasar-like radio AGN / HERG
cat['flag_unc'] = (cat['AGN_final']==-1) | (cat['RadioAGN_final']==-1)    #-> no secure classification 
#sys.exit()
        

#trms['Area']
#trms['rms']

# select only the 'proper' detected sources.... what to do about completeness of higher order wavelet scales??
print('Catalogue contains {0} sources'.format(len(cat)))
sel = (cat['match_code']==1)
cat = cat[sel]
print('Catalogue contains {0} sources after clean match selection'.format(len(cat)))
sel = (cat['Peak_flux_1'] / cat['Isl_rms_1']) > 5.
cat = cat[sel]
print('Catalogue contains {0} sources after S/N selection'.format(len(cat)))
sel = cat['GroupSize'].mask
cat = cat[sel]
print('Catalogue contains {0} sources after group==1 selection'.format(len(cat)))

tflux = cat['Total_flux_1'] 
pflux = cat['Peak_flux_1']
rms = cat['Isl_rms_1']


print('Catalogue contains {0} SFGs '.format(np.sum(cat['flag_SFG']==1)))
print('Catalogue contains {0} HERGs '.format(np.sum(cat['flag_HERG']==1)))
print('Catalogue contains {0} LERGs '.format(np.sum(cat['flag_LERG']==1)))
print('Catalogue contains {0} RQ AGN '.format(np.sum(cat['flag_RQ']==1)))
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
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_astro_qual_offset_dRA_pos1')

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
ax.set_xlim(-8,8) 
ax.set_ylim(-7,7)  
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
pp.set_attrib(ax, xlabel='RA (deg)', ylabel='DEC (deg)',xlim=[221.3,214.8], ylim=[31.6,37.06])
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_astro_qual_offset_dRA_pos')

patch_list, artist_list = r2.get_mpl_patches_texts(origin=0)
f,ax = pp.paper_single_ax()
for p in patch_list: ax.add_patch(p)
c = ax.scatter(cat['RA_1'][sel], cat['DEC_1'][sel], c=dDEC[sel], marker='.', vmin=-2*dDECs, vmax=2*dDECs)
cbar = plt.colorbar(c, extend='both')
cbar.set_label('dDEC (arcsec)')
#ax.hlines(1.4, 0.1, 1e3, color='k')
pp.invert_xlim(ax)
pp.set_attrib(ax, xlabel='RA (deg)', ylabel='DEC (deg)',xlim=[221.3,214.8], ylim=[31.6,37.06])
#pp.set_attrib(ax,xlim=[221.3,214.8], ylim=[31.6,37.06])
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_astro_qual_offset_dDEC_pos')

def med_in_bins(x,y,bins=10,startx=10,endx=-10,xlim=-100):
    y = y[x>xlim]
    x = x[x>xlim]
    xt = x.copy()
    xt.sort()
    xt = xt[startx:endx]
    xbins = np.linspace(np.min(xt), np.max(xt),bins+1)
    mx = np.zeros(bins)
    my = np.nan*np.zeros(bins)
    mn = np.zeros(bins)
    mys = np.zeros(bins)
    mxs = np.zeros((2,bins))
    for i in range(bins):
        mx[i] = (xbins[i]+xbins[i+1]) /2.
        s = (x > xbins[i]) & (x <= xbins[i+1])
        mn[i] = np.sum(s)
        if mn[i] >= 2:
            my[i] = np.median(y[s])
            mys[i] = 1.2533 * np.std(y[s]) / np.sqrt(np.sum(s))
        #print(xbins[i],xbins[i+1])
        #print('xrange, mean = ',xbins[i],xbins[i+1],mx[i])
        #ax.plot(1000*(10**x[s]),y[s],marker='x')
        #for xx,yy in zip(1000*(10**x[s]),y[s]):
        #print(f'{xx:.2f} {yy:.2f}')
        #print('median y,x =',my[i],mx[i])
    for i in range(bins):
        mxs[0,i] = xbins[i+1]-mx[i]
        mxs[1,i] = mx[i]-xbins[i]
    
    return mx,my,mys,mxs,mn

#cat = cat[sel]

nu1 = 54.
nu2 = 144.   # used by pepe, rohit
f1min = 0.000646 * 5  # min rms
f2min = 0.001116 * 5  # median rms
f1mindeep = 3.8158618e-05 * 5  # min rms
f2mindeep = 7.0344686e-05 * 5  # median rms
alpha = np.log10(cat['Total_flux_1']/cat['Total_flux_2']) / np.log10(nu1/nu2)
f,ax = pp.paper_single_ax()
ax.semilogx()
c = ax.scatter(1000*cat['Total_flux_1'], alpha, marker='.',alpha=0.7)
#c = ax.scatter(cat['Total_flux_1'][sel], alpha[sel], marker='.')
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_1']), alpha)
#mx = 10**mx
#ax.errorbar(1000*mx,my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0]))]), color='k',linestyle='none')
X,Y,Yerr,Xerr = 1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))])
ax.errorbar(X,Y,Yerr,Xerr, color='k',ls='none')
np.savez('beerze_plots/deephba_hbashift_flux_ratio.npz',X,Y,Yerr,Xerr)
xr = np.linspace(3e-2,20,10)
yr = np.log(xr/0.00015836431647618383 ) / np.log10(nu1/nu2)
ax.plot(1000*xr,yr)
xr = np.logspace(-4,3,10)
alphamin = np.log10(f1mindeep/xr) / np.log10(nu2/nu1)
alphamin2 = np.log10(f2mindeep/xr) / np.log10(nu2/nu1)
ax.fill_between(1000*xr,alphamin, -10*np.ones(len(xr)), color='gray', alpha=0.2)
ax.fill_between(1000*xr,alphamin2, -10*np.ones(len(xr)) , color='gray', alpha=0.2)
pp.set_attrib(ax, xlabel='$S_{54}$ (mJy)', ylabel='$\\alpha_{54}^{144}$', ylim=[-2.2,1.8], xlim=[3,20e3])
pp.format_log_axis10(ax, axis='x')
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_flux_ratio')

selbright = sel & (cat['Total_flux_1']>20e-3)
frat = cat['Total_flux_1'] /(cat['Total_flux_2'] * (nu1/nu2)**-0.7)
patch_list, artist_list = r2.get_mpl_patches_texts(origin=0)


f,ax = pp.paper_single_ax()
for p in patch_list: ax.add_patch(p)
#c = ax.scatter(lcat['RA'], lcat['DEC'], c='gray', marker='.', s=5)
c = ax.scatter(cat['RA_1'][selbright], cat['DEC_1'][selbright], c=frat[selbright], marker='.', vmin=0.66, vmax=1.5)
cbar = plt.colorbar(c, extend='both')
cbar.set_label('Flux density ratio')
#ax.hlines(1.4, 0.1, 1e3, color='k')
pp.invert_xlim(ax)
pp.set_attrib(ax, xlabel='RA (deg)', ylabel='DEC (deg)',xlim=[221.3,214.8], ylim=[31.6,37.06])
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_flux_ratio_pos')


f,ax = pp.paper_single_ax()
c = ax.scatter(np.log10(cat['Total_flux_2']), alpha, marker='.',alpha=0.7)
#c = ax.scatter(cat['Total_flux_1'][sel], alpha[sel], marker='.')
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_2']), alpha)
#mx = 10**mx
#mxs = 10**mxs
ax.errorbar(mx,my,mys,mxs, color='k',ls='none')
#ax.semilogx()
xr = np.linspace(1e-3,10,10)
alphamin = np.log10(f1min/xr) / np.log10(nu1/nu2)
alphamin2 = np.log10(f2min/xr) / np.log10(nu1/nu2)
ax.fill_between(np.log10(xr),alphamin, 10*np.ones(len(xr)), color='gray', alpha=0.2)
ax.fill_between(np.log10(xr),alphamin2, 10*np.ones(len(xr)) , color='gray', alpha=0.2)
pp.set_attrib(ax, xlabel='$\log_{10} S_{144}$ (Jy)', ylabel='$\\alpha_{54}^{144}$', ylim=[-2.2,1.8], xlim=[-3.2,0.9])
pp.format_log_axis10(ax, axis='x')

ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_flux_ratio1')

f,ax = pp.paper_single_ax()
c = ax.scatter(1000*cat['Total_flux_2'], alpha, marker='.',alpha=0.7)
ax.semilogx()
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_1']), alpha)
#mx = 10**mx
#mxs = 10**mxs
ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='k',ls='none')
#ax.semilogx()
xr = np.linspace(4e-5,10,10)
alphamin = np.log10(f1min/xr) / np.log10(nu1/nu2)
alphamin2 = np.log10(f2min/xr) / np.log10(nu1/nu2)
ax.fill_between(1000*xr,alphamin, 10*np.ones(len(xr)), color='gray', alpha=0.2)
ax.fill_between(1000*xr,alphamin2, 10*np.ones(len(xr)) , color='gray', alpha=0.2)
ax.hlines(-0.5,0.5,9000,'k',linestyle='dotted')
pp.set_attrib(ax, xlabel='$S_{144\mathrm{\,MHz}}$ (mJy)', ylabel='$\\alpha_{54}^{144}$', ylim=[-2.2,1.8], xlim=[0.5,9000])
pp.format_log_axis10(ax, axis='x')

ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_flux_ratio2')




f,ax = pp.paper_single_ax()
m1 = cat['flag_SFG']
m2 = (cat['flag_LERG']) | (cat['flag_HERG'])
c = ax.scatter(1000*cat['Total_flux_2'], alpha, marker='.',c='gray' ,alpha=0.1)
#c = ax.scatter(1000*cat['Total_flux_2'][m3], alpha[m3], marker='.',c='C2' ,alpha=0.7,label='HERG')
c = ax.scatter(1000*cat['Total_flux_2'][m2], alpha[m2], marker='.',c='C1' ,alpha=0.7,label='AGN')
c = ax.scatter(1000*cat['Total_flux_2'][m1], alpha[m1], marker='.',c='C0' ,alpha=1,label='SFG')
ax.semilogx()
#mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_2'][m1]), alpha[m1],endx=-1,bins=3,xlim=-2.2)
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_2'][m1]), alpha[m1],startx=1,endx=-1,bins=3,xlim=-2.2)
ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C0',marker='s',mec='k',ecolor='k',ls='none')
#mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_2'][m3]), alpha[m3],bins=3)
#ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C2',marker='s',mec='k')#ls='none')
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_2'][m2]), alpha[m2],startx=1,bins=7,endx=-5,xlim=-2.2)
ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C1',marker='s',mec='k',ecolor='k',ls='none')
#ax.semilogx()
xr = np.linspace(4e-5,10,10)
alphamin = np.log10(f1min/xr) / np.log10(nu1/nu2)
alphamin2 = np.log10(f2min/xr) / np.log10(nu1/nu2)
ax.fill_between(1000*xr,alphamin, 10*np.ones(len(xr)), color='gray', alpha=0.2)
ax.fill_between(1000*xr,alphamin2, 10*np.ones(len(xr)) , color='gray', alpha=0.2)
ax.hlines(-0.5,0.5,9000,'k',linestyle='dotted')
pp.set_attrib(ax, xlabel='$S_{144\mathrm{\,MHz}}$ (mJy)', ylabel='$\\alpha_{54}^{144}$', ylim=[-2.,1.5], xlim=[0.9,3000],ytick_spacing=0.5,ytick_min_spacing=0.1)
pp.format_log_axis10(ax, axis='x')
ax.legend()
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_flux_ratio2_type')





f,ax = pp.paper_single_ax()
m1 = cat['flag_SFG']
m2 = ((cat['flag_LERG']) | (cat['flag_HERG'])) & (cat['Size']*3600. >= 2)
m3 = ((cat['flag_LERG']) | (cat['flag_HERG'])) & (cat['Size']*3600. < 2)
c = ax.scatter(1000*cat['Total_flux_2'], alpha, marker='.',c='gray' ,alpha=0.1)
#c = ax.scatter(1000*cat['Total_flux_2'][m3], alpha[m3], marker='.',c='C2' ,alpha=0.7,label='HERG')
c = ax.scatter(1000*cat['Total_flux_2'][m3], alpha[m3], marker='.',c='C2' ,alpha=0.7,label="AGN, Size $<2''$")
c = ax.scatter(1000*cat['Total_flux_2'][m2], alpha[m2], marker='.',c='C1' ,alpha=0.7,label="AGN, Size $\geq2''$")
c = ax.scatter(1000*cat['Total_flux_2'][m1], alpha[m1], marker='.',c='C0' ,alpha=1,label='SFG')
ax.semilogx()
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_2'][m1]), alpha[m1],startx=1,endx=-1,bins=3,xlim=-2.2)
ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C0',marker='s',mec='k',ecolor='k',ls='none')
#mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_2'][m3]), alpha[m3],bins=3)
#ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C2',marker='s',mec='k')#ls='none')
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_2'][m2]), alpha[m2],endx=-5,bins=5,xlim=-2.2)
ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C1',marker='s',mec='k',ecolor='k',ls='none')
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_2'][m3]), alpha[m3],endx=-5,bins=5,xlim=-2.2)
ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C2',marker='s',mec='k',ecolor='k',ls='none')
#ax.semilogx()
xr = np.linspace(4e-5,10,10)
alphamin = np.log10(f1min/xr) / np.log10(nu1/nu2)
alphamin2 = np.log10(f2min/xr) / np.log10(nu1/nu2)
ax.fill_between(1000*xr,alphamin, 10*np.ones(len(xr)), color='gray', alpha=0.2)
ax.fill_between(1000*xr,alphamin2, 10*np.ones(len(xr)) , color='gray', alpha=0.2)
ax.hlines(-0.5,0.5,9000,'k',linestyle='dotted')
pp.set_attrib(ax, xlabel='$S_{144\mathrm{\,MHz}}$ (mJy)', ylabel='$\\alpha_{54}^{144}$', ylim=[-2.,1.5], xlim=[0.9,3000],ytick_spacing=0.5,ytick_min_spacing=0.1)
pp.format_log_axis10(ax, axis='x')
ax.legend()
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_flux_ratio2_type2')



f,ax = pp.paper_single_ax()
m1 = cat['flag_SFG']
m2 = ((cat['flag_LERG']) | (cat['flag_HERG'])) & (cat['Size']*3600. >= 2)
m3 = ((cat['flag_LERG']) | (cat['flag_HERG'])) & (cat['Size']*3600. < 2)
c = ax.scatter(1000*cat['Total_flux_1'], alpha, marker='.',c='gray' ,alpha=0.1)
#c = ax.scatter(1000*cat['Total_flux_2'][m3], alpha[m3], marker='.',c='C2' ,alpha=0.7,label='HERG')
c = ax.scatter(1000*cat['Total_flux_1'][m3], alpha[m3], marker='.',c='C2' ,alpha=0.7,label="AGN, Size $<2''$")
c = ax.scatter(1000*cat['Total_flux_1'][m2], alpha[m2], marker='.',c='C1' ,alpha=0.7,label="AGN, Size $\geq2''$")
c = ax.scatter(1000*cat['Total_flux_1'][m1], alpha[m1], marker='.',c='C0' ,alpha=1,label='SFG')
ax.semilogx()
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_1'][m1]), alpha[m1],startx=1,endx=-1,bins=3,xlim=-2.2)
ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C0',marker='s',mec='k',ecolor='k',ls='none')
#mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_2'][m3]), alpha[m3],bins=3)
#ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C2',marker='s',mec='k')#ls='none')
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_1'][m2]), alpha[m2],startx=1,endx=-5,bins=5,xlim=-2.2)
ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C1',marker='s',mec='k',ecolor='k',ls='none')
mx,my,mys,mxs,mn = med_in_bins(np.log10(cat['Total_flux_1'][m3]), alpha[m3],startx=1,endx=-5,bins=5,xlim=-2.2)
ax.errorbar(1000*(10**mx),my,mys,np.array([1000*(10**(mx))-1000*(10**(mx-mxs[0])),1000*(10**(mx+mxs[1]))-1000*(10**(mx))]), color='C2',marker='s',mec='k',ecolor='k',ls='none')
#ax.semilogx()
xr = np.logspace(-4,3,10)
alphamin = np.log10(f1mindeep/xr) / np.log10(nu2/nu1)
alphamin2 = np.log10(f2mindeep/xr) / np.log10(nu2/nu1)
ax.fill_between(1000*xr,alphamin, -10*np.ones(len(xr)), color='gray', alpha=0.2)
ax.fill_between(1000*xr,alphamin2, -10*np.ones(len(xr)) , color='gray', alpha=0.2)
ax.hlines(-0.5,0.5,9000,'k',linestyle='dotted')
pp.set_attrib(ax, xlabel='$S_{54\mathrm{\,MHz}}$ (mJy)', ylabel='$\\alpha_{54}^{144}$', ylim=[-2.,1.5], xlim=[0.9,3000],ytick_spacing=0.5,ytick_min_spacing=0.1)
pp.format_log_axis10(ax, axis='x')
ax.legend()
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pp.fig_save_many(f, 'beerze_plots/deephba_hbashift_flux_ratio_type2')

# no obv diff in z, SFR
'''
f,ax = pp.paper_single_ax()
c=ax.scatter(cat['Mass_conc'],cat['SFR_conc'],c=cat['Z_BEST'],s=20)
plt.colorbar(c)

f,ax = pp.paper_single_ax()
c=ax.scatter(cat['Total_flux_2'],cat['Maj_2']*3600.,c=alpha,s=30,vmin=-2,vmax=1.5)
ax.semilogy()
plt.colorbar(c,extend='both')
'''


'''
Catalogue contains 1364 sources after group==1 selection
Catalogue contains 34 SFGs 
Catalogue contains 74 HERGs 
Catalogue contains 404 LERGs 
Catalogue contains 18 RQ AGN 
'''


