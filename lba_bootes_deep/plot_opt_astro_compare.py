import os
import utils.plotting as pp
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table, Column, MaskedColumn
import astropy.units as u


def stats(x):
    N = len(x)
    m = np.mean(x)
    s = np.std(x)
    me = s/np.sqrt(N)
    return m,s,me


#clobber = True
clobber = 0
stilts = '/net/beerze/data2/wwilliams/software/bin/stilts '

lofar_cat_name = 'bootes_deep_lba.cat.fits'


lofar_dummy_cat_name = 'dummy.cat.fits'
m_lofar_dummy_cat_name = 'dummy.cat.match_{n}.fits'.format(n='pslocal')

dcat = Table.read(lofar_cat_name)
dcat['RA'] = dcat['RA'][np.random.shuffle(np.arange(len(dcat)))]
dcat.write(lofar_dummy_cat_name, overwrite=True)


## match to deep opt

m_lofar_cat_name = 'bootes_deep_lba.cat.match_{n}.fits'.format(n='pslocal')
cat_name = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/pslocal.fits'
mradius = 10.
if (not os.path.exists(m_lofar_cat_name)) or clobber:
    cmd = stilts+'tskymatch2 in1={in1} ra1=RA dec1=DEC ra2=ra dec2=dec in2={in2} error={rad:f} out={out} find=best1'.format(in1=lofar_cat_name, in2=cat_name, out=m_lofar_cat_name, rad=mradius)
    os.system(cmd)
if (not os.path.exists(m_lofar_dummy_cat_name)) or clobber:
    cmd = stilts+'tskymatch2 in1={in1} ra1=RA dec1=DEC ra2=ra dec2=dec in2={in2} error={rad:f} out={out} find=best1'.format(in1=lofar_dummy_cat_name, in2=cat_name, out=m_lofar_dummy_cat_name, rad=mradius)
    os.system(cmd)
    
    
lcat = Table.read(lofar_cat_name)



ocat = Table.read(cat_name)
    
cat = Table.read(m_lofar_cat_name)

rcat = Table.read(m_lofar_dummy_cat_name)
cat.sort('Source_id')


#trms['Area']
#trms['rms']

# select only the 'proper' detected sources.... what to do about completeness of higher order wavelet scales??
print('Catalogue contains {0} sources'.format(len(cat)))
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




dRA = 3600.*(rcat['RA_1'] - rcat['ra_2']) *np.cos(rcat['DEC_1']*np.pi/180.)
dDEC = 3600.*(rcat['DEC_1'] - rcat['dec_2'])
d = np.sqrt(dRA**2. + dDEC**2.)
dRAm, dRAs, dRAms = stats(dRA)
dDECm, dDECs, dDECms = stats(dDEC)
dm, ds, dms = stats(d)
print('random cat positions')
print('dRA = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dRAm, me=dRAms, s=dRAs))
print('dDEC = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dDECm, me=dDECms, s=dDECs))
print('d = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dm, me=dms, s=ds))




N = len(cat)
dRA = 3600.*(cat['RA_1'] - cat['ra_2']) *np.cos(cat['DEC_1']*np.pi/180.)
dDEC = 3600.*(cat['DEC_1'] - cat['dec_2'])
d = np.sqrt(dRA**2. + dDEC**2.)

dRAm, dRAs, dRAms = stats(dRA)
dDECm, dDECs, dDECms = stats(dDEC)
dm, ds, dms = stats(d)


print('dRA = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dRAm, me=dRAms, s=dRAs))
print('dDEC = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dDECm, me=dDECms, s=dDECs))
print('d = {m:.2f} \pm {me:.2f} (\sigma={s:.2f})'.format(m=dm, me=dms, s=ds))

dRAm, dRAs, dRAms = stats(dRA[sel])
dDECm, dDECs, dDECms = stats(dDEC[sel])
dm, ds, dms = stats(d[sel])

print ('after selection')
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
c = ax.scatter(ocat['ra'], ocat['dec'], c='gray', marker='.')
c = ax.scatter(lcat['RA'], lcat['DEC'], c='k', marker='.')
c = ax.scatter(cat['RA_1'], cat['DEC_1'], c='C0', marker='.')
cbar = plt.colorbar(c)
cbar.set_label('dRA (arcsec)')
#ax.hlines(1.4, 0.1, 1e3, color='k')
pp.set_attrib(ax, xlabel='RA (J2000)', ylabel='DEC (J2000)')
pp.fig_save_many(f, 'beerze_plots/pslocal_astro_qual_offset_dRA_pos')

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
pp.fig_save_many(f, 'beerze_plots/pslocal_astro_qual_offset')


f,ax = pp.paper_single_ax()
c = ax.scatter(cat['RA_1'][sel], cat['DEC_1'][sel], c=dRA[sel], marker='.', vmin=-2*dRAs, vmax=2*dRAs)
cbar = plt.colorbar(c, extend='both')
cbar.set_label('dRA (arcsec)')
#ax.hlines(1.4, 0.1, 1e3, color='k')
pp.invert_xlim(ax)
pp.set_attrib(ax, xlabel='RA (J2000)', ylabel='DEC (J2000)')
pp.fig_save_many(f, 'beerze_plots/pslocal_astro_qual_offset_dRA_pos')
f,ax = pp.paper_single_ax()
c = ax.scatter(cat['RA_1'][sel], cat['DEC_1'][sel], c=dDEC[sel], marker='.', vmin=-2*dDECs, vmax=2*dDECs)
cbar = plt.colorbar(c, extend='both')
cbar.set_label('dDEC (arcsec)')
#ax.hlines(1.4, 0.1, 1e3, color='k')

pp.invert_xlim(ax)
pp.set_attrib(ax, xlabel='RA (J2000)', ylabel='DEC (J2000)')
pp.fig_save_many(f, 'beerze_plots/pslocal_astro_qual_offset_dDEC_pos')



