import matplotlib.pyplot as plt
import numpy as np
import aplpy as ap
import utils.plotting as pp
import utils.make_subim as mm
import astropy.io.fits as pf
from astropy.table import Table, Column

dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7/fin_im/'
imfits = 'image_full_ampphase_m_m.NS.int.restored.fits'
rmsfits = 'bootes_deep_lba.rms.fits'
regfile = 'facets_fin.reg'
regfile2 = '../selfcal_sources.ds9.reg'

dat = pf.getdata(imfits)
head = pf.getheader(imfits)
dat = dat[np.isfinite(dat)]

ra = head['CRVAL1']
dec = head['CRVAL2']
ps = head['CDELT2']
t = mm.extract_subim(imfits, ra,dec, 5.)
trms = mm.extract_subim(rmsfits, ra,dec, 5.)

def std_sigclip(x, nit=10, nclip=5.):
    
    x = x[np.isfinite(x)]
    
    for i in range(nit):
        std = np.std(x)
        mn = np.mean(x)
        
        x = x[ (np.abs (x-mn) - nclip*std < 0)]
    
    std = np.std(x)
    
    
    return std

t[0].data = t[0].data * 1000.  
dat = t[0].data
#rms = np.std(dat)
datdA = dat[3887:4887,3887:4887]
rms = std_sigclip(datdA) 
vmin = -1.*rms
vmax = 25.*rms
vmid = 10.*rms

dA = np.sum(np.isfinite(datdA))*ps*ps
print("rms {rms:.2f} is mJy/bm in central {dA:5.2f} sq deg".format(rms=rms, dA=dA))

A = np.sum(np.isfinite(dat))*ps*ps
print("Area is {A:5.2f} sq deg".format(A=A))

pp.paper_single()

f = plt.figure(figsize=(13.28,12.))
ax1 = ap.FITSFigure(t[0],figure=f)
#ax1.show_colorscale(vmin=vmin, vmax=vmax, stretch='linear', cmap=plt.cm.cubehelix_r)
ax1.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=vmin, vmid=vmid, vmax=vmax, stretch='arcsinh')

ax1.add_colorbar()
ax1.colorbar.set_axis_label_text('Intensity (mJy/bm)') 
#ax1.recenter(ra,dec,width=6.,height=6.)

ax1.add_beam()
ax1.beam.set_color('k')
ax1.beam.set_facecolor('w')
ax1.ticks.set_color('k') 
#f.tight_layout()

ax1.show_regions(regfile)
ax1.show_regions(regfile2)
#ax1.recenter(ra,dec,width=5.6,height=5.6)
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, 'beerze_plots/bootes_deep_lba_image_with_facets_selfcaldir')


