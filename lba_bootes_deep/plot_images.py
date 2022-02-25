import matplotlib.pyplot as plt
import numpy as np
import aplpy as ap
import utils.plotting as pp
import utils.make_subim as mm
import astropy.io.fits as pf
from astropy.table import Table, Column



import configparser
config = configparser.ConfigParser()
config.read('/net/lofar2//data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/image.config')
dpath = config['SETUP']['dpath']
impath = '/net/lofar2//data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'
plotpath = config['SETUP']['plotpath']
plotpath='plots/'
imfits = impath+config['image']['imfits']

rmsfits = impath+config['image']['rmsfits']
regfile = impath+config['image']['regfile']
#regfile2 = config['image']['regfile2']


dat = pf.getdata(imfits)
head = pf.getheader(imfits)
dat = dat[np.isfinite(dat)]

ra = head['CRVAL1']
dec = head['CRVAL2']
ps = head['CDELT2']
t = mm.extract_subim(imfits, ra,dec, 2.85)
trms = mm.extract_subim(rmsfits, ra,dec, 2.85)

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
ax1.recenter(ra,dec,width=5.5,height=5.5)

Zx = 0.5
ax1aBBox = [0.59,0.11,0.3,0.3]

ax1.ticks.set_color('k') 
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image',types=['.png','.pdf'],dpi=300)

ax1a = ap.FITSFigure(t[0],figure=f, subplot=ax1aBBox)
ax1a.ticks.hide()
ax1a.tick_labels.hide()
ax1a.axis_labels.hide()
ax1a.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=vmin, vmid=vmid, vmax=vmax, stretch='arcsinh')
ax1a.recenter(ra,dec+0.25,width=Zx,height=Zx)
ax1.ax.indicate_inset_zoom(ax1a.ax)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_annotate',types=['.png','.pdf'],dpi=300)
#pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_annotate',types=['.png','.pdf'])
#f.



f = plt.figure(figsize=(13.28,12.))
ax1 = ap.FITSFigure(t[0],figure=f)
#ax1.show_colorscale(vmin=vmin, vmax=vmax, stretch='linear', cmap=plt.cm.cubehelix_r)
ax1.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=vmin, vmid=vmid, vmax=vmax, stretch='arcsinh')

ax1.add_colorbar()
ax1.colorbar.set_axis_label_text('Intensity (mJy/bm)') 
ax1.recenter(ra,dec,width=5.5,height=5.5)
ax1.ticks.set_color('k') 
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
ax1.recenter(ra,dec,width=1,height=1)
#f.tight_layout()
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_zoom1',types=['.png','.pdf'],dpi=300)

ax1.show_regions(facetregfile)
#ax1.show_regions(regfile2)
ax1.recenter(ra,dec,width=5.6,height=5.6)
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_with_facets')

ax1.show_regions(regfile)
#ax1.show_regions(regfile2)
ax1.recenter(ra,dec,width=5.6,height=5.6)
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_with_regions')

ax1.recenter(ra,dec,width=2,height=2)
#f.tight_layout()
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_zoom2')


ax1.recenter(ra+1,dec+1,width=2,height=2)
#f.tight_layout()
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_zoom2a')


ax1.recenter(ra-1,dec+1,width=2,height=2)
#f.tight_layout()
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_zoom2b')


ax1.recenter(ra-1,dec-1,width=2,height=2)
#f.tight_layout()
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_zoom2c')


ax1.recenter(ra+1,dec-1,width=2,height=2)
#f.tight_layout()
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_zoom2d')

ax1.recenter(ra,dec,width=1,height=1)
#f.tight_layout()
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.95)
pp.fig_save_many(f, plotpath+'bootes_deep_lba_image_zoom1')


