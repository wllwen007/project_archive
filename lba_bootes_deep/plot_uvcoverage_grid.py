#! /usr/bin/env python

import sys
import matplotlib.pyplot as pl
import matplotlib as mpl
import casacore.tables as pt
from utils import plotting as pp
import numpy as np
import matplotlib

#check usage
if ( len(sys.argv) < 2 ) :
    print("Usage:")
    print("\tuv_coverage.py <file/s>")
    print("\t<> required")
    print("")
    sys.exit(1)


MSs = sys.argv[1:]




Nu = Nv = 2001
u0 = -25
v0 = -30
u1 = 25
v1 = 30
ugrid = np.linspace(u0,u1,Nu)
vgrid = np.linspace(v0,v1,Nv)




griddeduv = None
name = MSs[0].split('/')[-1].split('_')[0]

figname = 'beerze_plots/uvcoverage.grid.'+name

for MS1 in MSs:

    maintable1 = pt.table(MS1)
    chantable1 = pt.table(MS1+'/SPECTRAL_WINDOW/')

    # get the UVW column
    fuvw_column1 = maintable1.getcol("UVW")
    # get the freq column
    ffreq_col1 = chantable1.getcol('CHAN_FREQ')
    ffreq_col1 = ffreq_col1.flatten()


    color_list = pl.cm.BrBG(np.linspace(0, 1, 2))
    color_list = pl.cm.cubehelix(np.linspace(0, 1, 4))
    cols = color_list[1:]

    print("{n} frequencies".format(n=len(ffreq_col1)))
    print("{n} times".format(n=len(fuvw_column1)))
    print()

    freq_samp = 2
    time_samp = 2
    freq_col1 = ffreq_col1[0:len(ffreq_col1):freq_samp]
    freq_col1 = freq_col1[(freq_col1>34.e6)&(freq_col1<75.e6)]
    uvw_column1 = fuvw_column1[0:len(fuvw_column1):time_samp,:]
    vdata1=uvw_column1[:,1]
    udata1=uvw_column1[:,0]

    udatawave = (udata1*freq_col1[:,np.newaxis]/(1000*2.998e8)).flatten()
    vdatawave = (vdata1*freq_col1[:,np.newaxis]/(1000*2.998e8)).flatten()

    uvgrid, uvgrid1, uvgrid2 = np.histogram2d( udatawave,vdatawave ,bins=(ugrid, vgrid))
    
    if griddeduv is None:
        griddeduv = uvgrid.T
        griddeduv += uvgrid[::-1,::-1].T
    else:
        griddeduv += uvgrid.T
        griddeduv += uvgrid[::-1,::-1].T

Npnts1 = len(uvw_column1)*len(freq_col1)
print("selection: {n}".format(n=figname))
print("{n} freqeuncies".format(n=len(freq_col1)))
print("{n} times".format(n=len(uvw_column1)))
print("{n} points to be plotted".format(n=Npnts1))


#f=pl.figure()
#ax = pl.subplot(111)
f,ax = pp.paper_single_ax(AR=1)
pl.minorticks_on()
#pl.axis('equal')
# plot the data
ax.set_xlabel("$u$ (k$\lambda$)")
ax.set_ylabel("$v$ (k$\lambda$)")
c=ax.imshow(np.log10(griddeduv),cmap=pl.cm.cubehelix_r, extent=[u0,u1,v0,v1])
cbar = pl.colorbar(c)
cbar.set_label('$\log$ Density')
ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(2))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(2))
pp.fig_save_many(f, figname)

pl.show()
