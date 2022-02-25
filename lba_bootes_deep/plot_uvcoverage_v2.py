#! /usr/bin/env python

import sys
#import pydal as dal
import matplotlib.pyplot as pl
import pyrap.tables as pt
import utils.plot_util as pp
import numpy as np
import matplotlib

#check usage
if ( len(sys.argv) < 2 ) or ( len(sys.argv) > 4 ):
    print "Usage:"
    print "\tuv_coverage.py <file> [antenna1] [antenna2]"
    print "\t<> required"
    print "\t[] optional"
    print ""
    sys.exit(1)


MS1 = sys.argv[1]
#MS2 = sys.argv[2]



#msds= dal.dalDataset()
#if ( True != msds.open(sys.argv[1]) ):
        #print "ERROR: Could not open file: " + sys.argv[1]
        #print "       Please check the file and try again."
        #sys.exit(1)

#tablename = "MAIN";

#if ( 3 == len(sys.argv) ):  # if a single antenna value is given
    #msds.setFilter("UVW","ANTENNA1 = " + sys.argv[2])
    #title("UV coverage for" + sys.argv[1] + "\nAntenna(" + sys.argv[2] \
      #+ ")")

#elif ( 4 == len(sys.argv) ):  # if two antenna values are given
    #msds.setFilter("UVW","ANTENNA1 = " + sys.argv[2] + " AND ANTENNA2 = " + sys.argv[3])
    #title("UV coverage for" + sys.argv[1] + "\nAntenna pair(" + sys.argv[2] \
      #+ "," + sys.argv[3] + ")")
  
#else:  # if no antenna values are given
    #flagger="FALSE"
   ## msds.setFilter("UVW")
    #msds.setFilter("UVW")
    #title("UV coverage for\n" + sys.argv[1] )

#maintable = msds.openTable( tablename )

#figname = MS.replace('.ms','').replace('.MS','')+'.uvcoverage.png'
#fignamewave = MS.replace('.ms','').replace('.MS','')+'.uvcoverage.wave.png'
figname = 'uvcoverage'
fignamewave = 'uvcoverage.wave'
fignamewave2 = 'uvcoverage.wave1freq'
fignamewave3 = 'uvcoverage.wavedens'
maintable1 = pt.table(MS1)
chantable1 = pt.table(MS1+'/SPECTRAL_WINDOW/')
#maintable2 = pt.table(MS2)
#chantable2 = pt.table(MS2+'/SPECTRAL_WINDOW/')

# get the UVW column
fuvw_column1 = maintable1.getcol("UVW")
# get the freq column
ffreq_col1 = chantable1.getcol('CHAN_FREQ')
ffreq_col1 = ffreq_col1.flatten()
# get the UVW column
#uvw_column2 = maintable2.getcol("UVW")
## get the freq column
#freq_col2 = chantable2.getcol('CHAN_FREQ')
#freq_col2 = freq_col2.flatten()


color_list = pl.cm.BrBG(np.linspace(0, 1, 2))
color_list = pl.cm.cubehelix(np.linspace(0, 1, 4))
cols = color_list[1:]

print "{n} freqeuncies".format(n=len(ffreq_col1))
print "{n} times".format(n=len(fuvw_column1))
print

freq_samp = 400
time_samp = 10
freq_col1 = ffreq_col1[0:len(ffreq_col1):freq_samp]
uvw_column1 = fuvw_column1[0:len(fuvw_column1):time_samp,:]
#freq_col2 = freq_col2[0:len(freq_col2):freq_samp]
#uvw_column2 = uvw_column2[0:len(uvw_column2):time_samp,:]
vdata1=uvw_column1[:,1]
udata1=uvw_column1[:,0]



Npnts1 = len(uvw_column1)*len(freq_col1)
print "selection: {n}".format(n=figname)
print "{n} freqeuncies".format(n=len(freq_col1))
print "{n} times".format(n=len(uvw_column1))
print "{n} points to be plotted".format(n=Npnts1)

# get the data from the column
#data = uvw_column.data()
#vdata2=uvw_column2[:,1]
#udata2=uvw_column2[:,0]
#print vdata[5]
#print udata[5]
#print uvw_column[5]

#f=pl.figure()
#ax = pl.subplot(111)
f,ax = pp.paper_single_ax()
pl.minorticks_on()
pl.axis('equal')
# plot the data
ax.set_xlabel("U [km]")
ax.set_ylabel("V [km]")
ax.plot(udata1/1000.,vdata1/1000.,'.',c=cols[0], markersize=2, alpha=0.1)
ax.plot(-udata1/1000.,-vdata1/1000.,'.',c=cols[1], markersize=2, alpha=0.1)
#ax.plot(udata2/1000.,vdata2/1000.,'b,')
#ax.plot(-udata2/1000.,-vdata2/1000.,'r,')
#ax.plot(-udata,-vdata,'b,')
pp.fig_save_many(f, figname)


freq_samp = 40
time_samp = 10
freq_col1 = ffreq_col1[0:len(ffreq_col1):freq_samp]
uvw_column1 = fuvw_column1[0:len(fuvw_column1):time_samp,:]
vdata1=uvw_column1[:,1]
udata1=uvw_column1[:,0]
Npnts1 = len(uvw_column1)*len(freq_col1)
print "selection: {n}".format(n=fignamewave)
print "{n} freqeuncies".format(n=len(freq_col1))
print "{n} times".format(n=len(uvw_column1))
print "{n} points to be plotted".format(n=Npnts1)

#f=pl.figure()
#ax = pl.subplot(111)
f,ax = pp.paper_single_ax()
pl.minorticks_on()
pl.axis('equal')
# plot the data
ax.set_xlabel("$u$ [k$\lambda$]")
ax.set_ylabel("$v$ [k$\lambda$]")
for freq in freq_col1:
    wave = 1000*2.998e8/freq  # freq in Hz
    udatawave = udata1/wave
    vdatawave = vdata1/wave
    ax.plot(udatawave,vdatawave,'.',c=cols[0], markersize=2, alpha=0.1)
    ax.plot(-udatawave,-vdatawave,'.',c=cols[1], markersize=2, alpha=0.1)
#for freq in freq_col2:
    #wave = 1000*2.998e8/freq  # freq in Hz
    #udatawave = udata2/wave
    #vdatawave = vdata2/wave
    #ax.plot(udatawave,vdatawave,'b,', alpha=0.1)
    #ax.plot(-udatawave,-vdatawave,'r,', alpha=0.1)
#ax.plot(-udata,-vdata,'b,')
pp.fig_save_many(f, fignamewave)



freq_samp = 10
time_samp = 10
freq_col1 = ffreq_col1[0:len(ffreq_col1):freq_samp]
uvw_column1 = fuvw_column1[0:len(fuvw_column1):time_samp,:]
vdata1=uvw_column1[:,1]
udata1=uvw_column1[:,0]
Npnts1 = len(uvw_column1)*len(freq_col1)
print "selection: {n}".format(n=fignamewave2)
print "{n} freqeuncies".format(n=len(freq_col1))
print "{n} times".format(n=len(uvw_column1))
print "{n} points to be plotted".format(n=Npnts1)
f,ax = pp.paper_single_ax()
pl.minorticks_on()
pl.axis('equal')
# plot the data
ax.set_xlabel("$u$ [k$\lambda$]")
ax.set_ylabel("$v$ [k$\lambda$]")
for freq in freq_col1:
    wave = 1000*2.998e8/freq  # freq in Hz
    udatawave = udata1/wave
    vdatawave = vdata1/wave
    ax.plot(udatawave,vdatawave,'.',c=cols[0], markersize=2, alpha=0.01)
    ax.plot(-udatawave,-vdatawave,'.',c=cols[1], markersize=2, alpha=0.01)
#for freq in freq_col2:
    #wave = 1000*2.998e8/freq  # freq in Hz
    #udatawave = udata2/wave
    #vdatawave = vdata2/wave
    #ax.plot(udatawave,vdatawave,'b,', alpha=0.1)
    #ax.plot(-udatawave,-vdatawave,'r,', alpha=0.1)
#ax.plot(-udata,-vdata,'b,')
pp.fig_save_many(f, fignamewave2)

#freq_samp = 1
#time_samp = 1
#freq_col1 = ffreq_col1[0:len(ffreq_col1):freq_samp]
#uvw_column1 = fuvw_column1[0:len(fuvw_column1):time_samp,:]
#vdata1=uvw_column1[:,1]
#udata1=uvw_column1[:,0]
#Npnts1 = len(uvw_column1)*len(freq_col1)
#print "selection: {n}".format(n=fignamewave3)
#print "{n} freqeuncies".format(n=len(freq_col1))
#print "{n} times".format(n=len(uvw_column1))
#print "{n} points to be plotted".format(n=Npnts1)

##f=pl.figure()
##ax = pl.subplot(111)
#f,ax = pp.paper_single_ax()
#pl.minorticks_on()
#pl.axis('equal')
## plot the data
#ax.set_xlabel("$u$ [k$\lambda$]")
#ax.set_ylabel("$v$ [k$\lambda$]")
#uu = []
#vv = []
#for freq in freq_col1:
    #wave = 1000*2.998e8/freq  # freq in Hz
    #udatawave = udata1/wave
    #vdatawave = vdata1/wave
    ##ax.plot(udatawave,vdatawave,'b,', alpha=0.5)
    ##ax.plot(-udatawave,-vdatawave,'r,', alpha=0.5)
    #uu.append(udatawave)
    #uu.append(-udatawave)
    #vv.append(vdatawave)
    #vv.append(-vdatawave)
#uu = np.array(uu).flatten()
#vv = np.array(vv).flatten()

#c= ax.hexbin(uu,vv,bins='log',cmap=matplotlib.cm.cubehelix_r)
#cbar=colorbar(c)
#cbar.set_label('\log density')
##for freq in freq_col2:
    ##wave = 1000*2.998e8/freq  # freq in Hz
    ##udatawave = udata2/wave
    ##vdatawave = vdata2/wave
    ##ax.plot(udatawave,vdatawave,'b,', alpha=0.1)
    ##ax.plot(-udatawave,-vdatawave,'r,', alpha=0.1)
##ax.plot(-udata,-vdata,'b,')
#pp.fig_save_many(f, fignamewave3)


pl.show()
