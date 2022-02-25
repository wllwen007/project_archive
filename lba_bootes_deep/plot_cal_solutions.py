import utils.plotting as pp
from losoto.h5parm import h5parm
from losoto.lib_operations import *
import matplotlib.pyplot as plt

import numpy as np

import seaborn as sns
import matplotlib as mpl
pcmap = mpl.colors.ListedColormap(sns.color_palette("husl", 256))

pcmap = mpl.cm.hsv

egobs = 'L667916'

calfile = '/net/lofar2/data1/wwilliams/surveys/piperun/Pre-Facet-Calibrator-{obs}-v3-jul19/results/cal_values/cal_solutions.h5'.format(obs=egobs)

ioncalfile = '/net/bovenrijn/data1/wwilliams/surveys/pipework/Pre-Facet-Calibrator-L667910-sep19/instrument.h5imp_cal_ion'

instrument.h5imp_cal_FR
instrument.h5imp_cal_bandpass

H = h5parm(calfile, readonly=False)
Hion = h5parm(ioncalfile, readonly=False)


refAnt = 'CS001LBA'

pltAnts = ['CS002LBA', 'RS305LBA', 'RS307LBA', 'RS310LBA']
Nplt = len(pltAnts)

solsetion = Hion.getSolset('sol000')
#['amplitude000',
 #'clock',
 #'phase000',
 #'phaseOrig',
 #'phase_offset',
 #'tec',
 #'tec3rd']


solset = H.getSolset('calibrator')
#['phaseResid000',
 #'tec000',
 #'clock000',
 #'phase000',
 #'modelphase000',
 #'clockSmooth000']

bptab = solset.getSoltab('bandpass')
ftab = solset.getSoltab('faraday')
phtab = solset.getSoltab('phaseOrig')
patab = solset.getSoltab('polalign')

resphasesoltab = solsetion.getSoltab('phaseOrig')
tecsoltab = solsetion.getSoltab('tec')
clocksoltab = solsetion.getSoltab('clock')

### pass - not implemented

phase, phase_axVals = phtab.getValues(reference=refAnt)
phase = normalize_phase(phase)
pa, pa_axVals = patab.getValues(reference=refAnt)
fr, fr_axVals = ftab.getValues(reference=refAnt)
bp, bp_axVals = bptab.getValues(reference=refAnt)
phaseres, phaseres_axVals = resphasesoltab.getValues(reference=refAnt)
#phaseres = normalize_phase(phaseres)
#modphase, modphase_axVals = modphasesoltab.getValues(reference=refAnt)
#modphase = normalize_phase(modphase)
tec, tec_axVals = tecsoltab.getValues(reference=refAnt)
clock, clock_axVals = clocksoltab.getValues(reference=refAnt)

def subplotsd(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True,
             subplot_kw=None, gridspec_kw=None, dummy=None, shareall=None,
             xlab=None, ylab=None, **fig_kw):
    """
    Create a figure and a set of subplots

    This utility wrapper makes it convenient to create common layouts of
    subplots, including the enclosing figure object, in a single call.

    modified from amtplotlibs one
    """
    fig = plt.figure(**fig_kw)
    # dummy bg axis
    if dummy:
        axd = fig.subplots(1,1)
        # Turn off axis lines and ticks of the big subplot
        axd.spines['top'].set_color('none')
        axd.spines['bottom'].set_color('none')
        axd.spines['left'].set_color('none')
        axd.spines['right'].set_color('none')
        #axd.minorticks_off()
        ##axd.majorticks_off()
        #axd.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
        axd.tick_params(axis='both',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        color='none',
                        labelcolor='none')
                        #bottom=False,      # ticks along the bottom edge are off
                        #top=False,         # ticks along the top edge are off
                        #labelbottom=False) # labels along the bottom edge are off

    axs = fig.subplots(nrows=nrows, ncols=ncols, sharex=sharex, sharey=sharey,
                       squeeze=squeeze, subplot_kw=subplot_kw,
                       gridspec_kw=gridspec_kw)
    fig.subplots_adjust(left=0.05, right=0.9, wspace=0.1)
    
    if shareall:
        axd.get_shared_x_axes().join(axd, axs[0])
        
    if xlab:
        # Set common labels
        axd.set_xlabel(xlab)
    if ylab:
        axd.set_ylabel(ylab)
        
    
    return fig, axs, axd


plt.minorticks_on()
    
fig, axs, axd = subplotsd(1, Nplt, sharex=True, sharey=True, dummy=True, figsize=(14,4), xlab='time (hr)', ylab='Freq (MHz)')

for iAnt, pltAnt in enumerate(pltAnts):
    ai = np.where(phase_axVals['ant'] == pltAnt)[0][0]
    time = phase_axVals['time']
    freq = phase_axVals['freq']
    
    time = (time-time[0]) / (3600.)
    freq = freq / 1e6
    
    phasei = phase[:, :, ai, 0]
    
    axi = axs[iAnt]
    
    c = axi.imshow(phasei.T, origin='lower', extent=[time[0],time[-1], freq[0], freq[-1]], aspect='auto', cmap=pcmap, vmin=-3.14, vmax=3.14)
    
    axi.set_title(pltAnt) #, pad=-20)

rect = 0.925, 0.1, 0.025, 0.8
cax = fig.add_axes(rect)
cbar = plt.colorbar(c,cax=cax)
cbar.set_label('Phase (rad)')
pp.fig_save_many(fig, 'beerze_plots/eg_cal_sols_{obs:s}_phase'.format(obs=egobs))
    


    
fig, axs, axd = subplotsd(1, Nplt, sharex=True, sharey=True, dummy=True, figsize=(14,4), xlab='time (hr)', ylab='Freq (MHz)')

for iAnt, pltAnt in enumerate(pltAnts):
    ai = np.where(pa_axVals['ant'] == pltAnt)[0][0]
    time = pa_axVals['time']
    freq = pa_axVals['freq']
    
    time = (time-time[0]) / (3600.)
    freq = freq / 1e6
    
    pai = pa[:, :, ai, 1]
    
    axi = axs[iAnt]
    
    c = axi.imshow(pai.T, origin='lower', extent=[time[0],time[-1], freq[0], freq[-1]], aspect='auto', cmap=pcmap, vmin=-3.14, vmax=3.14)
    
    axi.set_title(pltAnt) #, pad=-20)

rect = 0.925, 0.1, 0.025, 0.8
cax = fig.add_axes(rect)
cbar = plt.colorbar(c,cax=cax)
cbar.set_label('Phase (rad)')
pp.fig_save_many(fig, 'beerze_plots/eg_cal_sols_{obs:s}_pa'.format(obs=egobs))
    





fig, axs, axd = subplotsd(1, Nplt, sharex=True, sharey=True, dummy=True, figsize=(14,4), xlab='time (hr)', ylab='Freq (MHz)')

for iAnt, pltAnt in enumerate(pltAnts):
    ai = np.where(phase_axVals['ant'] == pltAnt)[0][0]
    time = phase_axVals['time']
    freq = phase_axVals['freq']
    
    time = (time-time[0]) / (3600.)
    freq = freq / 1e6
    
    phasei = phaseres[:, :, ai, 0]
    
    axi = axs[iAnt]
    
    c = axi.imshow(phasei.T, origin='lower', extent=[time[0],time[-1], freq[0], freq[-1]], aspect='auto', cmap=pcmap, vmin=-3.14, vmax=3.14)
    
    axi.set_title(pltAnt) #, pad=-20)

rect = 0.925, 0.1, 0.025, 0.8
cax = fig.add_axes(rect)
cbar = plt.colorbar(c,cax=cax)
cbar.set_label('Residual Phase (rad)')
pp.fig_save_many(fig, 'beerze_plots/eg_cal_sols_{obs:s}_phaseres'.format(obs=egobs))    




fig, axs, axd = subplotsd(1, Nplt, sharex=True, sharey=True, dummy=True,  xlab='time (hr)', ylab='TEC (TECu)', figsize=(14,4))

for iAnt, pltAnt in enumerate(pltAnts):
    ai = np.where(tec_axVals['ant'] == pltAnt)[0][0]
    time = tec_axVals['time']
    
    time = (time-time[0]) / (3600.)
    
    teci = tec[:, ai]
    
    axi = axs[iAnt]
    
    axi.scatter(time, teci, s=1)
    
    axi.set_title(pltAnt) #, pad=-20)
    
pp.fig_save_many(fig, 'beerze_plots/eg_cal_sols_{obs:s}_tec'.format(obs=egobs))    
    
fig, axs, axd = subplotsd(1, Nplt, sharex=True, sharey=True, dummy=True, figsize=(14,4), xlab='time (hr)', ylab='Delay (ns)')


clock=clock*1e9

for iAnt, pltAnt in enumerate(pltAnts):
    ai = np.where(tec_axVals['ant'] == pltAnt)[0][0]
    time = clock_axVals['time']
    
    time = (time-time[0]) / (3600.)
    
    clocki = clock[:, ai]
    
    axi = axs[iAnt]
    
    axi.scatter(time, clocki, s=1)
    
    axi.set_title(pltAnt) #, pad=-20)
    
    
pp.fig_save_many(fig, 'beerze_plots/eg_cal_sols_{obs:s}_delay'.format(obs=egobs))    
    


# example
#f,ax = pp.paper_single_ax()
