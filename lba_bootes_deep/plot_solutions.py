import utils.plotting as pp
from losoto.h5parm import h5parm
from losoto.lib_operations import *
import matplotlib.pyplot as plt

import numpy as np

import seaborn as sns
import matplotlib as mpl
pcmap = mpl.colors.ListedColormap(sns.color_palette("husl", 256))

pcmap = mpl.cm.hsv



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

pp.paper_single()
#'L667882'
#for egobs in [ 'L667924']:
for egobs in ['L667876', 'L667894', 'L667900', 'L667906', 'L667912', 'L667918', 'L667924']:
    
    for tc in range(8):
        t0 = 99e9
        print(tc)
        
        print(egobs)
        refAnt = 'CS001LBA'

        pltAnts = ['CS002LBA', 'RS305LBA', 'RS307LBA', 'RS310LBA']
        pltAnts = ['CS002LBA', 'RS106LBA', 'RS307LBA', 'RS508LBA']
        Nplt = len(pltAnts)
        plt.minorticks_on()
        figtec, axstec, axdtec = subplotsd(1, Nplt, sharex=True, sharey=True, dummy=True,  xlab='time (hr)', ylab='TEC (TECu)', figsize=(14,4))
        
            
        figdel, axsdel, axddel = subplotsd(1, Nplt, sharex=True, sharey=True, dummy=True, figsize=(14,4), xlab='time (hr)', ylab='Delay (ns)')  
                
                
        figph, axsph, axdph = subplotsd(1, Nplt, sharex=True, sharey=True, dummy=True, figsize=(14,4), xlab='time (hr)', ylab='Freq (MHz)')
        
        figrph, axsrph, axdrph = subplotsd(1, Nplt, sharex=True, sharey=True, dummy=True, figsize=(14,4), xlab='time (hr)', ylab='Freq (MHz)')
        figmph, axsmph, axdmph = subplotsd(1, Nplt, sharex=True, sharey=True, dummy=True, figsize=(14,4), xlab='time (hr)', ylab='Freq (MHz)')

    
    
        egt = f'TC{tc:02d}'
        egh5 = '/net/lofar2/data1/wwilliams/surveys/postcal/cal-ph-c0-{obs:s}_{t:s}.h5'.format(obs=egobs,t=egt)
        H = h5parm(egh5, readonly=False)


        solset = H.getSolset('sol000')
        #['phaseResid000',
        #'tec000',
        #'clock000',
        #'phase000',
        #'modelphase000',
        #'clockSmooth000']

        clocksoltab = solset.getSoltab('clock000')
        modphasesoltab = solset.getSoltab('modelphase000')
        phasesoltab = solset.getSoltab('phase000')
        resphasesoltab = solset.getSoltab('phaseResid000')
        tecsoltab = solset.getSoltab('tec000')

        freq = phasesoltab.getAxisValues('freq')
        phase, phase_axVals = phasesoltab.getValues(reference=refAnt)
        phase = normalize_phase(phase)
        phaseres, phaseres_axVals = resphasesoltab.getValues(reference=refAnt)
        phaseres = normalize_phase(phaseres)
        modphase, modphase_axVals = modphasesoltab.getValues(reference=refAnt)
        modphase = normalize_phase(modphase)
        tec, tec_axVals = tecsoltab.getValues(reference=refAnt)
        clock, clock_axVals = clocksoltab.getValues(reference=refAnt)



    

        for iAnt, pltAnt in enumerate(pltAnts):
            if pltAnt not in tec_axVals['ant']:
                print (pltAnt,'missing')
                continue
            ai = np.where(tec_axVals['ant'] == pltAnt)[0][0]
            time = tec_axVals['time']
            
            t0 = np.min((t0, time[0]))
            
            time = (time-t0) / (3600.)
            
            teci = tec[:, ai]
            
            axi = axstec[iAnt]
            
            axi.scatter(time, teci, s=1)
            
            axi.set_title(pltAnt) #, pad=-20)
            

        clock=clock*1e9

        for iAnt, pltAnt in enumerate(pltAnts):
            ai = np.where(tec_axVals['ant'] == pltAnt)[0][0]
            time = clock_axVals['time']
            
            time = (time-t0) / (3600.)
            
            clocki = clock[:, ai]
            
            axi = axsdel[iAnt]
            
            axi.scatter(time, clocki, s=1)
            
            axi.set_title(pltAnt) #, pad=-20)
            
            

        for iAnt, pltAnt in enumerate(pltAnts):
            ai = np.where(phase_axVals['ant'] == pltAnt)[0][0]
            time = phase_axVals['time']
            freq = phase_axVals['freq']
            
            time = (time-t0) / (3600.)
            freq = freq / 1e6
            
            phasei = phase[:, :, ai, 0]
            
            axi = axsph[iAnt]
            
            c = axi.imshow(phasei.T, origin='lower', extent=[time[0],time[-1], freq[0], freq[-1]], aspect='auto', cmap=pcmap, vmin=-3.14, vmax=3.14)
            
            axi.set_title(pltAnt) #, pad=-20)

        rect = 0.925, 0.1, 0.025, 0.8
        cax = figph.add_axes(rect)
        cbar = plt.colorbar(c,cax=cax)
        cbar.set_label('Phase (rad)')
            


        for iAnt, pltAnt in enumerate(pltAnts):
            ai = np.where(phase_axVals['ant'] == pltAnt)[0][0]
            time = phase_axVals['time']
            freq = phase_axVals['freq']
            
            time = (time-t0) / (3600.)
            freq = freq / 1e6
            
            phasei = phaseres[:, :, ai, 0]
            
            axi = axsrph[iAnt]
            
            c = axi.imshow(phasei.T, origin='lower', extent=[time[0],time[-1], freq[0], freq[-1]], aspect='auto', cmap=pcmap, vmin=-3.14, vmax=3.14)
            
            axi.set_title(pltAnt) #, pad=-20)

        rect = 0.925, 0.1, 0.025, 0.8
        cax = figrph.add_axes(rect)
        cbar = plt.colorbar(c,cax=cax)
        cbar.set_label('Residual Phase (rad)')   


        for iAnt, pltAnt in enumerate(pltAnts):
            ai = np.where(phase_axVals['ant'] == pltAnt)[0][0]
            time = phase_axVals['time']
            freq = phase_axVals['freq']
            
            time = (time-t0) / (3600.)
            freq = freq / 1e6
            
            phasei = modphase[:, :, ai, 0]
            
            axi = axsmph[iAnt]
            
            c = axi.imshow(phasei.T, origin='lower', extent=[time[0],time[-1], freq[0], freq[-1]], aspect='auto', cmap=pcmap, vmin=-3.14, vmax=3.14)
            
            axi.set_title(pltAnt) #, pad=-20)

        cax = figmph.add_axes(rect)
        cbar = plt.colorbar(c,cax=cax)
        cbar.set_label('Model phase (rad)')
            
   
        pp.fig_save_many(figtec, 'beerze_plots/eg_sols_{obs:s}_{t:s}_tec'.format(obs=egobs,t=egt))    
        pp.fig_save_many(figdel, 'beerze_plots/eg_sols_{obs:s}_{t:s}_delay'.format(obs=egobs,t=egt))  
        pp.fig_save_many(figph, 'beerze_plots/eg_sols_{obs:s}_{t:s}_phase'.format(obs=egobs,t=egt))
        pp.fig_save_many(figrph, 'beerze_plots/eg_sols_{obs:s}_{t:s}_phaseres'.format(obs=egobs,t=egt)) 
        pp.fig_save_many(figmph, 'beerze_plots/eg_sols_{obs:s}_{t:s}_phasemod'.format(obs=egobs,t=egt)) 

    #pp.fig_save_many(figtec, 'beerze_plots/eg_sols_{obs:s}_tec'.format(obs=egobs,t=egt))    
    #pp.fig_save_many(figdel, 'beerze_plots/eg_sols_{obs:s}_delay'.format(obs=egobs,t=egt))  
    #pp.fig_save_many(figph, 'beerze_plots/eg_sols_{obs:s}_phase'.format(obs=egobs,t=egt))
    #pp.fig_save_many(figrph, 'beerze_plots/eg_sols_{obs:s}_phaseres'.format(obs=egobs,t=egt)) 
    #pp.fig_save_many(figmph, 'beerze_plots/eg_sols_{obs:s}_phasemod'.format(obs=egobs,t=egt)) 
