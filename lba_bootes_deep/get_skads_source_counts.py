from utils import plotting as pp
import numpy as np
import pylab as pl
import matplotlib as mpl
#from utils.fits import *
from astropy.table import Table

def  get_skads_counts(ttype=[[0,1,2],[0,1,2,3,4]], postprocess=False, nu=151):
  print('SKADS')
  sft = ['None', 'Quiescent', 'Starburst']
  agnt = ['None', 'RQ',' FRI','FRII','GPS']
  
  sftype = ttype[0]
  agntype = ttype[1]
  
  selsft = [ sft[i] for i in sftype ]
  selagnt = [ agnt[i] for i in agntype ]
  print(' agn : %s : %s ' %(agntype, ','.join(selagnt)))  #None [0] Radio-Quiet [1] FRI [2] FRII [3] GPS [4]
  print(' sf : %s : %s' %(sftype, ','.join(selsft)))   #None [0] Quiescent [1] Starburst [2]
  
  data = Table.read('/net/beerze//data2/wwilliams/projects/skads/wilman_cat_all.fits')
  
  z = data['redshift']
  sftypes =data['sftype']
  agntypes = data['agntype']
  ff = data['itot_151']
  ff2 = data['itot_610']
  
  ff = 10**ff
  ff2 = 10**ff2
  
  if nu != 151.:
    alp = np.log(ff/ff2)/np.log(151./610) 
    ff = ff*(nu/151)**alp
  
  
  
  m = pl.zeros(len(agntypes),dtype=int)
  for agnt in agntype:
      m += (agntypes == agnt)
  for sft in sftype:
      m += (sftypes == sft)
  mask = m>0
  
  #mask = pl.zeros(len(agntypes),dtype=bool)
  #for i in range(len(agntypes)):
    #if agntypes[i] in agntype:
      #mask[i] = 1
    #if sftypes[i] in sftype:
      #mask[i] = 1
  fluxes = ff
  flux = fluxes[mask]
  #flux = pl.ma.masked_where(mask, fluxes).compressed()
  print(' n sources = %i' %(len(flux)))
  #b1 = 1.e-4
  #b2 = 4.
  #bins = pl.logspace(pl.log10(b1),pl.log10(b2),35.)
  #bin_centres = pl.zeros( len(bins) - 1 )
  #for i in range(len(bins) - 1):
    #bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  b1 = 1.e-4
  b2 = 4
  ##b2 = 500
  ##b2 = 100
  #b1 = 1e-3
  #b2 = 100

  A = 51.18
  beta = -1.436
  Nperbin = 100.
  bins=[b1]
  while bins[-1] < b2: 
    bins.append(bins[-1]+(Nperbin/A)*bins[-1]**(-beta))
    
    
  #bins = np.linspace(b1,b2,100)
    
  bin_centres = pl.zeros( len(bins) - 1 )
  for i in range(len(bins) - 1):
    bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  nbins =  len(bin_centres) 
  bincount = pl.zeros(nbins)
  #dnbins = pl.zeros( nbins )
  counts = pl.zeros( nbins )
  count_errors = pl.zeros( nbins )
  #print bins
  #print nbins
  #print len(bins)
  #print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %('BinL', 'BinC', 'BinW', 'Count', 'Count/3', 'Norm Count', 'Area deg^2', 'Area st') 
  Area = 20.*20.
  degtosterad = (pl.pi/180.)**2.
  Area_st = Area * degtosterad
  for i in range( nbins ): 
      
    if postprocess:
        
        mask = ((flux>bins[i]) * (flux<bins[i+1]))
        fluxM = pl.ma.masked_where(1-mask,flux).compressed()
        zM = pl.ma.masked_where(1-mask,z).compressed()
        
        agntypeM = pl.ma.masked_where(1-mask,agntypes).compressed()
        sftypeM = pl.ma.masked_where(1-mask,sftypes).compressed()
                
        weights = pl.ones(len(fluxM))
        
        for j in range(len(fluxM)):
            if agntypeM[j] == 1:
                if zM[j] >1.9:
                    weights[j] = (1+z[j])**-3.27
            elif agntypeM[j] == 2:
                if zM[j] > 2.5:
                    weights[j] = (1+z[j])**-2.5
            elif sftypeM[j] == 3:
                if zM[j] > 4.8:
                    weights[j] = (1+z[j])**-7.29
            cnt = pl.sum(weights )
            
    else:
    
        cnt = pl.sum( (flux>bins[i]) * (flux<bins[i+1]) )
    # area of sim
    bincount[ i ] = cnt
    nbin_err = pl.sqrt( cnt )
    dnbin = bins[i+1] - bins[i]
    counts[ i ] = ((bincount[i]/dnbin)/Area_st)*(bin_centres[i]**2.5)
    count_errors[ i ] = ((nbin_err/dnbin)/Area_st)*(bin_centres[i]**2.5)
    #print '%.4f\t%.4f\t%.4f\t%i\t%i\t%.3f\t%.2f\t%.5f' %(bins[i], bin_centres[ i ], dnbin, cnt, cnt/3., counts[i], Area, Area_st) 
  
  return bin_centres ,counts, count_errors


def  get_skads_counts_sfonly(majorlim=0, less=True):
  print('SKADS SF')
    
  data = Table.read('/net/beerze//data2/wwilliams/projects/skads/wilman_cat_sf_comp.fits')
  
  z = data['redshift']
  ff = data['itot_151']
  
  major = data['major_axis']
  
  if less:
    mask = np.where(major<majorlim)
  else:
    mask = np.where(major>=majorlim)
      
  fluxes = 10.**ff
  flux = fluxes[mask]
  
  print(' n sources = %i' %(len(flux)))
  #b1 = 1.e-4
  #b2 = 4.
  #bins = pl.logspace(pl.log10(b1),pl.log10(b2),35.)
  #bin_centres = pl.zeros( len(bins) - 1 )
  #for i in range(len(bins) - 1):
    #bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  b1 = 1.e-4
  b2 = 4

  A = 51.18
  beta = -1.436
  Nperbin = 200.
  bins=[b1]
  while bins[-1] < b2: 
    bins.append(bins[-1]+(Nperbin/A)*bins[-1]**(-beta))
  bin_centres = pl.zeros( len(bins) - 1 )
  for i in range(len(bins) - 1):
    bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  nbins =  len(bin_centres) 
  bincount = pl.zeros(nbins)
  #dnbins = pl.zeros( nbins )
  counts = pl.zeros( nbins )
  count_errors = pl.zeros( nbins )
  #print bins
  #print nbins
  #print len(bins)
  #print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %('BinL', 'BinC', 'BinW', 'Count', 'Count/3', 'Norm Count', 'Area deg^2', 'Area st') 
  Area = 20.*20.
  degtosterad = (pl.pi/180.)**2.
  Area_st = Area * degtosterad
  for i in range( nbins ): 
      
    cnt = pl.sum( (flux>bins[i]) * (flux<bins[i+1]) )
    # area of sim
    bincount[ i ] = cnt
    nbin_err = pl.sqrt( cnt )
    dnbin = bins[i+1] - bins[i]
    counts[ i ] = ((bincount[i]/dnbin)/Area_st)*(bin_centres[i]**2.5)
    count_errors[ i ] = ((nbin_err/dnbin)/Area_st)*(bin_centres[i]**2.5)
    #print '%.4f\t%.4f\t%.4f\t%i\t%i\t%.3f\t%.2f\t%.5f' %(bins[i], bin_centres[ i ], dnbin, cnt, cnt/3., counts[i], Area, Area_st) 
  
  return bin_centres ,counts, count_errors

def  get_skads_counts_agnonly():
  print('SKADS AGN')
    
  data = Table.read('/net/beerze//data2/wwilliams/projects/skads/wilman_cat_agn.fits')
  
  z = data['redshift']
  ff = data['itot_151']
  
  
  #mask = np.where(major>majorlim)
  flux = 10.**ff
  #flux = fluxes[mask]
  
  print(' n sources = %i' %(len(flux)))
  #b1 = 1.e-4
  #b2 = 4.
  #bins = pl.logspace(pl.log10(b1),pl.log10(b2),35.)
  #bin_centres = pl.zeros( len(bins) - 1 )
  #for i in range(len(bins) - 1):
    #bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  b1 = 1.e-4
  b2 = 4

  A = 51.18
  beta = -1.436
  Nperbin = 100.
  bins=[b1]
  while bins[-1] < b2: 
    bins.append(bins[-1]+(Nperbin/A)*bins[-1]**(-beta))
  bin_centres = pl.zeros( len(bins) - 1 )
  for i in range(len(bins) - 1):
    bin_centres[ i ] = (bins[i] + bins[i+1])/2.
    
  nbins =  len(bin_centres) 
  bincount = pl.zeros(nbins)
  #dnbins = pl.zeros( nbins )
  counts = pl.zeros( nbins )
  count_errors = pl.zeros( nbins )
  #print bins
  #print nbins
  #print len(bins)
  #print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %('BinL', 'BinC', 'BinW', 'Count', 'Count/3', 'Norm Count', 'Area deg^2', 'Area st') 
  Area = 20.*20.
  degtosterad = (pl.pi/180.)**2.
  Area_st = Area * degtosterad
  for i in range( nbins ): 
    cnt = pl.sum( (flux>bins[i]) * (flux<bins[i+1]) )
    # area of sim
    bincount[ i ] = cnt
    nbin_err = pl.sqrt( cnt )
    dnbin = bins[i+1] - bins[i]
    counts[ i ] = ((bincount[i]/dnbin)/Area_st)*(bin_centres[i]**2.5)
    count_errors[ i ] = ((nbin_err/dnbin)/Area_st)*(bin_centres[i]**2.5)
    #print '%.4f\t%.4f\t%.4f\t%i\t%i\t%.3f\t%.2f\t%.5f' %(bins[i], bin_centres[ i ], dnbin, cnt, cnt/3., counts[i], Area, Area_st) 
  
  return bin_centres ,counts, count_errors

if __name__ == "__main__":

    name='SKADS'
    nu = 60.
    unit='mJy'

    scaleF = 1000.

    fig, ax1 = pp.paper_double_ax()
    pl.minorticks_on()

    #skads_bins, skads_counts, skads_count_errors = get_skads_counts_agnonly()
    #np.save('skads_agn.npy', np.array([skads_bins, skads_counts, skads_count_errors])) # AGN
    #ax1.plot(skads_bins*scaleF, skads_counts, color='k', linestyle = '--',label='Wilman et al. 2008 - AGN')
    #skads_bins, skads_counts, skads_count_errors = get_skads_counts_sfonly(majorlim=7) # SF
    #np.save('skads_sf.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
    #ax1.plot(skads_bins*scaleF, skads_counts, color='g', linestyle = ':',label='Wilman et al. 2008 - SF $<7$"')
    #skads_bins, skads_counts, skads_count_errors = get_skads_counts_sfonly(majorlim=7,less=False) # SF
    #np.save('skads_sf.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
    #ax1.plot(skads_bins*scaleF, skads_counts, color='b', linestyle = ':',label='Wilman et al. 2008 - SF $>7$"')
    #skads_bins, skads_counts, skads_count_errors = get_skads_counts_sfonly(majorlim=0,less=False) # SF
    #np.save('skads_sf.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
    #ax1.plot(skads_bins*scaleF, skads_counts, color='k', linestyle = ':',label='Wilman et al. 2008 - SF')
    
    skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[1,2],[1,2,3,4]], postprocess=True, nu=nu) # All
    np.save('skads_all.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
    ax1.plot(skads_bins*scaleF, skads_counts, color='gray', linestyle = '-',label='Wilman et al. 2008 - All pp')

    if 1:
        skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[1],[]], postprocess=False, nu=nu) # All
        np.save('skads_all.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
        ax1.errorbar(skads_bins*scaleF, skads_counts, skads_count_errors, color='blue', linestyle = '-',label='Wilman et al. 2008 - quiescent')
        skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[2],[]], postprocess=False, nu=nu) # All
        np.save('skads_all.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
        ax1.plot(skads_bins*scaleF, skads_counts, color='cyan', linestyle = '-',label='Wilman et al. 2008 - SB')
        skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[],[1]], postprocess=False, nu=nu) # All
        np.save('skads_all.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
        ax1.plot(skads_bins*scaleF, skads_counts, color='yellow', linestyle = '-',label='Wilman et al. 2008 - AGN RQ')
        skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[],[2]], postprocess=False, nu=nu) # All
        np.save('skads_all.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
        ax1.plot(skads_bins*scaleF, skads_counts, color='orange', linestyle = '-',label='Wilman et al. 2008 - AGN FRI')
        skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[],[3]], postprocess=False, nu=nu) # All
        np.save('skads_all.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
        ax1.plot(skads_bins*scaleF, skads_counts, color='red', linestyle = '-',label='Wilman et al. 2008 - AGN FRII')
        skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[],[4]], postprocess=False, nu=nu) # All
        np.save('skads_all.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
        ax1.plot(skads_bins*scaleF, skads_counts, color='green', linestyle = '-',label='Wilman et al. 2008 - AGN GPS')

    skads_bins, skads_counts, skads_count_errors = get_skads_counts(ttype=[[1,2],[1,2,3,4]], postprocess=False, nu=nu) # All
    np.save('skads_all.npy', np.array([skads_bins, skads_counts, skads_count_errors]))
    ax1.plot(skads_bins*scaleF, skads_counts, color='k', linestyle = '-',label='Wilman et al. 2008 - All')


    y1,y2 = ax1.get_ylim()
    #ax1.vlines(5.*noise.min()*scaleF, y1,y2, linestyles='dotted')
    #ax1.vlines(5.*noise.max()*scaleF, y1,y2, linestyles='dotted')

    #pl.title(' Source Counts ')
    pl.ylabel(r'$S^{5/2} dN / dS$ [Jy$^{3/2}$ sr$^{-1}$]')
    pl.xlabel(r'%i MHz flux [%s]' %(nu,unit))
    pl.rcParams['legend.fontsize'] ='medium'
    pl.legend(loc=4,numpoints=1,frameon=False)
    #pl.minorticks_on()
    #pl.savefig('%s_source_counts.png' %(name))
    #pl.savefig('%s_source_counts.eps' %(name))

    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax1.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    
    pp.set_attrib(ax1, xlim=(6e-1,1.2e4), ylim=(2e0,7e4))

    pp.fig_save_many(fig, 'beerze_plots/%s_source_counts' %(name))
