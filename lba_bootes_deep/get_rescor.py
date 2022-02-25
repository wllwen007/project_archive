import os
import astropy.wcs as pw
from astropy.table import Table,Column
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

import utils.plotting as pp
from scipy.optimize import leastsq

from math import *
import scipy.special

from matplotlib.ticker import MultipleLocator, LogFormatter, LogLocator, FuncFormatter



def flux_ratio_env_res (p,x,y):
  res = y - flux_ratio_env (p,x)
  return res
def flux_ratio_env (p,x):
  f =  p[0] + p[1]/(x**p[2])
  return f
def flux_ratio_env_inv (p,x):
  f =  p[0] - p[1]/(x**p[2])
  return f

    
    
    
    
def fit_env(snr, ffp, envlim=0.9):
    
    
    #SNRrange = np.logspace(np.log10(5), np.log10(snr.max()/2.), 15)
    SNRrange = np.array([6,7,8,10,12,15,20,30,40,60,80,100,200,500])
    env = np.ones(len(SNRrange)-1)*1.2
    #nsrci = np.sum(detect)
    for k in range(len(SNRrange)-1):
        SiSp = np.ma.masked_where(1-((snr>SNRrange[k])*(snr<SNRrange[k+1])),ffp).compressed()
        #print '#',SNRrange[k], SNRrange[k+1], len(SiSp)
        if len(SiSp) > 1:
            SiSpi = 10.0
            #print i,SiSpi, np.sum(SiSp<SiSpi), np.sum(SiSp<SiSpi)/len(SiSp)
            while float(np.sum(SiSp<SiSpi))/len(SiSp) > 1-envlim:
                #print i,SiSpi, np.sum(SiSp<SiSpi), np.sum(SiSp<SiSpi)/len(SiSp)
                #print k, SNRrange[k], len(SiSp),np.sum(SiSp<SiSpi)
                SiSpi-= 0.01
            env[k] = SiSpi
    #env[ k+1 ] = env[ k ]
    SNRrange = np.sqrt(SNRrange[1:]*SNRrange[:-1])
    
    p0 = [1, 1, 1]
    print('fitting {f:.0f} percent envelope'.format(f=100*envlim))
    pfit,cov,info,mesg,success = leastsq(flux_ratio_env_res, p0, args=(SNRrange, env), full_output=True)
    if success==1:
        # calculate final chi square
        chisq=sum(info["fvec"]*info["fvec"])
        dof=len(SNRrange)-len(p0)
        # chisq, sqrt(chisq/dof) agrees with gnuplot
        print("Converged with chi squared ",chisq)
        print("degrees of freedom, dof ", dof)
        print("RMS of residuals (i.e. sqrt(chisq/dof)) ", np.sqrt(chisq/dof))
        print("Reduced chisq (i.e. variance of residuals) ", chisq/dof)
        # uncertainties are calculated as per gnuplot, "fixing" the result
        # for non unit values of the reduced chisq.
        # values at min match gnuplot
        print("Fitted parameters at minimum, with 68% C.I.:")
        for i,pmin in enumerate(pfit):
            print("%2i %12f +/- %10f"%(i,pmin,np.sqrt(cov[i,i])*np.sqrt(chisq/dof)))
        print("Correlation matrix")
        # correlation matrix close to gnuplot
        for i in range(len(pfit)):
            print("%10i"%i, end=' ')
            for j in range(i+1):
                print("%10f"%(cov[i,j]/np.sqrt(cov[i,i]*cov[j,j]),), end=' ')
            print()
    else: 
        print("Not converged", mesg)
    #-----------------------------------------------
    s = ['%.3f' %(p) for p in pfit]
    print(', '.join(s))
    return SNRrange, env, pfit

def fit_env_up(snr, ffp, envlim=0.9):
    
    
    #SNRrange = np.logspace(np.log10(5), np.log10(snr.max()/2.), 15)
    SNRrange = np.array([4.5,5,5.5,6,6.5,7,7.5,8,9,10,11,12,15,20,30,40,60,80,100,200])
    env = np.ones(len(SNRrange)-1)*1.2
    #nsrci = np.sum(detect)
    for k in range(len(SNRrange)-1):
        SiSp = np.ma.masked_where(1-((snr>SNRrange[k])*(snr<SNRrange[k+1])),ffp).compressed()
        #print '#',SNRrange[k], SNRrange[k+1], len(SiSp)
        if len(SiSp) > 1:
            SiSpi = 1.0
            #print i,SiSpi, np.sum(SiSp<SiSpi), np.sum(SiSp<SiSpi)/len(SiSp)
            while float(np.sum(SiSp<SiSpi))/len(SiSp) <envlim:
                #print i,SiSpi, np.sum(SiSp<SiSpi), np.sum(SiSp<SiSpi)/len(SiSp)
                #print k, SNRrange[k], len(SiSp),np.sum(SiSp<SiSpi)
                SiSpi+= 0.01
            env[k] = SiSpi
    #env[ k+1 ] = env[ k ]
    SNRrange = np.sqrt(SNRrange[1:]*SNRrange[:-1])
    
    p0 = [1, 1, 1]
    print('fitting {f:.0f} percent envelope'.format(f=100*envlim))
    pfit,cov,info,mesg,success = leastsq(flux_ratio_env_res, p0, args=(SNRrange, env), full_output=True)
    if success==1:
        # calculate final chi square
        chisq=sum(info["fvec"]*info["fvec"])
        dof=len(SNRrange)-len(p0)
        # chisq, sqrt(chisq/dof) agrees with gnuplot
        print("Converged with chi squared ",chisq)
        print("degrees of freedom, dof ", dof)
        print("RMS of residuals (i.e. sqrt(chisq/dof)) ", np.sqrt(chisq/dof))
        print("Reduced chisq (i.e. variance of residuals) ", chisq/dof)
        # uncertainties are calculated as per gnuplot, "fixing" the result
        # for non unit values of the reduced chisq.
        # values at min match gnuplot
        print("Fitted parameters at minimum, with 68% C.I.:")
        for i,pmin in enumerate(pfit):
            print("%2i %12f +/- %10f"%(i,pmin,np.sqrt(cov[i,i])*np.sqrt(chisq/dof)))
        print("Correlation matrix")
        # correlation matrix close to gnuplot
        for i in range(len(pfit)):
            print("%10i"%i, end=' ')
            for j in range(i+1):
                print("%10f"%(cov[i,j]/np.sqrt(cov[i,i]*cov[j,j]),), end=' ')
            print()
    else: 
        print("Not converged", mesg)
    #-----------------------------------------------
    s = ['%.3f' %(p) for p in pfit]
    print(', '.join(s))
    return SNRrange, env, pfit



def mapvalue(fits, ra, dec):
    """determine if a point is inside a given coverage fits map
args
ra,dec - position coordinates (float)
    """
    
    # Parse the WCS keywords in the primary HDU

    #x,y,f,c = wcs.wcs_sky2pix(ra, dec,ra*0, ra*0, 1)
    
    hdulist = pf.open(fits)
    wcs = pw.WCS(hdulist[0].header)
    
    data = pf.getdata(fits)
    
    skycrd = np.array([ra,dec,0*ra,0*dec]).transpose()
    #skycrd = np.array([ra,dec]).transpose()
    pixcrd = wcs.wcs_sky2pix(skycrd, 0)
    x = pixcrd[:,0]
    y = pixcrd[:,1]
    x = np.array(x, dtype=int) #-1
    y = np.array(y, dtype=int) #-1
    
    nc,nf,ny,nx = data.shape
    
    


    return data[0,0,y,x]



scaleF = 1e3


LOFARcatfile = 'bootes_deep_lba_hbashift.cat.fits'

LOFARcatfileR = 'bootes_deep_lba_hbashift.cat.resolved.fits'



clobber=False
sr=3./3600.
firstbdsmcat = LOFARcatfile.replace('.fits','_first.fits')
#vo_url = 'http://vizier.u-strasbg.fr/viz-bin/votable/-A?-source=VIII/92&'
vo_url = 'VIII/92/first14'
if not os.path.exists(firstbdsmcat) or clobber:
    cmd = "/net/beerze/data2/wwilliams/software/bin/stilts  cdsskymatch cdstable='{url}' in={infits1} out={outfits} radius={radius:f} ra=RA dec=DEC find=best".format(url=vo_url, infits1=LOFARcatfile, outfits=firstbdsmcat, radius=sr)
    print(cmd)
    os.system(cmd)
    
#LOFARcatfile = firstbdsmcat


noiselim1 = 0.6   #in mJy
noiselim2 = 1.   #in mJy

Fscale = 1e3   ## scale factor to apply to cat fluxes plot in mJy
BMAJ = 15./3600.
BMIN = 15./3600.

LOFARcatall = Table.read(LOFARcatfile)
LOFARcat = Table.read(LOFARcatfile)
LOFARcat = LOFARcat[LOFARcat['Source_id'] >= 0]
LOFARcat = LOFARcat[LOFARcat['Peak_flux']/LOFARcat['Isl_rms'] >= 5]
LOFARcat_s = LOFARcat[LOFARcat['S_Code'] == 'S']


loss2 = 1

loss = 1
peakcor = LOFARcat['Peak_flux'] / loss

#peakcor = LOFARcat['Peak_flux'] / (loss*loss2)
#sel = peakcor/LOFARcat['Isl_rms'] > 5

#loss_sel = loss[sel]
#LOFARcat_sel = LOFARcat[sel]

#loss = loss[LOFARcat['S_Code=='S']
#LOFARcat= LOFARcat[LOFARcat['S_Code=='S']

nu = 54.
f = LOFARcat['Total_flux']
ef = LOFARcat['E_Total_flux']
fp = LOFARcat['Peak_flux']  
fpcor = LOFARcat['Peak_flux']   / (loss)
efp = LOFARcat['E_Peak_flux']
rms = LOFARcat['Isl_rms']
ffp = f/fp
effp = np.sqrt((efp/fp)**2.+(ef/f)**2.)
ffpcor = f/fpcor
effpcor = np.sqrt((efp/fpcor)**2.+(ef/f)**2.)
snr = fp/rms


fall = LOFARcatall['Total_flux']
efall = LOFARcatall['E_Total_flux']
fpall = LOFARcatall['Peak_flux']  
fpcorall = LOFARcatall['Peak_flux']   / (loss)
efpall = LOFARcatall['E_Peak_flux']
rmsall = LOFARcatall['Isl_rms']
ffpall = fall/fpall
effpall = np.sqrt((efpall/fpall)**2.+(efall/fall)**2.)
ffpcorall = fall/fpcorall
effpcorall = np.sqrt((efpall/fpcorall)**2.+(efall/fall)**2.)
snrall = fpall/rmsall



#rad_f = angular_separation(LOFARcat_f['RA'], LOFARcat_f['DEC'],RA0,DEC0)
#bloss_f = bandwidth_smearing2(97.5e3,150.e6,7.,rad_f*3600.)
#tloss_f = time_smearing2(8, 7.,rad_f*3600)
#loss_f = bloss_f*tloss_f*loss2
#f_f = LOFARcat_f.Total_flux
#ef_f = LOFARcat_f.E_Total_flux
#fp_f = LOFARcat_f.Peak_flux  
#fpcor_f = LOFARcat_f.Peak_flux   / (loss_f)
#efp_f = LOFARcat_f.E_Peak_flux
#rms_f = LOFARcat_f.Isl_rms
#ffpcor_f = f_f/fpcor_f
#snr_f = fp_f/rms_f


#snrcor = fpcor/rms
#snr2 = fp/rms2



def S14(S60, alpha=-0.8, freq=54.):
    # S60 in Jy, S14 in mJy
    return 1000*S60*(1400./freq)**alpha

def phimed(S14):
    # S14 in mJy
    return 2*(S14)**0.3

def phimed2(S14):
    # S14 in mJy
    return 2.*(S14 < 1.) +  (2*(S14)**0.3)*(S14>=1.)

def h(phi,S60):                  
    S = S14(S60)
    #return 2**(-1*(phi/phimed(S))**0.62)
    return np.exp(-1*np.log(2)*(phi/phimed(S))**0.62)

def h2(phi,S60):                  
    S = S14(S60)
    #return 2**(-1*(phi/phimed(S))**0.62)
    return np.exp(-1*np.log(2)*(phi/phimed2(S))**0.62)

#def h(phi,S):    
    ##return 2**(-1*(phi/phimed(S))**0.62)
    #return np.exp(-np.log(2)*(phi/phimed(S))**0.62)

#def h2(theta):
    #f = (1./1.6**theta)*(theta<4) +  (theta**-1.3 -0.01)*(theta>=4)
    #f[f<0] = 0
    #return f




if 0:
    
    fig,ax = pp.paper_single_ax()
    #ax.loglog(snr, ffp, 'k.', alpha=0.1)
    ax.scatter(snrall, ffpcorall,  marker='o',s=1,c='gray')
    ax.scatter(snr, ffpcor, marker='o',s=7,c='C0')
    ax.loglog()
    
    #SNRrange, env, pfit = fit_env(snr, ffpcor, envlim=0.66)
    #ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env(snr, ffpcor, envlim=0.90)
    #ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env(snr, ffpcor, envlim=0.95)
    #SNRrange1, env1, pfit1 = fit_env(snr, ffpcor, envlim=0.84)
    sel = (LOFARcat['Maj']*3600. < 25.)
    ax.scatter(snr[sel], ffpcor[sel], marker='o',s=7,c='C1')
    SNRrange, env, pfit = fit_env(snr[sel], ffpcor[sel], envlim=0.05)
    #ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env_up(snr, ffpcor, envlim=0.95)
    
    ax.plot(SNRrange, env, 'C2', marker='o',linestyle='none')
    
    SNRrange_env = np.logspace(np.log10(1), np.log10(2000), 100)
    
    envelope = flux_ratio_env (pfit, SNRrange_env )
    envelope_inv = flux_ratio_env_inv (pfit, SNRrange_env )
    #envelope_inv2 = flux_ratio_env ([pfit1[0],3*pfit1[1],pfit1[0]], SNRrange_env )
    
    ax.semilogx(SNRrange_env, np.ones_like(SNRrange_env),'k')
    ax.semilogx(SNRrange_env, envelope,'C2')
    ax.semilogx(SNRrange_env, envelope_inv,'C1')
    #ax.semilogx(SNRrange_env, envelope_inv2,'C3')

    
    

    x1,x2 = ax.get_xlim()
    ax.set_xlim(1, 2000)
    ax.set_ylim(0.5, 30)
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xlabel(r'$S_{\rm peak} / \sigma_{\rm local}$')
    ax.set_ylabel(r'$S_{\rm total}/S_{\rm peak}$')

        
    pp.fig_save_many(fig,'beerze_plots/total_peak_fit_new')  
    
    Res = (ffpcorall > flux_ratio_env (pfit, snrall ))
    LOFARcatall.add_column(Column(name='Resolved',data=Res))
    LOFARcatall.write(LOFARcatfileR,overwrite=True)  
    
if 1:
    
    fig,ax = pp.paper_single_ax()
    #ax.loglog(snr, ffp, 'k.', alpha=0.1)
    ax.scatter(snrall, ffpcorall,  marker='o',s=1,c='gray')
    ax.scatter(snr, ffpcor, marker='o',s=7,c='C0')
    ax.loglog()
    
    #SNRrange, env, pfit = fit_env(snr, ffpcor, envlim=0.66)
    #ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env(snr, ffpcor, envlim=0.90)
    #ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env(snr, ffpcor, envlim=0.95)
    #SNRrange1, env1, pfit1 = fit_env(snr, ffpcor, envlim=0.84)
    SNRrange, env, pfit = fit_env(snr, ffpcor, envlim=0.9875)
    #ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env_up(snr, ffpcor, envlim=0.95)
    
    ax.plot(SNRrange, env, 'C2', marker='o',linestyle='none')
    
    SNRrange_env = np.logspace(np.log10(1), np.log10(2000), 100)
    
    envelope = flux_ratio_env (pfit, SNRrange_env )
    envelope_inv = flux_ratio_env_inv (pfit, SNRrange_env )
    envelope_inv2 = flux_ratio_env ([2.0,-1*pfit[1],pfit[2]], SNRrange_env )
    
    ax.semilogx(SNRrange_env, np.ones_like(SNRrange_env),'k')
    ax.semilogx(SNRrange_env, envelope,'C2')
    #ax.semilogx(SNRrange_env, envelope_inv,'C1')
    ax.semilogx(SNRrange_env, envelope_inv2,'C3')

    x1,x2 = ax.get_xlim()
    ax.set_xlim(1, 2000)
    ax.set_ylim(0.5, 30)
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xlabel(r'$S_{\rm peak} / \sigma_{\rm local}$')
    ax.set_ylabel(r'$S_{\rm total}/S_{\rm peak}$')

        
    pp.fig_save_many(fig,'beerze_plots/total_peak_fit')  
    Res = (ffpcorall > flux_ratio_env ([2.0,-1*pfit[1],pfit[2]], snrall ))
    '''
    
    fig,ax = pp.paper_single_ax()
    #ax.loglog(snr, ffp, 'k.', alpha=0.1)
    ax.loglog(snr_f, ffpcor_f, 'k.', alpha=0.1)
    
    SNRrange, env, pfit = fit_env(snr_f, ffpcor_f, envlim=0.90)
    ax.plot(SNRrange, env, 'g.')
    SNRrange, env, pfit = fit_env(snr_f, ffpcor_f, envlim=0.95)
    ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env(snr, ffpcor, envlim=0.99)
    #ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env_up(snr, ffpcor, envlim=0.95)
    
    ax.plot(SNRrange, env, 'g.')
    
    SNRrange_env = np.logspace(np.log10(5), np.log10(snr.max()), 100)
    
    envelope = flux_ratio_env (pfit, SNRrange_env )
    envelope_inv = flux_ratio_env_inv (pfit, SNRrange_env )
    
    ax.semilogx(SNRrange_env, np.ones_like(SNRrange_env),'k:')
    ax.semilogx(SNRrange_env, envelope,'g-')
    ax.semilogx(SNRrange_env, envelope_inv,'g-')

    
    
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

    x1,x2 = ax.get_xlim()
    ax.set_xlim(3, x2)
    ax.set_ylim(0.5, 30)
    ax.set_xlabel(r'$S_{\rm peak} / \sigma$')
    ax.set_ylabel(r'$S_{\rm int}/S_{\rm peak}$')

        
    pp.fig_save_many(fig,'beerze_plots/total_peak_fit_f')  
    '''
    
if 0:
    fig,ax = pp.paper_single_ax()
    #ax.loglog(snr, ffp, 'k.', alpha=0.1)
    ax.loglog(snr_s, ffpcor_s, 'k.', alpha=0.1)
    
    SNRrange, env, pfit = fit_env(snr_s, ffpcor_s, envlim=0.90)
    ax.plot(SNRrange, env, 'g.')
    SNRrange, env, pfit = fit_env(snr_s, ffpcor_s, envlim=0.95)
    ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env(snr, ffpcor, envlim=0.99)
    #ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env_up(snr, ffpcor, envlim=0.95)
    
    ax.plot(SNRrange, env, 'g.')
    
    SNRrange_env = np.logspace(np.log10(5), np.log10(snr.max()), 100)
    
    envelope = flux_ratio_env (pfit, SNRrange_env )
    envelope_inv = flux_ratio_env_inv (pfit, SNRrange_env )
    
    ax.semilogx(SNRrange_env, np.ones_like(SNRrange_env),'k:')
    ax.semilogx(SNRrange_env, envelope,'g-')
    ax.semilogx(SNRrange_env, envelope_inv,'g-')

    
    
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

    x1,x2 = ax.get_xlim()
    ax.set_xlim(3, x2)
    ax.set_ylim(0.5, 30)
    ax.set_xlabel(r'$S_{\rm peak} / \sigma$')
    ax.set_ylabel(r'$S_{\rm int}/S_{\rm peak}$')

        
    pp.fig_save_many(fig,'beerze_plots/total_peak_fit_s')  
    
    fig,ax = pp.paper_single_ax()
    #ax.loglog(snr, ffp, 'k.', alpha=0.1)
    ax.loglog(snr_g, ffpcor_g, 'k.', alpha=0.1)
    
    SNRrange, env, pfit = fit_env(snr_g, ffpcor_g, envlim=0.90)
    ax.plot(SNRrange, env, 'g.')
    SNRrange, env, pfit = fit_env(snr_g, ffpcor_g, envlim=0.95)
    ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env(snr, ffpcor, envlim=0.99)
    #ax.plot(SNRrange, env, 'g.')
    #SNRrange, env, pfit = fit_env_up(snr, ffpcor, envlim=0.95)
    
    ax.plot(SNRrange, env, 'g.')
    
    SNRrange_env = np.logspace(np.log10(5), np.log10(snr.max()), 100)
    
    envelope = flux_ratio_env (pfit, SNRrange_env )
    envelope_inv = flux_ratio_env_inv (pfit, SNRrange_env )
    
    ax.semilogx(SNRrange_env, np.ones_like(SNRrange_env),'k:')
    ax.semilogx(SNRrange_env, envelope,'g-')
    ax.semilogx(SNRrange_env, envelope_inv,'g-')

    
    
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

    x1,x2 = ax.get_xlim()
    ax.set_xlim(3, x2)
    ax.set_ylim(0.5, 30)
    ax.set_xlabel(r'$S_{\rm peak} / \sigma$')
    ax.set_ylabel(r'$S_{\rm int}/S_{\rm peak}$')

        
    pp.fig_save_many(fig,'beerze_plots/total_peak_fit_g')  
    
    #facets = np.unique(LOFARcat['Facet_id'])
    #for i in range(len(facets)):
        
        #ind = (LOFARcat['Facet_id']==facets[i])
        
        #SNRrange, env, pfit = fit_env(snr[ind], ffpcor[ind], envlim=0.95)
        
        #fig,ax = pp.paper_single_ax()
        ##ax.loglog(snr, ffp, 'k.', alpha=0.1)
        #ax.loglog(snr[ind], ffpcor[ind], 'k.', alpha=0.1)
        #ax.plot(SNRrange, env, 'g.')
        
        #SNRrange_env = np.logspace(np.log10(5), np.log10(snr.max()), 100)
        
        #envelope = flux_ratio_env (pfit, SNRrange_env )
        #envelope_inv = flux_ratio_env_inv (pfit, SNRrange_env )
        
        #ax.semilogx(SNRrange_env, np.ones_like(SNRrange_env),'k:')
        #ax.semilogx(SNRrange_env, envelope,'g-')
        #ax.semilogx(SNRrange_env, envelope_inv,'g-')

        
        
        #ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        #ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

        #x1,x2 = ax.get_xlim()
        #ax.set_xlim(3, x2)
        #ax.set_ylim(0.5, 30)
        #ax.set_xlabel(r'$S_{\rm peak} / \sigma$')
        #ax.set_ylabel(r'$S_{\rm int}/S_{\rm peak}$')
        #ax.set_title("s{f}".format(f=facets[i]))

            
        #pp.fig_save_many(fig,'beerze_plots/total_peak_fit_s{f:.0f}'.format(f=facets[i])) 

    
    #rad_bins = [0.,0.5,1.,1.5,2.,3.]
    #for i in range(len(rad_bins)-1):
        
        #ind = (rad>=rad_bins[i]) & (rad<rad_bins[i+1])
        
        #SNRrange, env, pfit = fit_env(snr[ind], ffpcor[ind], envlim=0.95)
        
        #fig,ax = pp.paper_single_ax()
        ##ax.loglog(snr, ffp, 'k.', alpha=0.1)
        #ax.loglog(snr[ind], ffpcor[ind], 'k.', alpha=0.1)
        #ax.plot(SNRrange, env, 'g.')
        
        #SNRrange_env = np.logspace(np.log10(5), np.log10(snr.max()), 100)
        
        #envelope = flux_ratio_env (pfit, SNRrange_env )
        #envelope_inv = flux_ratio_env_inv (pfit, SNRrange_env )
        
        #ax.semilogx(SNRrange_env, np.ones_like(SNRrange_env),'k:')
        #ax.semilogx(SNRrange_env, envelope,'g-')
        #ax.semilogx(SNRrange_env, envelope_inv,'g-')

        
        
        #ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        #ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

        #x1,x2 = ax.get_xlim()
        #ax.set_xlim(3, x2)
        #ax.set_ylim(0.5, 30)
        #ax.set_xlabel(r'$S_{\rm peak} / \sigma$')
        #ax.set_ylabel(r'$S_{\rm int}/S_{\rm peak}$')
        #ax.set_title("${r1:.1f} < R < {r2:.1f}$ deg".format(r1=rad_bins[i],r2=rad_bins[i+1]))

            
        #pp.fig_save_many(fig,'beerze_plots/total_peak_fit_{i:d}'.format(i=i))  
        
fig,ax = pp.paper_single_ax()   
n,b,p = ax.hist(LOFARcat['Maj']*3600., bins=50, histtype='step')
y1,y2 = ax.get_ylim()
ax.vlines(7.4,y1,y2)
ax.set_ylim(y1,y2)
pp.fig_save_many(fig,'beerze_plots/total_maj_dist_r')  
        
fig,ax = pp.paper_single_ax()   
n,b,p = ax.hist(LOFARcat['Min']*3600., bins=50, histtype='step')
y1,y2 = ax.get_ylim()
ax.vlines(5.6,y1,y2)
ax.set_ylim(y1,y2)
pp.fig_save_many(fig,'beerze_plots/total_min_dist_r')  


def snr_env2(snr, a, b,c):
    return a/snr**c + b

# get unresolved from envelope
unres_ind = (ffpcor <= flux_ratio_env_inv(pfit, snr))
unres_ind2 = (ffpcor <= snr_env2(snr,50,1.05,1.7))
res_ind = ~unres_ind

# SIZE - FLUX
if 1:
    fig,ax = pp.paper_single_ax()
    #ax.set_yscale('log')
    ax.set_xscale('log')
    #size = np.sqrt(LOFARcat['Maj']*LOFARcat['Min'])
    size = np.sqrt(LOFARcat['DC_Maj']*LOFARcat['DC_Min'])
    ax.plot(f[res_ind]*Fscale, size[res_ind]*3600., 'k.',alpha=0.1)
    ax.plot(f[unres_ind]*Fscale, size[unres_ind]*3600., 'k.',alpha=0.1)

    x1,x2 = ax.get_xlim()
    srange = np.logspace(np.log10(x1),np.log10(x2), 50)/Fscale   # in Jy
    ax.plot(srange*Fscale, phimed(S14(srange/Fscale)), 'r')
    #ax.plot(srange*Fscale, phimed(S14(srange,alpha=-0.5)), 'r:')
    l, = ax.plot(srange*Fscale, phimed2(S14(srange)), 'r', label='Median angular size')
    #ax.plot(srange*Fscale, phimed2(S14(srange,alpha=-0.5)), 'r:')
    #lim1 = 3600*(BMAJ*BMIN*srange/(5*noiselim1))/BMAJ
    #lim2 = 3600*(BMAJ*BMIN*srange/(5*noiselim1))/BMIN
    #minbeamsize = np.min(LOFARcat['Facet_BMAJ']*LOFARcat['Facet_BMIN'])
    #maxbeamsize = np.max(LOFARcat['Facet_BMAJ']*LOFARcat['Facet_BMIN'])
    minbeamsize = np.min(BMAJ*BMIN)
    maxbeamsize = np.max(BMAJ*BMIN)
    limA1a = 3600*np.sqrt(minbeamsize*srange*Fscale/(5*noiselim1))
    limA1b = 3600*np.sqrt(maxbeamsize*srange*Fscale/(5*noiselim1))
    limA1 = 3600*np.sqrt(BMAJ*BMIN*srange*Fscale/(5*noiselim1))
    limA2 = 3600*np.sqrt(BMAJ*BMIN*srange*Fscale/(5*noiselim2))
    limB1 = 3600*np.sqrt(BMAJ*BMIN*srange*Fscale/(5*noiselim2))
    
    meanbeam = 3600.*sqrt(BMAJ*BMIN)
    DC_limA1 = np.sqrt(limA1**2. - meanbeam**2.)
    DC_limA2 = np.sqrt(limA2**2. - meanbeam**2.)
    #limA = np.max((lim1,lim2,lim3), axis=0)
    #ax.fill_between(srange*Fscale, limA1, limA2, color='b', alpha=0.2)
    ax.fill_between(srange*Fscale, DC_limA1, DC_limA2, color='b', alpha=0.2, label='Maximum size')

    pb = plt.Rectangle((0, 0), 1, 1, fc='b', ec='b',  alpha=0.2)
    ax.legend([l, pb], ['Median angular size', 'Maximum size'],loc='upper left')
    
    #ax.plot(srange*Fscale, limA1a, 'b')
    #ax.plot(srange*Fscale, limA1b, 'b')
    #ax.plot(srange*Fscale, lim2, 'b')
    #ax.legend()

    y1,y2 = ax.get_ylim()
    #ax.semilogx(f150_lim_range,alpha_lim1, ':' , c='gray' )
    #ax.semilogx(f150_lim_range,np.ones_like(f150_lim_range)*mean_alpha, '-', c='gray' )
    ax.set_ylim(0,100)
    ax.set_xlim(3, 500)
    ax.set_xlabel(r'$S_{\rm int}$ [mJy]')
    ax.set_ylabel(r'$ \theta$ [arcsec]')
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

    pp.fig_save_many(fig, 'beerze_plots/size_flux')

    #sys.exit()

    fig,ax = pp.paper_single_ax()
    ax.set_xscale('log')
    thetalim1 = DC_limA1
    thetalim2 = DC_limA2
    ax.fill_between(srange*Fscale, h(thetalim1, srange), h(thetalim2, srange), color='b', alpha=0.35, label='variable')
    ax.fill_between(srange*Fscale, h2(thetalim1, srange), h2(thetalim2, srange), color='g', edgecolor='g', facecolor='none', hatch='X', label='fixed')
    pb = plt.Rectangle((0, 0), 1, 1, edgecolor='g', facecolor='none',  hatch='X')
    pg = plt.Rectangle((0, 0), 1, 1, facecolor='b', edgecolor='b', alpha=0.35)
    ax.legend([pg, pb], [r'$\theta_{\mathrm{med}}$ fixed', r'$\theta_{\mathrm{med}}$ variable'])
    #ax.legend()
    ax.set_xlabel(r'$S_{\rm int}$ [mJy]')
    ax.set_ylabel(r'$h(>\theta_{\mathrm{max}})$')
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    pp.fig_save_many(fig, 'beerze_plots/res_cor_h')
    
    
    fig,ax = pp.paper_single_ax()
    ax.set_xscale('log')
    thetalim1 = DC_limA1
    thetalim2 = DC_limA2
    ax.fill_between(srange*Fscale, 1/(1.-h(thetalim1, srange)), 1./(1-h(thetalim2, srange)), color='b', alpha=0.35, label='variable')
    ax.fill_between(srange*Fscale, 1/(1.-h2(thetalim1, srange)), 1/(1.-h2(thetalim2, srange)), color='g', edgecolor='g', facecolor='none', hatch='X',  label='fixed')
    pb = plt.Rectangle((0, 0), 1, 1, edgecolor='g', facecolor='none', hatch='X')
    pg = plt.Rectangle((0, 0), 1, 1, facecolor='b', edgecolor='b', alpha=0.35)
    ax.legend([pg, pb], [r'$\theta_{\mathrm{med}}$ fixed', r'$\theta_{\mathrm{med}}$ variable'])
    #ax.legend()
    ax.set_xlabel(r'$S_{\rm int}$ [mJy]')
    ax.set_ylabel(r'$1/[1-h(>\theta_{\mathrm{max}})]$')
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    pp.fig_save_many(fig, 'beerze_plots/res_cor')
    
    RC1low = 1/(1.-h(thetalim1, srange))
    RC1high = 1./(1-h(thetalim2, srange))
    RC2low = 1/(1.-h2(thetalim1, srange))
    RC2high = 1/(1.-h2(thetalim2, srange))
     
    
    
    np.save('rescor.npy', np.array((srange, RC1low, RC1high, RC2low, RC2high)))


def get_res_cor(S):
    srange, RC1low, RC1high, RC2low, RC2high = np.load('rescor.npy')
    si = sum(srange<S)
    return RC1low[si], RC1high[si], RC2low[si], RC2high[si]

# SIZE - PEAK
if 1:
    fig,ax = pp.paper_single_ax()
    #ax.set_yscale('log')
    ax.set_xscale('log')
    #size = np.sqrt(LOFARcat['Maj']*LOFARcat['Maj'])  # geometric mean
    size = np.sqrt(LOFARcat['DC_Maj']*LOFARcat['DC_Maj'])  # geometric mean
    ax.plot(fp*Fscale, size*3600., 'k.',alpha=0.1)

    x1,x2 = ax.get_xlim()
    peakrange = np.logspace(np.log10(x1),np.log10(x2), 50)/Fscale
    ax.plot(peakrange*Fscale, phimed(S14(peakrange)), 'r')
    ax.plot(peakrange*Fscale, phimed(S14(peakrange,alpha=-0.5)), 'r:')
    lim1 = 3600*np.sqrt(BMAJ*BMIN*5*noiselim1/(peakrange*Fscale))
    ax.plot(peakrange*Fscale, lim1, 'b')

    y1,y2 = ax.get_ylim()
    #ax.semilogx(f150_lim_range,alpha_lim1, ':' , c='gray' )
    #ax.semilogx(f150_lim_range,np.ones_like(f150_lim_range)*mean_alpha, '-', c='gray' )
    ax.set_ylim(0,100)
    ax.set_xlim(2, 500)
    ax.set_xlabel(r'$S_{\rm peak}$ [mJy]')
    ax.set_ylabel(r'$ \theta$ [arcsec]')
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

    pp.fig_save_many(fig, 'beerze_plots/size_peak')


fig,ax = pp.paper_single_ax()
Size = LOFARcat['Maj'][(f*Fscale > 1) & (f*Fscale < 5)]*3600.
#ax.set_yscale('log')
#ax.set_xscale('log')
#ax.plot(f*Fscale, , 'k.',alpha=0.1)
Size = LOFARcat['Maj'][(f*Fscale > 1) & (f*Fscale <= 2.5)]*3600.
ax.hist(Size, density=True, cumulative=-1, histtype='step')
Size = LOFARcat['Maj'][(f*Fscale > 2.5) & (f*Fscale <= 5)]*3600.
ax.hist(Size, density=True, cumulative=-1, histtype='step')
Size = LOFARcat['Maj'][(f*Fscale > 5) & (f*Fscale <= 7.5)]*3600.
ax.hist(Size, density=True, cumulative=-1, histtype='step')
Size = LOFARcat['Maj'][(f*Fscale > 7.5) & (f*Fscale <= 10)]*3600.
ax.hist(Size, density=True, cumulative=-1, histtype='step')
x1,x2 = ax.get_xlim()
y1,y2 = ax.get_ylim()
#ax.semilogx(f150_lim_range,alpha_lim1, ':' , c='gray' )
#ax.semilogx(f150_lim_range,np.ones_like(f150_lim_range)*mean_alpha, '-', c='gray' )
#ax.set_ylim(y1,y2)
#ax.set_xlim(3, x2)
ax.set_ylabel(r'$ h(>\theta_{a})$')
ax.set_xlabel(r'$ \theta_{a}$ [arcsec]')
#ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

pp.fig_save_many(fig, 'beerze_plots/size_distrib_1mJy_5mJy')


#fig,ax = pp.paper_single_ax()
#ax.set_yscale('log')
#ax.set_xscale('log')
#ax.plot(snr_s, ffpcor_s, 'k.',alpha=0.1)
#snrange = np.logspace(0,3,50)
#def snr_env(snr, a, beta):
    #return a/(snr**beta)
#def snr_env2(snr, a, b,c):
    #return a/snr**c + b
##ax.plot(snrange, 1+snr_env(snrange,25,2.5), 'r')
##ax.plot(snrange, 1-snr_env(snrange,25,2.5), 'r')
##ax.plot(snrange, snr_env2(snrange,10,1.05,1.3), 'g')
##ax.plot(snrange, snr_env2(snrange,10,1.05,1.7), 'g:')
##ax.plot(snrange, snr_env2(snrange,50,1.15,1.7), 'b:')
#ax.plot(snrange, snr_env2(snrange,250,1.25,2.7), 'b:')
##ax.plot(snrange, 1-snr_env2(snrange,25,2.5), 'g')
#y1,y2 = ax.get_ylim()
##ax.semilogx(f150_lim_range,alpha_lim1, ':' , c='gray' )
##ax.semilogx(f150_lim_range,np.ones_like(f150_lim_range)*mean_alpha, '-', c='gray' )
##ax.set_ylim(y1,y2)
#x1,x2 = ax.get_xlim()
#ax.set_xlim(3, x2)
#ax.set_ylim(0.5, 30)
#ax.set_xlabel(r'$S_{\rm peak} / \sigma$')
#ax.set_ylabel(r'$S_{\rm int}/S_{\rm peak}$')
#ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
#ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
#pp.fig_save_many(fig, 'beerze_plots/total_peak_snr_s')


#fig,ax = pp.paper_single_ax()
#ax.set_yscale('log')
#ax.set_xscale('log')
#ax.plot(snr_g, ffpcor_g, 'k.',alpha=0.1)
#snrange = np.logspace(0,3,50)
#def snr_env(snr, a, beta):
    #return a/(snr**beta)
#def snr_env2(snr, a, b,c):
    #return a/snr**c + b
##ax.plot(snrange, 1+snr_env(snrange,25,2.5), 'r')
##ax.plot(snrange, 1-snr_env(snrange,25,2.5), 'r')
##ax.plot(snrange, snr_env2(snrange,10,1.05,1.3), 'g')
##ax.plot(snrange, snr_env2(snrange,10,1.05,1.7), 'g:')
##ax.plot(snrange, snr_env2(snrange,50,1.15,1.7), 'b:')
#ax.plot(snrange, snr_env2(snrange,250,1.25,2.7), 'b:')
##ax.plot(snrange, 1-snr_env2(snrange,25,2.5), 'g')
#y1,y2 = ax.get_ylim()
##ax.semilogx(f150_lim_range,alpha_lim1, ':' , c='gray' )
##ax.semilogx(f150_lim_range,np.ones_like(f150_lim_range)*mean_alpha, '-', c='gray' )
##ax.set_ylim(y1,y2)
#x1,x2 = ax.get_xlim()
#ax.set_xlim(3, x2)
#ax.set_ylim(0.5, 30)
#ax.set_xlabel(r'$S_{\rm peak} / \sigma$')
#ax.set_ylabel(r'$S_{\rm int}/S_{\rm peak}$')
#ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
#ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
#pp.fig_save_many(fig, 'beerze_plots/total_peak_snr_g')


fig,ax = pp.paper_single_ax()
ax.set_yscale('log')
ax.set_xscale('log')
ax.plot(snr, ffpcor, 'k.',alpha=0.1)
snrange = np.logspace(0,4,50)
def snr_env(snr, a, beta):
    return a/(snr**beta)
def snr_env2(snr, a, b,c):
    return a/snr**c + b
def snr_env2p(p, x):
    return p[1]/x**p[2] + p[0]
#ax.plot(snrange, 1+snr_env(snrange,25,2.5), 'r')
#ax.plot(snrange, 1-snr_env(snrange,25,2.5), 'r')
#ax.plot(snrange, snr_env2(snrange,10,1.05,1.3), 'g')
#ax.plot(snrange, snr_env2(snrange,10,1.05,1.7), 'g:')
#ax.plot(snrange, snr_env2(snrange,50,1.15,1.7), 'b:')
#ax.plot(snrange, snr_env2(snrange,250,1.25,2.7), 'r', ls='dashed', lw=2, label='Envelope')
#ax.plot(snrange, snr_env2(snrange,150,1.1,2.7), 'r', ls='dashed', lw=2, label='Envelope')

#ax.plot(snrange, snr_env2p([1.02, 14.36, 1.55],snrange), 'r', ls='dashed', lw=2, label='Fitted envelope')
ax.plot(snrange, snr_env2p([1.055, 74.6, 2.02],snrange), 'r', ls='dashed', lw=2, label='Fitted envelope')

ax.legend()

#ax.plot(snrange, 1-snr_env2(snrange,25,2.5), 'g')
y1,y2 = ax.get_ylim()
#ax.semilogx(f150_lim_range,alpha_lim1, ':' , c='gray' )
#ax.semilogx(f150_lim_range,np.ones_like(f150_lim_range)*mean_alpha, '-', c='gray' )
#ax.set_ylim(y1,y2)
x1,x2 = ax.get_xlim()
ax.xaxis.set_major_locator(LogLocator(subs=[1,2,5]))
ax.set_xlim(4.5, 1200)
#ax.set_xlim(3, 2000)
ax.set_ylim(0.51, 25.)
ax.set_xlabel(r'$S_{\rm peak} / \sigma$')
ax.set_ylabel(r'$S_{\rm int}/S_{\rm peak}$')
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pp.fig_save_many(fig, 'beerze_plots/total_peak_snr')



#fig,ax = pp.paper_single_ax()
#ax.set_yscale('log')
#ax.set_xscale('log')
#ax.plot(snr, (LOFARcat['Maj'] * LOFARcat['Min'])/(BMAJ*BMIN), 'k.',alpha=0.1)
#y1,y2 = ax.get_ylim()
##ax.semilogx(f150_lim_range,alpha_lim1, ':' , c='gray' )
##ax.semilogx(f150_lim_range,np.ones_like(f150_lim_range)*mean_alpha, '-', c='gray' )
##ax.set_ylim(y1,y2)
#x1,x2 = ax.get_xlim()
#ax.set_xlim(3, x2)
#ax.set_xlabel(r'$S_{\rm peak} / \sigma$')
#ax.set_ylabel(r'$\theta_{a} \theta_{b}/ b_{a} b_{b}$')
#pp.fig_save_many(fig, 'beerze_plots/source_size_area_snr')


fig,ax = pp.paper_single_ax()
ax.set_yscale('log')
ax.set_xscale('log')
ax.plot(fp*scaleF, f*scaleF, 'k.',alpha=0.1)
y1,y2 = ax.get_ylim()
#ax.semilogx(f150_lim_range,alpha_lim1, ':' , c='gray' )
#ax.semilogx(f150_lim_range,np.ones_like(f150_lim_range)*mean_alpha, '-', c='gray' )
#ax.set_ylim(y1,y2)
x1,x2 = ax.get_xlim()
#ax.set_xlim(3, x2)
ax.set_xlabel(r'$S_{\rm peak}$ [mJy beam$^{-1}$]')
ax.set_ylabel(r'$S_{\rm int}$ [mJy]')
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pp.fig_save_many(fig, 'beerze_plots/total_peak')

fig,ax = pp.paper_single_ax()
ax.set_yscale('log')
ax.set_xscale('log')
ax.plot(snr, f*scaleF, 'k.',alpha=0.1)
y1,y2 = ax.get_ylim()
#ax.semilogx(f150_lim_range,alpha_lim1, ':' , c='gray' )
#ax.semilogx(f150_lim_range,np.ones_like(f150_lim_range)*mean_alpha, '-', c='gray' )
#ax.set_ylim(y1,y2)
x1,x2 = ax.get_xlim()
#ax.set_xlim(3, x2)
ax.set_xlabel(r'$S_{\rm peak} / \sigma$')
ax.set_ylabel(r'$S_{\rm int}$ [mJy]')
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pp.fig_save_many(fig, 'beerze_plots/total_snr')









####################################

#def flux_ratio_env_res (p,x,y):
  #res = y - (1 + p[2]) + (p[0]**2. + (p[1]/x)**2.)**0.5
  #return res
#def flux_ratio_env (p,x):
  #f = (1+p[2]) - (p[0]**2. + (p[1]/x)**2.)**0.5
  #return f
#def flux_ratio_env_inv (p,x):
  #f = (1+p[2]) + (p[0]**2. + (p[1]/x)**2.)**0.5
  #return f
  
