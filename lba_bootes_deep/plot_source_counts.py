import matplotlib as mpl
#mpl.use('agg')
import utils.plotting as pp
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table, Column
import astropy.coordinates as ac
import astropy.units as u

def get_res_cor(S):
    srange, RC1low, RC1high, RC2low, RC2high = np.load('/net/lofar2//data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/rescor.npy')
    si = sum(srange<S)
    if si == len(srange):
        si -= 1
    return np.array([RC1low[si], RC1high[si], RC2low[si], RC2high[si]])


import configparser
config = configparser.ConfigParser()
config.read('image.config')
plotpath = config['SETUP']['plotpath']
deeptable = config['deeptable']['deeptable']
deepids = config['deeptable']['deepids']
cattable = config['srccnts']['cattable']
rmstable = config['srccnts']['rmstable']
compltable = config['srccnts']['compltable']


deephba = Table.read(deeptable)
sF = (60./150.)**-0.6
deephba['Total_flux'] = deephba['Total_flux']*sF
deephba['Peak_flux'] = deephba['Peak_flux']*sF
deephba['DC_Maj'] = deephba['DC_Maj']*3600.
deephba.add_column(Column(data=np.zeros(len(deephba)),name='Size'))
deephba['Size'][deephba['S_Code']!='Z'] = deephba['DC_Maj'][deephba['S_Code']!='Z']
deephba['Size'][deephba['S_Code']=='Z'] = deephba['LGZ_Size'][deephba['S_Code']=='Z']

#cat = Table.read('bootes_deep_lba.cat.mergedeblend.fits')
#cat = Table.read('bootes_deep_lba.gcat.fits')
cat = Table.read(cattable)
trms = Table.read(rmstable)
trms['rms'] = 1e-3*trms['rms']  # put in Jy


trms_starmask = Table.read('bootes_deep_lba_hbashift.Area_rms_starmask.fits')
trms_rin = Table.read('bootes_deep_lba_hbashift.Area_rms_rin2.fits')
trms_rout = Table.read('bootes_deep_lba_hbashift.Area_rms_rout2.fits')
trms_rin['rms'] = 1e-3*trms_rin['rms']  # put in Jy
trms_rout['rms'] = 1e-3*trms_rout['rms']  # put in Jy
trms_starmask['rms'] = 1e-3*trms_starmask['rms']  # put in Jy

#/net/lofar2//data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/bootes_deep_lba.compl_flux.fits
area = 23.6 #sq deg
area = 23.443604861110927  # max of trmstab
deeparea = 9.62237 # sq deg
deeparea = 8.609220312500138 # sq deg - sum of finite pixels in starmask

rad_circ = 1.
area_circ = np.pi*rad_circ*rad_circ
c0 = ac.SkyCoord(218.0,34.5,unit='deg')

#trms['Area']
#trms['rms']

# select only the 'proper' detected sources.... what to do about completeness of higher order wavelet scales??
print('Catalogue contains {0} sources'.format(len(cat)))
sel = (cat['Peak_flux'] / cat['Isl_rms']) > 5.
cat = cat[sel]
print('Catalogue contains {0} sources after selection'.format(len(cat)))

tflux = cat['Total_flux'] *1000.
pflux = cat['Peak_flux']*1000.
rms = cat['Isl_rms']*1000.




rms_area_interp = interp1d(trms['rms'],trms['Area_frac'],kind='linear', bounds_error=False, fill_value=(0.,trms['Area_frac'][-1]))
rms_area_interp_rout = interp1d(trms_rout['rms'],trms_rout['Area_frac'],kind='linear', bounds_error=False, fill_value=(0.,trms_rout['Area_frac'][-1]))
rms_area_interp_rin = interp1d(trms_rin['rms'],trms_rin['Area_frac'],kind='linear', bounds_error=False, fill_value=(0.,trms_rin['Area_frac'][-1]))

do_compl_cor = False
if compltable != 'None':
    do_compl_cor = True
    tcompl = Table.read(compltable)
    compl_area_interp = interp1d(tcompl['flux'],tcompl['detF'],kind='linear', bounds_error=False, fill_value=(0.,tcompl['detF'][-1]))
    compl_area_s_interp = interp1d(tcompl['flux'],tcompl['detF_err'],kind='linear', bounds_error=False, fill_value=(0.,tcompl['detF_err'][-1]))

    compl_area_circ_interp = interp1d(tcompl['flux'],tcompl['detFcent'],kind='linear', bounds_error=False, fill_value=(0.,tcompl['detFcent'][-1]))
    compl_area_s_circ_interp = interp1d(tcompl['flux'],tcompl['detFcent_err'],kind='linear', bounds_error=False, fill_value=(0.,tcompl['detFcent_err'][-1]))
    
smear_factor = 1.4
sigfactor = 5.
Avis = rms_area_interp(tflux / sigfactor)   
Avis = rms_area_interp(pflux * smear_factor / sigfactor)   

#Compl = compl_area_interp(tflux)

fmin = np.min(tflux)
fmax = np.max(tflux)

fbins = np.linspace(fmin, fmax, 20.)
fbins = np.logspace(np.log2(fmin), np.log2(fmax), 20., base=2)
fbinc = np.sqrt((fbins[1:]*fbins[:-1] ))
dfbinc = fbins[1:] - fbins[:-1]
nbins,_ = np.histogram(tflux,bins=fbins) 

sigval = 5. *1.4
## note in practice check also Stot/Speak...
rmsbins = fbinc / sigval
Abins = []
for rms in rmsbins:
    Ai = np.min((np.sum(trms['rms'] < rms), len(trms)-1))
    A = trms['Area_frac'][Ai]
    Abins.append(A)
Abins = np.array(Abins)

print (nbins, rmsbins, Abins)
nbins_vcor,_ = np.histogram(tflux,bins=fbins, weights=1./Avis) 
nbins_weights,_ = np.histogram(np.ones_like(tflux),bins=fbins, weights=1./Avis) 

f,ax = pp.paper_single_ax()
ax.plot(fbinc, nbins, label='raw')
ax.plot(fbinc, nbins_vcor, label='visibility area cor')
ax.loglog()
pp.set_attrib(ax, xlabel='$S_{\mathrm{total}}$ (mJy)', ylabel='$N(S_{\mathrm{total}})$')
pp.fig_save_many(f, plotpath+'raw_counts')



def get_deZotti_counts(f, alpha=-0.8, value='all'):
    dat = np.genfromtxt('/net/beerze/data2/wwilliams/projects/src_counts/deZotti2010_compilation.dat', skip_header=22, dtype=[('logS','<f8'),('cnt','<f8'),('u_e_cnt','<f8'),('l_e_cnt','<f8'),('var','<f8'),('ref','S5')])
    # col. 1 =      log(S [Jy])
# col. 2 =      S^2.5dN/dS[Jy^1.5/sr]
# col. 3-4 =    positive and negative errors on the counts.
# col. 5 =      For surveys covering less than 25 deg^2 the contribution to such errors
#               due to the sampling variance [eq.(6) in De Zotti et al. (2009)];
#               such contributions are negligible for larger areas.
# col. 6 =      references (see the list below)
#    bo08: 2008ApJ...681.1129B
#    br72: 1972AJ.....77..405B
#    ci99: 1999MNRAS.302..222C
#    fo06: 2006ApJS..167..103F
#    gr99: 1999MNRAS.305..297G
#    ho03: 2003AJ....125..465H
#    ke08: 2008ApJS..179...71K
#    ow08: 2008AJ....136.1889O
#    ri00: 2000ApJ...533..611R
#    se08A:2008MNRAS.386.1695S
#    wi97: 1997ApJ...475..479W
    
    if value != 'all':
        dat = dat[dat['ref']==value]
    
    bin_centres = 10**(dat['logS'])  # to Jy
    counts = dat['cnt']
    count_errors = np.array([dat['l_e_cnt'], dat['u_e_cnt']])
    
    #counts = counts*(1e3)**2.5


    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors 


def get_bonato_counts(f, alpha=-0.8, value='all'):
    ## S [mJy], Count [Jy^1.5/sr], Err superior Count [Jy^1.5/sr], Err inferior Count [Jy^1.5/sr], Reference
    dat = np.genfromtxt('/net/beerze/data2/wwilliams/projects/src_counts/Bonato_counts.txt', skip_header=22, dtype=[('S','<f8'),('cnt','<f8'),('u_e_cnt','<f8'),('l_e_cnt','<f8'),('ref','S5')])
    # col. 1 =      log(S [Jy])
# col. 2 =      S^2.5dN/dS[Jy^1.5/sr]
# col. 3-4 =    positive and negative errors on the counts.
# col. 5 =      For surveys covering less than 25 deg^2 the contribution to such errors
#               due to the sampling variance [eq.(6) in De Zotti et al. (2009)];
#               such contributions are negligible for larger areas.
# col. 6 =      references (see the list below)
#    bo08: 2008ApJ...681.1129B
#    br72: 1972AJ.....77..405B
#    ci99: 1999MNRAS.302..222C
#    fo06: 2006ApJS..167..103F
#    gr99: 1999MNRAS.305..297G
#    ho03: 2003AJ....125..465H
#    ke08: 2008ApJS..179...71K
#    ow08: 2008AJ....136.1889O
#    ri00: 2000ApJ...533..611R
#    se08A:2008MNRAS.386.1695S
#    wi97: 1997ApJ...475..479W
    
    if value != 'all':
        dat = dat[dat['ref']==value]
    
    bin_centres = dat['S']/1e3  # to Jy
    counts = dat['cnt']
    count_errors = np.array([dat['l_e_cnt'], dat['u_e_cnt']])
    
    #counts = counts*(1e3)**2.5


    frq = 1400.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors 


def get_deephba_counts(f, alpha=-0.8):
    dat = Table.read('/net/beerze/data2/wwilliams/projects/lba_bootes_deep/python/Boo_v10_final_d1_f1_7sd.counts', format='ascii')
    
    
    bin_centres = dat['x']/1e3 # to Jy
    counts = dat['yfid']
    count_errors = np.array([dat['zsysl'], dat['zsysu']])
    
    #counts = counts*(1e3)**2.5


    frq = 144.
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors 



def get_lolss_counts(f, alpha=-0.8):
    dat = np.load('/net/beerze/data2/wwilliams/projects/lofar/lolss/lolss_src_counts.npz')
    
    
    bin_centres = dat['binc'] # to Jy
    counts = dat['cnt']
    count_errors = np.array([dat['l_e_cnt'], dat['u_e_cnt']])
    
    #counts = counts*(1e3)**2.5


    frq = 54.393
    
    # scale to required freq
    bin_centres = bin_centres*(f/frq)**alpha
    counts = counts*(f/frq)**(1.5*alpha)
    count_errors = count_errors*(f/frq)**(1.5*alpha)
    
    return frq, bin_centres ,counts, count_errors 


def PE(cnt):
    
    # correct poisson errors
    PE = {"S": 1, "beta": 0., "gamma": 1.}
    
    
    # poisson errors ! low counts
    if cnt > 0:
        cnt_u_err = (cnt + PE['S']*np.sqrt(cnt+0.75) + (PE['S']**2+3)/4.) - cnt
        cnt_l_err = cnt - ( cnt*(1. - 1./(9.*cnt) - PE['S']/(3.*np.sqrt(cnt)) + PE['beta']*cnt**PE['gamma'] )**3. )
    else:
        cnt_u_err = 0
        cnt_l_err = 0
    return cnt_l_err, cnt_u_err

def get_source_counts(fluxes, area, compl_area_interp=None,  compl_area_s_interp=None, vis_area=None, scale_vis_rms=None, pfluxes=None, pfluxerrs=None, rmss=None, sigfactor=5., smear_factor=1., fbins=None, nbins=20, verbose=False, rescor=False, mode='each'):
    '''
    fluxes in Jy, rms in Jy/bm
    area in sqdeg
    
    visarea should be a table of rms vs Area_frac
    scale_vis_rms = factor by which we should scale the rms in the vis_area tbale to be in Jy/bm
    
    mode is avg for average per bin, each for vis area and peak flux per source
    '''
    # area -> sr
    degtosterad = (np.pi/180.)**2.
    area = area*degtosterad
        
    # sort by flux
    ind = fluxes.argsort()
    fluxes = fluxes[ind]
    if pfluxes is not None:
        pfluxes = pfluxes[ind]
    rmss = rmss[ind]
    
    Nsrcs = len(fluxes)
    
    
    if vis_area is not None:
        print ('using vis_area')
        assert isinstance(vis_area, Table), "vis_area should be an astropy.tables.Table instance"
        
        if scale_vis_rms is not None:
            vis_area['rms'] = vis_area['rms'] * scale_vis_rms
        
        #Avis = []
        ##for rmsi in rms:
            ##Ai = np.min((np.sum(trms['rms'] < rmsi), len(trms)-1))
        #maxind = len(trms)-1
        #for pflux in pfluxes:
            #rmsi = pflux / sigfactor
            #Ai = np.min((np.sum(vis_area['rms'] < rmsi), maxind))
            #A = vis_area['Area_frac'][Ai]
            #Avis.append(A)
        #Avis = np.array(Avis)
        
        def get_avis_err(rms_area_interp, pfluxes, pfluxerrs,smear_factor,sigfactor, Nsamp=100):
            
            Avis_arr = np.zeros((len(pfluxes),Nsamp))
            pfluxes_samp = np.zeros((len(pfluxes),Nsamp))
            for i in range(len(pfluxes)):
                pfluxes_samp[i,:] = np.random.normal(pfluxes[i], pfluxerrs[i], size=Nsamp)
            for i in range(Nsamp):
                Avis_arr[:,i] = rms_area_interp(pfluxes_samp[:,i] *smear_factor/ sigfactor) 
            
            Avis_err = np.std(Avis_arr, axis=1) 
            #print ('got Avis errors')
            #print (Avis_err.shape)
            return Avis_err
        
        
        rms_area_interp1 = interp1d(vis_area['rms'],vis_area['Area_frac'],kind='linear', bounds_error=False, fill_value=(0.,vis_area['Area_frac'][-1]))
        Avis = rms_area_interp1(pfluxes *smear_factor/ sigfactor)   
        Avis_err =  get_avis_err(rms_area_interp1,pfluxes, 0.5*pfluxerrs,smear_factor,sigfactor)
            
        #Compl = np.ones_like(fluxes) 
    else:
        print('no vis_area')
        Avis = np.ones(Nsrcs)
        Avis_err = np.zeros(Nsrcs)
    
    if compl_area_interp:
        Compl = compl_area_interp(fluxes)
        if compl_area_s_interp: 
            sCompl = compl_area_s_interp(fluxes) 
        else:
            sCompl = np.ones(Nsrcs)
    else:
        print('no completeness v flux')
        Compl = np.ones(Nsrcs)
        sCompl = np.ones(Nsrcs)
    
    
    # calculate the flux bins
    if fbins is None:
        
        
        # they're sorted - use the second min and fourth max
        fmin = fluxes[1]
        fmax = fluxes[-4]

        fbins = np.logspace(np.log2(fmin), np.log2(fmax), nbins+1, base=2)
        ## test divide first bin in 2
        #fbins=np.hstack((fbins[0], (fbins[0]+fbins[1])/2, fbins[1:]))
    else:
        assert isinstance(fbins, np.ndarray), 'fbins should be specified as a numpy array'
        # tbc should also be monotonically increasing
        assert np.all(np.diff(fbins) > 0), 'fbins should be monotonically increasing'
        nbins = len(fbins) - 1
    
    outtable = Table([Column(name='bin', data=np.zeros(nbins,dtype=int)),
                      Column(name='flux_mid', data=np.zeros(nbins,dtype=float)),
                      Column(name='flux_low', data=np.zeros(nbins,dtype=float)),
                      Column(name='flux_high', data=np.zeros(nbins,dtype=float)),
                      Column(name='flux_mean', data=np.zeros(nbins,dtype=float)),
                      Column(name='dn', data=np.zeros(nbins,dtype=int)),
                      Column(name='dn_l', data=np.zeros(nbins,dtype=float)),
                      Column(name='dn_u', data=np.zeros(nbins,dtype=float)),
                      Column(name='area_mean', data=np.zeros(nbins,dtype=float)),
                      Column(name='compl_mean', data=np.zeros(nbins,dtype=float)),
                      Column(name='compl_err_mean', data=np.zeros(nbins,dtype=float)),
                      Column(name='rescor', data=np.zeros(nbins,dtype=float)),
                      Column(name='rescor_err', data=np.zeros(nbins,dtype=float)),
                      Column(name='dnds', data=np.zeros(nbins,dtype=float)),
                      Column(name='E_src_cnt_raw', data=np.zeros(nbins,dtype=float)),
                      Column(name='E_src_cnt_raw_l', data=np.zeros(nbins,dtype=float)),
                      Column(name='E_src_cnt_raw_u', data=np.zeros(nbins,dtype=float)),
                      Column(name='E_src_cnt', data=np.zeros(nbins,dtype=float)),
                      Column(name='E_src_cnt_l', data=np.zeros(nbins,dtype=float)),
                      Column(name='E_src_cnt_u', data=np.zeros(nbins,dtype=float)),])
    
    # goemetric mean
    fbinc = np.sqrt((fbins[1:]*fbins[:-1] ))
    dfbinc = fbins[1:] - fbins[:-1]
    
    dnds = np.zeros_like(fbins)
    
    if verbose:
        print("flow-fhigh area compl rescor count  Escnt")
    for i in range(nbins):
        rescori = 1.
        
        selind = (fluxes >= fbins[i]) & (fluxes < fbins[i+1])
        
        fi = fbinc[i]
        df = dfbinc[i]
        
        
        if rescor:
            
            COR_RES_arr = get_res_cor(fi)
            
            COR_RES = COR_RES_arr[2]
            minCOR_RES = COR_RES_arr.min()
            maxCOR_RES = COR_RES_arr.max()
            
            counts_corres = COR_RES
            counts_corres_err = np.sqrt((0.1*COR_RES)**2.+ (maxCOR_RES-minCOR_RES)**2.)
            
        else:
            counts_corres = 1.
            counts_corres_err = 0.
        
        fmean = np.average(fluxes[selind])
        
        dn = np.sum(selind)
        #dncor = np.sum(1./Avis[selind])
        #dncor = np.sum(1./Compl[selind])
        dncor = counts_corres *np.sum(1./(Compl[selind]*Avis[selind]))
        
        
        
        Avis_avg = np.average(Avis[selind])
        Avis_err_avg = np.average(Avis_err[selind])
        Compl_avg = np.average(Compl[selind])
        sCompl_avg = np.average(sCompl[selind])
        if verbose: print('**',Avis[selind], Avis_err[selind])
        if verbose: print(Avis_avg, Avis_err_avg)
        Compl_avg = np.average(Compl[selind])
        
        dn_l, dn_u = PE(dn)
        #dn_cor_l, dn_cor_u = PE(dncor)
        dn_cor_l1, dn_cor_u1 = PE(dn)
        
        #dn_cor_l = np.sqrt((dncor**2)*(np.sum((Avis_err[selind]/Avis[selind])**2)) + np.sum((sCompl[selind]/Compl[selind])**2)) +  counts_corres_err**2.)
        #dn_cor_u = np.sqrt( (dncor**2)*(np.sum((Avis_err[selind]/Avis[selind])**2)) + np.sum((sCompl[selind]/Compl[selind])**2)) +  counts_corres_err**2.)
        dn_cor_l = np.sqrt((dncor**2)*((Avis_err_avg/Avis_avg)**2) + (sCompl_avg/Compl_avg)**2 +  counts_corres_err**2.)
        dn_cor_u = np.sqrt((dncor**2)*((Avis_err_avg/Avis_avg)**2) + (sCompl_avg/Compl_avg)**2 +  counts_corres_err**2. )
    
        if mode =='avg':
            dncor = counts_corres *dn / Compl_avg
            Nsamp = 1000
            dncor_arr = np.random.normal(loc=counts_corres,scale=counts_corres_err,size=Nsamp) * np.random.poisson(lam=dn,size=Nsamp) / np.random.normal(loc=Compl_avg,scale=sCompl_avg,size=Nsamp)
            dn_cor_l,dn_cor_u =  np.percentile(dncor_arr, [16,84])
            dn_cor_u = dn_cor_u - dncor
            dn_cor_l = dncor - dn_cor_l
            #dn_cor_u = (dn+dn_u)/(Compl_avg-sCompl_avg)   - dncor
            #dn_cor_l =  dncor - (dn-dn_l)/(Compl_avg-sCompl_avg) 
            #dn_cor_u = np.sqrt(dn_cor_u**2 + counts_corres_err**2.)
            #dn_cor_l = np.sqrt(dn_cor_l**2 + counts_corres_err**2.)
            
            #dn_cor_l = np.sqrt(dn_l**2 + (sCompl_avg/Compl_avg)**2 + counts_corres_err**2.)
            #dn_cor_u = np.sqrt(dn_u**2 + (sCompl_avg/Compl_avg)**2 + counts_corres_err**2.)
        
        dnds_cor =  dncor / df
        Escnt = dnds_cor * (fi**2.5) / area
        Escnt_l = (dn_cor_l/df)* (fi**2.5) / area
        Escnt_u = (dn_cor_u/df)* (fi**2.5) / area
        
        
        dnds = dn / df
        Escnt_raw = dnds * (fi**2.5) / area
        Escnt_raw_l = (dn_l/df)* (fi**2.5) / area
        Escnt_raw_u = (dn_u/df)* (fi**2.5) / area
        
        
        if verbose:
            print("{f1:5.2f}-{f2:5.2f} {a:5.2f} {c:5.2f} {r:5.2f}+/-{re:5.2f} {dn:5d}+{dnu:3.0f}-{dnl:3.0f}  {dnc:5.0f}+{dncu:3.0f}-{dncl:3.0f} {ec:5.2f}+{ecu:5.2f}-{ecl:5.2f}".format(f1=fbins[i], f2=fbins[i+1], a=Avis_avg, c=Compl_avg, r=counts_corres,re=counts_corres_err, dn=dn,dnl=dn_l,dnu=dn_u, dnc=dncor, dncl=dn_cor_l, dncu=dn_cor_u, dnds=dnds, ec=Escnt, ecl=Escnt_l, ecu=Escnt_u))
            
        outtable['bin'][i] = i
        outtable['flux_low'][i] = fbins[i]
        outtable['flux_high'][i] = fbins[i+1]
        outtable['flux_mid'][i] = fi
        outtable['flux_mean'][i] = fmean
        outtable['dn'][i] = dn
        outtable['dn_l'][i] = dn_l
        outtable['dn_u'][i] = dn_u
        outtable['area_mean'][i] = Avis_avg
        outtable['compl_mean'][i] = Compl_avg
        outtable['compl_err_mean'][i] = sCompl_avg
        outtable['rescor'][i] = counts_corres
        outtable['rescor_err'][i] = counts_corres_err
        outtable['dnds'][i] = dnds
        outtable['E_src_cnt_raw'][i] = Escnt_raw
        outtable['E_src_cnt_raw_l'][i] = Escnt_raw_l
        outtable['E_src_cnt_raw_u'][i] = Escnt_raw_u
        outtable['E_src_cnt'][i] = Escnt
        outtable['E_src_cnt_l'][i] = Escnt_l
        outtable['E_src_cnt_u'][i] = Escnt_u
        
        for col in outtable.colnames:
            if outtable[col].info.dtype is np.dtype('float64'):
                outtable[col].format = '%.2f'
        
            
    return outtable


c = ac.SkyCoord(cat['RA'],cat['DEC'],unit='deg')
c0 = ac.SkyCoord(218.0,34.5,unit='deg')
seps = c0.separation(c).to(u.deg).value
smear_factors = np.ones(len(cat))
smear_factors[seps<3.] = 2.1
smear_factors[seps<2.5] = 1.9
smear_factors[seps<2.] = 1.7
smear_factors[seps<1.5] = 1.6
smear_factors[seps<1.0] = 1.5
smear_factors[seps<0.5] = 1.4
smear_factors = 1./smear_factors
    

# do rescor 
#srctabr = get_source_counts(cat['Total_flux'], area, vis_area=trms, scale_vis_rms=1, pfluxes=cat['Peak_flux'], pfluxerrs=cat['E_Peak_flux'], rmss=cat['Isl_rms'], sigfactor=5., smear_factor=1., fbins=None, nbins=25, rescor=deephba, verbose=True)

srctab = get_source_counts(cat['Total_flux'], area, vis_area=trms, scale_vis_rms=1., pfluxes=cat['Peak_flux'], pfluxerrs=cat['E_Peak_flux'], rmss=cat['Isl_rms'], sigfactor=5., smear_factor=1, fbins=None, nbins=24, verbose=False, rescor=True)
mcirc = (c0.separation(c) < rad_circ*u.deg)
mcirc2 = (c0.separation(c) >= rad_circ*u.deg)
#srctab_circ = get_source_counts(cat['Total_flux'][mcirc], area_circ, vis_area=trms, scale_vis_rms=1, pfluxes=cat['Peak_flux'][mcirc], pfluxerrs=cat['E_Peak_flux'][mcirc], rmss=cat['Isl_rms'][mcirc], sigfactor=5., smear_factor=1., fbins=None, nbins=25, verbose=False, rescor=True)
srctab_circ = get_source_counts(cat['Total_flux'][mcirc], area_circ, vis_area=trms_rin, scale_vis_rms=1, pfluxes=cat['Peak_flux'][mcirc], pfluxerrs=cat['E_Peak_flux'][mcirc], rmss=cat['Isl_rms'][mcirc], sigfactor=5., smear_factor=1., fbins=srctab['flux_low'], verbose=False, rescor=True, nbins=20)
srctab_circ_comp = get_source_counts(cat['Total_flux'][mcirc], area_circ, scale_vis_rms=1, rmss=cat['Isl_rms'][mcirc], sigfactor=5., smear_factor=1., fbins=srctab['flux_low'], verbose=False, rescor=True, nbins=20, compl_area_interp=compl_area_circ_interp, compl_area_s_interp=compl_area_s_circ_interp,mode='avg')
srctab_circ2 = get_source_counts(cat['Total_flux'][mcirc2], area-area_circ, vis_area=trms_rout, scale_vis_rms=1, pfluxes=cat['Peak_flux'][mcirc2], pfluxerrs=cat['E_Peak_flux'][mcirc2], rmss=cat['Isl_rms'][mcirc2], sigfactor=5., smear_factor=1, fbins=srctab['flux_low'], verbose=False, rescor=True)
mopt = cat['opt_flag']
srctab_opt = get_source_counts(cat['Total_flux'][mopt], deeparea, vis_area=trms_starmask, scale_vis_rms=1, pfluxes=cat['Peak_flux'][mopt], pfluxerrs=cat['E_Peak_flux'][mopt], rmss=cat['Isl_rms'][mopt], sigfactor=5., smear_factor=1., fbins=None, nbins=15, verbose=False, rescor=True)
magn = cat['flag_AGN']
srctab_agn = get_source_counts(cat['Total_flux'][magn], deeparea, vis_area=trms_starmask, scale_vis_rms=1, pfluxes=cat['Peak_flux'][magn], pfluxerrs=cat['E_Peak_flux'][magn], rmss=cat['Isl_rms'][magn], sigfactor=5., smear_factor=1., fbins=None, nbins=15, verbose=False, rescor=True)
msf = cat['flag_SFG']
srctab_sf = get_source_counts(cat['Total_flux'][msf], deeparea, vis_area=trms_starmask, scale_vis_rms=1, pfluxes=cat['Peak_flux'][msf], pfluxerrs=cat['E_Peak_flux'][msf], rmss=cat['Isl_rms'][msf], sigfactor=5., smear_factor=1., fbins=None, nbins=7, verbose=False, rescor=True)

if do_compl_cor:
    srctabc = get_source_counts(cat['Total_flux'], area, scale_vis_rms=1, rmss=cat['Isl_rms'], sigfactor=5., smear_factor=1., fbins=None, nbins=25, verbose=False, compl_area_interp=compl_area_interp, compl_area_s_interp=compl_area_s_interp, mode='avg', rescor=True)

# use smear_factor
srctabs = get_source_counts(cat['Total_flux'], area, vis_area=trms, scale_vis_rms=1, pfluxes=cat['Peak_flux'], pfluxerrs=cat['E_Peak_flux'], rmss=cat['Isl_rms'], sigfactor=5., smear_factor=smear_factors, fbins=None, nbins=25, verbose=False)



deepsrctab = get_source_counts(deephba['Total_flux'], deeparea, vis_area=None, scale_vis_rms=1, pfluxes=deephba['Peak_flux'], rmss=deephba['Isl_rms'], sigfactor=5., smear_factor=1.0, fbins=None, nbins=30, verbose=False)

Hfrq, Hbinc , Hcnts, Hcnt_errs = get_deephba_counts(54.393, alpha=-0.8)
Lfrq, Lbinc , Lcnts, Lcnt_errs = get_lolss_counts(54.393, alpha=-0.8)
Zfrq, Zbinc , Zcnts, Zcnt_errs = get_deZotti_counts(54.393, alpha=-0.8)
Bfrq, Bbinc , Bcnts, Bcnt_errs = get_bonato_counts(54.393, alpha=-0.8)
Bfrq1, Bbinc1 , Bcnts1, Bcnt_errs1 = get_bonato_counts(54.393, alpha=-0.6)
Zfrq1, Zbinc1 , Zcnts1, Zcnt_errs1 = get_deZotti_counts(54.393, alpha=-0.6)

f,ax = pp.paper_single_ax()
for fi in range(len(srctab)):
    fe1 = srctab['flux_low'][fi]
    fe2 = srctab['flux_high'][fi]
    selind = (deephba['Total_flux'] > fe1 ) & (deephba['Total_flux'] <= fe2)
    sizes = deephba['Size'][selind]
    h,b = np.histogram(sizes, bins=20, range=(0, 50))
    h0 = np.sum(sizes < 15.) / len(sizes)
    A = h[np.newaxis,:]
    c = ax.imshow(np.log10(A.T),extent=(1e3*fe1,1e3*fe2,0,50),origin='lower',vmin=0,vmax=10,aspect='auto')
    #ax.step(b[:-1],h,label='bin'+str(fi), where='pre')
    print ("{f0:.2f} {f1:.2f} {n:d} {f:.3f}".format(f0=1e3*fe1, f1=1e3*fe2,n=len(sizes), f=np.sum(sizes >15.)/len(sizes)))
ax.hlines(20,1e-2, 1e4)
cbar = plt.colorbar(c)
cbar.set_label('$\log N$')
ax.semilogx()
pp.set_attrib(ax, xlim=(1e3*deepsrctab['flux_low'][0],1e3*deepsrctab['flux_high'][-1]), xlabel='flux (mJy)',ylabel='size (arcsec)')
pp.fig_save_many(f, plotpath+'size_dist_by_flux')


f,ax = pp.paper_single_ax()
ax.errorbar(Bbinc*1e3, Bcnts, [Bcnt_errs[0], Bcnt_errs[1]], marker='.', color='cyan', linestyle='none', label=r'Bonato (scaled $\alpha=-0.8$)')
ax.loglog()
#pp.format_log_axis10(ax, axis='x')
pp.set_attrib(ax, xlabel='$S_{54\mathrm{\,MHz}}$ (mJy)', ylabel='$S^{2.5} dN/dS$ (Jy$^{1.5}$ sr$^{-1}$)', xlim=(6e-1,1.2e4), ylim=(2e0,7e4))
ax.legend()
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
pp.fig_save_many(f, plotpath+'esource_counts_bonato')

f,ax = pp.paper_single_ax()
lL1 = ax.errorbar(Zbinc*1e3, Zcnts, [Zcnt_errs[0], Zcnt_errs[1]], marker='.', color='gray', linestyle='none', label=r'1.4\,GHz counts (scaled $\alpha=-0.8$)',zorder=-10)
lL2 = ax.errorbar(Zbinc1*1e3, Zcnts1, [Zcnt_errs1[0], Zcnt_errs1[1]], marker='.', color='k', linestyle='none', label=r'1.4\,GHz counts (scaled $\alpha=-0.6$)',zorder=-10)
ax.errorbar(Bbinc*1e3, Bcnts, [Bcnt_errs[0], Bcnt_errs[1]], marker='.', color='gray', linestyle='none')
ax.errorbar(Bbinc1*1e3, Bcnts1, [Bcnt_errs1[0], Bcnt_errs1[1]], marker='.', color='k', linestyle='none')
lLoLSS = ax.errorbar(Lbinc*1e3, Lcnts, [Lcnt_errs[0], Lcnt_errs[1]], marker='.', color='C2', linestyle='none', label=r'LoLSS',zorder=-10)
#ax.errorbar(Hbinc*1e3, Hcnts, [Hcnt_errs[0], Hcnt_errs[1]], marker='.', color='C3', linestyle='none', label=r'Deep HBA')

ax.loglog()
#pp.format_log_axis10(ax, axis='x')
pp.set_attrib(ax, xlabel='$S_{54\mathrm{\,MHz}}$ (mJy)', ylabel='$S^{2.5} dN/dS$ (Jy$^{1.5}$ sr$^{-1}$)', xlim=(6e-1,1.2e4), ylim=(2e0,7e4))
leg1 = ax.legend(handles=[lL1, lL2, lLoLSS],loc='upper left')
ax.add_artist(leg1)
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())

S = srctab
lSr = ax.errorbar(S['flux_mid']*1e3, S['E_src_cnt_raw'], [S['E_src_cnt_raw_l'], S['E_src_cnt_raw_u']], marker='o', color='C0', mfc='none',linestyle='none', label='Raw')
#lSv = ax.errorbar(S['flux_mid']*1e3, S['E_src_cnt'], [S['E_src_cnt_l'], S['E_src_cnt_u']], marker='x', color='C0',linestyle='none', label='Corrected (vis area)')
pp.fig_save_many(f, plotpath+'esource_counts')
if do_compl_cor:
    S = srctabc
    lSc = ax.errorbar(S['flux_mid']*1e3, S['E_src_cnt'], [S['E_src_cnt_l'], S['E_src_cnt_u']], marker='o', color='C0',linestyle='none', label='Corrected')
    
    #leg2 = ax.legend(handles=[lSr,lSv,lSc],loc='lower right')
    leg2 = ax.legend(handles=[lSr,lSc],loc='lower right')
    pp.fig_save_many(f, plotpath+'esource_counts_withcompl')

#unsing inner circle only
plot_inner = True
if plot_inner:
    S = srctab_circ
    l1 = ax.errorbar(S['flux_mid']*1e3, S['E_src_cnt_raw'], [S['E_src_cnt_raw_l'], S['E_src_cnt_raw_u']], marker='o', color='C4', mfc='none',linestyle='none', label='Raw')
    #l1v = ax.errorbar(S['flux_mid']*1e3, S['E_src_cnt'], [S['E_src_cnt_l'], S['E_src_cnt_u']], marker='x', color='C4',linestyle='none', label='Corrected (vis area)')
    S = srctab_circ_comp
    l1c = ax.errorbar(S['flux_mid']*1e3, S['E_src_cnt'], [S['E_src_cnt_l'], S['E_src_cnt_u']], marker='o', color='C4',linestyle='none', label='Corrected')
    #leg2 = ax.legend(handles=[lSr,lSv,lSc,l1,l1v,l1c],loc='lower right')
    leg2 = ax.legend(handles=[lSr,lSc,l1,l1c],loc='lower right')
    pp.fig_save_many(f, plotpath+'esource_counts_withinner')

## using outer circle
#S = srctab_circ2
#ax.errorbar(S['flux_mid']*1e3, S['E_src_cnt_raw'], [S['E_src_cnt_raw_l'], S['E_src_cnt_raw_u']], marker='o', color='C3', mfc='none',linestyle='none', label='Raw')
#ax.errorbar(S['flux_mid']*1e3, S['E_src_cnt'], [S['E_src_cnt_l'], S['E_src_cnt_u']], marker='o', color='C3',linestyle='none', label='Corrected (vis area)')

plot_source_types = False
if plot_source_types:
    S = srctab_sf
    lsf =ax.errorbar(S['flux_mid']*1e3, S['E_src_cnt'], [S['E_src_cnt_l'], S['E_src_cnt_u']], marker='o', color='C1',linestyle='none', label='SF')
    S = srctab_agn
    lagn = ax.errorbar(S['flux_mid']*1e3, S['E_src_cnt'], [S['E_src_cnt_l'], S['E_src_cnt_u']], marker='o', color='C2',linestyle='none', label='AGN')
    #ax.errorbar(srctab_opt['flux_mid']*1e3, srctab_opt['E_src_cnt'], [srctab_opt['E_src_cnt_l'], srctab_opt['E_src_cnt_u']], marker='o', color='C3',linestyle='none', label='opt flag')

    leg2 = ax.legend(handles=[lSr,lSv,lSc,l1,l1v,l1c,lsf,lagn],loc='lower right')
    pp.fig_save_many(f, plotpath+'esource_counts_types')
#ax.errorbar(deepsrctab['flux_mid']*1e3, deepsrctab['E_src_cnt'], [deepsrctab['E_src_cnt_l'], deepsrctab['E_src_cnt_u']], marker='.', color='C1', label='DeepHBA')
#ax.errorbar(srctabs['flux_mid']*1e3, srctabs['E_src_cnt'], [srctabs['E_src_cnt_l'], srctabs['E_src_cnt_u']], marker='.', color='C2', label='vis area smeared')


tf = np.logspace(np.log10(cat['Total_flux'].min()), np.log10(cat['Total_flux'].max()), 1000)
pf = np.logspace(np.log10(cat['Peak_flux'].min()), np.log10(cat['Peak_flux'].max()), 1000)

tf = np.logspace(np.log10(cat['Total_flux'].min()), 1, 1000)
pf = np.logspace(np.log10(cat['Peak_flux'].min()), 1, 1000)

f,ax = pp.paper_single_ax()
ax.semilogx()
#pp.fig_save_many(f, 'plots/compl_flux')
ax.plot(pf*1e3, rms_area_interp(pf/5.), c='C0', label='RMS area')
ax.plot(tf*1e3, compl_area_interp(tf), c='C1', label='Detected Fraction')
ax.semilogx()
ax.legend()
#pp.set_attrib(ax, xlabel='RMS (mJy)', ylabel='Area Fraction')
pp.set_attrib(ax, xlabel='$S_{54\mathrm{\,MHz}}$ (mJy)', ylabel='Fraction',xlim=(4,100))
pp.fig_save_many(f, plotpath+'compl_flux')

rms = np.logspace(np.log10(1), np.log10(20), 100)*1e-3
f,ax = pp.paper_single_ax()
ax.plot(rms*1e3, rms_area_interp(rms),c='C1',label='interp')
ax.plot(trms['rms']*1e3, trms['Area_frac'],'.',c='C0',label='raw')
pp.set_attrib(ax, xlabel='RMS (mJy)', ylabel='Area Fraction',xlim=(1,30))
ax.semilogx()
ax.legend()
pp.fig_save_many(f, plotpath+'vis_area')


## in flux bins compare the LBA counts with the de Zotti counts to find what spec ind fits best...
fluxes = cat['Total_flux'].copy()
ind = fluxes.argsort()
fluxes = fluxes[ind]
fmin = fluxes[1]
fmax = fluxes[-4]
nbins = 5

Zfrq, Zbinc , Zcnts, Zcnt_errs = get_deZotti_counts(54.393, alpha=-0.8)
alpha = np.arange(-0.85, -0.61, 0.01)
Zbinc = np.nan*np.ones((len(Zbinc), len(alpha)))
Zcnts = np.nan*np.ones((len(Zbinc), len(alpha)))
Zcnt_errs = np.nan*np.ones((2, len(Zbinc), len(alpha)))
for ai in range(len(alpha)):
    _, Zbinc[:,ai] , Zcnts[:,ai], Zcnt_errs[:,:,ai] = get_deZotti_counts(54.393, alpha=alpha[ai])


SS = srctabc.copy()
#Zcnts = Zcnts / (Zbinc**2.5)
#SS['E_src_cnt'] = SS['E_src_cnt'] / (SS['flux_mid']**2.5)

fbins = np.logspace(np.log2(fmin), np.log2(fmax), nbins+1, base=2)
fbinc = np.sqrt((fbins[1:]*fbins[:-1] ))
dfbinc0 = fbins[1:] - fbinc
dfbinc1 = fbinc - fbins[:-1]
min_alpha = np.nan*np.ones_like(fbinc)
for fi in range(nbins):
    si = (SS['flux_mid']>=fbins[fi]) &(SS['flux_mid']<fbins[fi+1])
    
    fmin = np.min(SS['flux_low'][si])
    fmax = np.max(SS['flux_high'][si])
    print(fi, np.sum(si), fbinc[fi], fmin, fmax)
    
    sep =  np.nan*np.ones(len(alpha))
    for ai in range(len(alpha)):
        zi = (Zbinc[:,ai]>=fmin)&(Zbinc[:,ai]<fmax)
        print(ai, np.sum(zi))
        sep[ai] = (np.average(Zcnts[zi,ai],weights=1./(Zcnt_errs[0,zi,ai]**2.+Zcnt_errs[1,zi,ai]**2.)) - np.average(SS['E_src_cnt'][si],weights=1./(SS['E_src_cnt_l'][si]**2+SS['E_src_cnt_u'][si]**2.)))**2.
        sep[ai] = (np.average(Zcnts[zi,ai]) - np.average(SS['E_src_cnt'][si]))**2.
    print(sep, alpha[sep.argmin()])
    
    min_alpha[fi] = alpha[sep.argmin()]



nalp = np.load(plotpath+'deephba_hbashift_flux_ratio.npz')

f,ax = pp.paper_single_ax()
#ax.plot(fbinc*1e3, min_alpha)
ax.errorbar(fbinc*1e3, min_alpha, yerr=None, xerr=[dfbinc1*1e3,dfbinc0*1e3],marker='o',c='C0',linestyle='none', label='$\\nu=1400\mathrm{\,MHz}$')
ax.errorbar(nalp['arr_0'],nalp['arr_1'],nalp['arr_2'],nalp['arr_3'],marker='o',c='C2',linestyle='none', label='$\\nu=144\mathrm{\,MHz}$')
#ax.errorbar(nalp['X'],nalp['Y'],nalp['Yerr'],nalp['Xerr'],marker='o',c='C2',linestyle='none', label='$\\nu=144\mathrm{\,MHz}$')
ax.semilogx()
ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2) )
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.05) )
pp.set_attrib(ax, xlabel='$S_{54\mathrm{\,MHz}}$ (mJy)', ylabel='$\\alpha_{54}^{\\nu}$',xlim=(6e-1,1.2e4),ylim=(-1., -0.4))
ax.legend()
pp.fig_save_many(f, plotpath+'src_cnts_alpha_scaling')




## at 130-200 mJy there is a shift - from being consistent with -0.8 to -0.6 ?? real ?? possibly at 30mJy at 150MHz
## below 13 mJy incompleteness??

