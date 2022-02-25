import glob
import os
import sys
import numpy as np
import astropy.units as u
from astropy.table import Table, Column, vstack
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.io.fits as fits
import mocpy
from utils.make_subim import flatten    
from process_cat import Make_Shape, sourcename


clobber = True



import configparser
config = configparser.ConfigParser()
config.read('image.config')
plotpath = config['SETUP']['plotpath']
dpath = config['SETUP']['dpath']
deeptable = config['deeptable']['deeptable']
deepids = config['deeptable']['deepids']



cat =  config['fluxes']['cattable']
gcat =  config['fluxes']['gcattable']
tcat = Table.read(dpath+cat)
tcompcat = Table.read(dpath+cat)
tgcat = Table.read(dpath+gcat)

catout =  config['finalcat']['catout']
compcatout =  config['finalcat']['compcat']


deepcat = config['deeptable']['deepcat']
deepcato = config['deeptable']['deeptable']
starmask = config['deeptable']['starmask']

#omoc = mocpy.MOC.from_fits('/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes_merged_optical/final_mocs/star_mask/Bootes_smaskb_asec_MOC.fits')
#smoc = mocpy.MOC.from_fits('/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes_merged_optical/final_mocs/star_mask/Bootes_smaskb_chi2s_MOC.fits')
#optmoc = mocpy.MOC.from_fits('/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes_merged_optical/final_mocs/Bootes_i_MOC.fits')
#starmask = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/data_release/bootes//image_full_ampphase_di_m.NS_shift.blanked.scaled.rms_spmask.fits'



tdeepcat = Table.read(deepcat)
tdeepcato = Table.read(deepcato)



def filter_table(t,mask):

# arguments: input cat, mask fits, 

    mask=fits.open(mask)
    if 'CTYPE3' in mask[0].header:
        print('Flattening the mask')
        mask=flatten(mask,None,None,None,None,None)
    
    w=WCS(mask[0].header)

    pos=w.wcs_world2pix(t['RA'],t['DEC'],0)
    maxy,maxx=mask[0].data.shape
    filter=[]
    for i,r in enumerate(t):
        if (i % 5000)==0:
            print(i)
        x=int(pos[0][i])
        y=int(pos[1][i])
        if x<0 or y<0 or x>=maxx or y>=maxy:
            inmask=False
        else:
            inmask=~np.isnan(mask[0].data[y,x])

        filter.append(inmask)
    return np.array(filter)


Ncat = len(tcat)

# change to arcsec
for col in  ['Maj','Min','DC_Maj','DC_Min']:
    tcat[col] = 3600.*tcat[col]
    tgcat[col] = 3600.*tgcat[col]
    tcompcat[col] = 3600.*tcompcat[col]

if 'mcode' not in tcat.colnames:
    tcat.add_column(Column(name='mcode',data=-99*np.ones(Ncat, dtype=int)))
if 'Artefact' not in tcat.colnames:
    tcat.add_column(Column(name='Artefact',data=np.zeros(Ncat, dtype=bool)))
if 'handled' not in tcat.colnames:
    tcat.add_column(Column(name='handled',data=np.zeros(Ncat, dtype=bool)))
if 'confused' not in tcat.colnames:
    tcat.add_column(Column(name='confused',data=np.zeros(Ncat, dtype=bool)))

if 'HBA_Deep_flag' not in tcat.colnames:
    tcat.add_column(Column(name='HBA_Deep_flag',data=np.ones(Ncat, dtype=bool)))
if 'opt_flag' not in tcat.colnames:
    tcat.add_column(Column(name='opt_flag',data=np.zeros(Ncat, dtype=bool)))
#if 'opt_star_flag' not in tcat.colnames:
    #tcat.add_column(Column(name='opt_star_flag',data=np.zeros(Ncat, dtype=bool)))
#if 'ir_star_flag' not in tcat.colnames:
    #tcat.add_column(Column(name='ir_star_flag',data=np.zeros(Ncat, dtype=bool)))
    
if 'Deblended_From' not in tcat.colnames:
    tcat.add_column(Column(name='Deblended_From',data=np.zeros(Ncat, dtype='S25')))
if 'HBA_Deep_Name' not in tcat.colnames:
    tcat.add_column(Column(name='HBA_Deep_Name',data=np.zeros(Ncat, dtype='S22')))
    tcat.add_column(Column(name='HBA_Deep_Sep',data=np.nan*np.zeros(Ncat, dtype=float)))
    tcat.add_column(Column(name='HBA_Deep_ALPHA_J2000',data=np.nan*np.zeros(Ncat, dtype=float)))
    tcat.add_column(Column(name='HBA_Deep_DELTA_J2000',data=np.nan*np.zeros(Ncat, dtype=float)))
if 'HBA_DeepWide_id' not in tcat.colnames:
    tcat.add_column(Column(name='HBA_DeepWide_id',data=-1*np.ones(Ncat, dtype=int)))
    tcat.add_column(Column(name='HBA_DeepWide_Sep',data=np.nan*np.zeros(Ncat, dtype=float)))

if 'Size' not in tcat.colnames:
    tcat.add_column(Column(name='Size',data=np.nan*np.zeros(Ncat, dtype=float)))
    tcat['Size'] = 2*tcat['DC_Maj']  # in arcsec
    tcat.add_column(Column(name='Width',data=np.nan*np.zeros(Ncat, dtype=float)))
    tcat['Width'] = 2*tcat['DC_Min']  # in arcsec
    tcat.add_column(Column(name='Angle',data=np.nan*np.zeros(Ncat, dtype=float)))
    tcat['Angle'] = tcat['DC_PA']
     
     
if 'Component_Name' not in tcompcat.colnames:
    tcompcat.add_column(Column(name='Component_Name',data=tcompcat['Source_Name']))
     
     
fcodelist = 'lba_match_deephba_codes.txt'
if not os.path.isfile(fcodelist):
    print('need match codes: run match_to_deep_hba.py')
    
tmcode = Table.read(fcodelist,format='ascii')
#print('(0) no match')
#print('(1) simple nearest match ok - default')
#print('(2) complex')
#print('(3) lba source is artefact')
#print('(4) no match possible')
#print('(5) multiple match, no gaus components')
#print('(6) special case HBA is wrong and misses source assoc and id')
#print('(7) special case part of LBA/HBA is missing/ needs by hand id')

for i in range(len(tcat)):
    mcode = tmcode[tmcode['Source_Name'] == tcat[i]['Source_Name']]['match_code'][0]
    tcat[i]['mcode'] = mcode


# split into gauss components first, then merge, because a few new sources are made by remerging the gaus components
    

#colnames = ['Source_id',  'Isl_id',  'RA',  'E_RA',  'DEC',  'E_DEC',  'Total_flux',  'E_Total_flux',  'Peak_flux',  'E_Peak_flux',  ,  'Maj',  'E_Maj',  'Min',  'E_Min',  'PA',  'E_PA',  'Maj_img_plane',  'E_Maj_img_plane',  'Min_img_plane',  'E_Min_img_plane',  'PA_img_plane',  'E_PA_img_plane',  'DC_Maj',  'E_DC_Maj',  'DC_Min',  'E_DC_Min',  'DC_PA',  'E_DC_PA',  'DC_Maj_img_plane',  'E_DC_Maj_img_plane',  'DC_Min_img_plane',  'E_DC_Min_img_plane',  'DC_PA_img_plane',  'E_DC_PA_img_plane',  'Isl_Total_flux',  'E_Isl_Total_flux',  'Isl_rms',  'Isl_mean',  'Resid_Isl_rms',  'Resid_Isl_mean',  'S_Code',  'Source_Name',  'Ngaus',  'mcode',  'Artefact',  'handled',  'HBA_Deep_flag',  'opt_flag',  'opt_star_flag',  'ir_star_flag',  'HBA_Deep_Name',  'HBA_Deep_Sep',  'Size',  'Width',  'Angle']

# remove these cols from the source catalogue - we don't have them from the gaus cat and don't use them
remove_cols = ['RA_max',  'E_RA_max',  'DEC_max',  'E_DEC_max']
tcat.remove_columns(remove_cols)
tcompcat.remove_columns(remove_cols)

tgcat.add_column(Column(name='Component_Name',data=np.zeros(len(tgcat), dtype='S25')))
for ig in range(len(tgcat)):
    tgcat['Component_Name'][ig] = sourcename(tgcat['RA'][ig],tgcat['DEC'][ig])

catcols = tcat.colnames
gcatcols = tgcat.colnames

print('cat now has ',len(tcat),'entries')
# get the sources to split into component gaussians
gauslist = []
newgaussrcs = []
glist = glob.glob('LBABOO*gaussplit')
for g in glist:
    gauslist.append(g.replace('.gaussplit',''))
print('To split into gaussian components')
print(gauslist)
for g in gauslist:
    print(g)
    catind = np.where(tcat['Source_Name'] == g)[0][0]
    catkeep = np.where(tcat['Source_Name'] != g)[0]
    gcatind = np.where(tgcat['Source_id'] == tcat['Source_id'][catind])[0]
    print(catind,'->',gcatind)
    
    tnewi = Table(tcat[catind])
    tnewi.meta['comments'] = []   # weirdly duplicating the comments in the metadata many times which then gets hstacked for each source and hangs on writing the final table!
    tcat = tcat[catkeep]
    for gcati in gcatind:
        tnew = tnewi.copy()
        for col in catcols:
            if col == 'Source_Name':
                print(gcati, tgcat['RA'][gcati],tgcat['DEC'][gcati], sourcename(tgcat['RA'][gcati],tgcat['DEC'][gcati]))
                tnew[col] = sourcename(tgcat['RA'][gcati],tgcat['DEC'][gcati])
            elif col == 'Ngaus':
                tnew[col] = 1
            elif col == 'Size':
                tnew[col] = 2.*tgcat['DC_Maj'][gcati]
            elif col == 'Width':
                tnew[col] = 2.*tgcat['DC_Min'][gcati]
            elif col == 'Angle':
                tnew[col] = tgcat['DC_PA'][gcati]
            elif col in ['mcode',  'Artefact',  'handled', 'confused', 'HBA_Deep_flag',  'opt_flag',  'HBA_Deep_Name',  'HBA_Deep_Sep', 'HBA_Deep_ALPHA_J2000','HBA_Deep_DELTA_J2000','HBA_DeepWide_id','HBA_DeepWide_Sep','Deblended_From']:
                # ok to keep defaults here
                pass
            elif col in gcatcols:
                tnew[col] = tgcat[col][gcati]
            else:
                print ('col ',col,' not handled, quitting')
                sys.exit()
        tnew['S_Code'] = 'G'
        tnew['Deblended_From'] = g
        
        tcat = vstack((tcat, tnew)) 
        print('Added source',tnew['Source_Name'][0])
        newgaussrcs.append(tnew['Source_Name'][0])
            
    
    
print('cat now has ',len(tcat),'entries')
print('added srcs:',newgaussrcs)
    
    
    
    
    
    
# get the sources to merge
mergelist = []
mlist = glob.glob('LBABOO*mergelist')
for m in mlist:
    with open(m) as f:
        merge_names = f.readline().split()
        merge_names.sort()
    if merge_names not in mergelist:
        mergelist.append(merge_names)
mergelist = np.array(mergelist)
if len(np.unique(mergelist)) != len(mergelist.flatten()):
    print (mergelist)
    print ('mergelist is not unique')
    sys.exit(1)
    
    
    
    
# merge sources
for mm in mergelist:
    selind = []
    for mmi in mm:
        selind.append(np.where(tcat['Source_Name']==mmi)[0][0])
        
    clist = tcat[selind]
    
    
    tfluxsum = np.sum(clist['Total_flux'])
    ra = np.sum(clist['RA']*clist['Total_flux'])/tfluxsum
    dec = np.sum(clist['DEC']*clist['Total_flux'])/tfluxsum
    sname = sourcename(ra,dec)
    print('New sourcename is',sname, end=' ')
    tcat['Source_Name'][selind] = sname
    tcat['RA'][selind] = ra
    tcat['DEC'][selind] = dec
    tcat['E_RA'][selind] = np.sqrt(np.mean(clist['E_RA']**2.0))
    tcat['E_DEC'][selind] = np.sqrt(np.mean(clist['E_DEC']**2.0))
    tcat['Total_flux'][selind] = np.sum(clist['Total_flux'])
    tcat['E_Total_flux'][selind] = np.sqrt(np.sum(clist['E_Total_flux']**2.0))
    maxpk = np.argmax(clist['Peak_flux'])
    tcat['Peak_flux'][selind] = clist[maxpk]['Peak_flux']
    tcat['E_Peak_flux'][selind] = clist[maxpk]['E_Peak_flux']
    tcat['S_Code'][selind] = 'Z'
    tcat['Ngaus'][selind] = np.sum(clist['Ngaus'])
    
    if len(np.unique(tcat['mcode'][selind])) > 1:
        print('setting to complex match')
        tcat['mcode'][selind] = 2
    
    tcat['Isl_rms'][selind] = np.mean(clist['Isl_rms'])
    ms = Make_Shape(clist)
    tcat['Size'][selind] = ms.length()
    if tcat['Size'][0]>1000:
        print(clist['Source_Name'])
        print('Unreasonable source size detected')
    tcat['Width'][selind] = ms.width()
    tcat['Angle'][selind] = ms.pa()
    for k in ['Maj','Min','PA','E_Maj','E_Min','E_PA','DC_Maj','DC_Min','DC_PA','E_DC_Maj','E_DC_Min','E_DC_PA']:
        tcat[k][selind] = np.nan
        
    with open(sname+'.components','w') as f:
        for s in clist['Source_Name'][:]:
            f.write(s+'\n')
            
            # update source name in component catalogue
            if s in tcompcat['Component_Name']:
                icomp = np.where(tcompcat['Component_Name']==s)[0][0]
                tcompcat['Source_Name'][icomp] = sname
            elif s in tgcat['Component_Name']:
                icomp = np.where(tgcat['Component_Name']==s)[0][0]
                
                tcompcat.add_row(tcompcat[0])
                for col in tcompcat.colnames:
                    if col in gcatcols:
                        tcompcat[col][-1] = tgcat[col][icomp]
                    #elif col =='Source_id':
                        #tcompcat[col][-1] = max(tcompcat[col])+1
                    #elif col =='Isl_id':
                        #tcompcat[col][-1] = -1
                    elif col =='Ngaus':
                        tcompcat[col][-1] = 1
                    else:
                        print (col)
                        sys.exit()
            else:
                print ('where does this component come from?', s)
                sys.exit()
    
    
    print(selind)
    
print(len(tcat), len(np.unique(tcat['Source_Name'])))
_,keepind = np.unique(tcat['Source_Name'],return_index=True)
tcat = tcat[keepind]
print(len(tcat), len(np.unique(tcat['Source_Name'])))

#sys.exit()



Cdeepopt = SkyCoord(tdeepcato['RA'],tdeepcato['DEC'],unit='deg')
Cdeep = SkyCoord(tdeepcat['RA'],tdeepcat['DEC'],unit='deg')




C = SkyCoord(tcat['RA'],tcat['DEC'],unit='deg')
#tcat['opt_flag'] = optmoc.contains(C.ra,C.dec)
#tcat['opt_star_flag'] = omoc.contains(C.ra,C.dec)
#tcat['ir_star_flag'] = smoc.contains(C.ra,C.dec)

## change opt flag to false where source lies in both the star masks
#tcat['opt_flag'][tcat['opt_star_flag']  & tcat['ir_star_flag'] ] = False

tcat['opt_flag'] = filter_table(tcat,starmask)

print ('getting HBA matches for each source')
for i in range(len(tcat)):
    
    tcati = tcat[i]
    
    
    ra = tcati['RA']
    dec = tcati['DEC']
    
    C = SkyCoord(ra,dec,unit='deg')
    
    name = tcati['Source_Name']
    sid = tcati['Source_id']
    mcode = tcati['mcode']
    
    name = name.replace('LBABOO ','')
    print(name,end=' ')
    



    if mcode == 3:
        print('is artefact')
        tcat['Artefact'][i] = True
        tcat['handled'][i] = True
    elif mcode == 4:
        print('no possible match')
        tcat['HBA_Deep_flag'][i] = False
        tcat['handled'][i] = True
    elif mcode == 0:
        print('no match')
        tcat['handled'][i] = True
    elif mcode == 5:
        print('multiple matches')
        tcat['handled'][i] = True
        tcat['confused'][i] = True
        
        # in the deep optical area we record all the possible ids within the size of the LBA source in a hba_deep_ids file for each source
        if tcat['opt_flag'][i]:
            csep = C.separation(Cdeepopt)
            imin = np.where(csep <  tcat['Maj'][i]*u.deg)[0]
            with open(name+'.hba_deep_ids','w') as f:
                for ii in imin:
                    f.write(tdeepcato['Source_Name'][ii]+'\t'+str(csep[ii].to(u.arcsec).value)+'\n')
        
    elif mcode == 1: # simple match
        if tcat['opt_flag'][i]:
            csep = C.separation(Cdeepopt)
            if os.path.isfile(name+'.hba_deep_id'):
                with open(name+'.hba_deep_id','r') as f:
                    m = f.readline().strip()
                    if m not in tdeepcato['Source_Name']:
                        print ('By-hand DEPP HBA NAME', m,' not found in catalouge')
                        sys.exit()
                imin = np.where(tdeepcato['Source_Name'] == m)[0][0]
                tcat['HBA_Deep_Name'][i]  = tdeepcato['Source_Name'][imin]
                tcat['HBA_Deep_Sep'][i]  = csep[imin].to(u.arcsec).value
                tcat['HBA_Deep_ALPHA_J2000'][i]  = tdeepcato['ALPHA_J2000'][imin]
                tcat['HBA_Deep_DELTA_J2000'][i]  = tdeepcato['DELTA_J2000'][imin]
                print ('by hand match to HBA DEEP opt source', tcat['HBA_Deep_Name'][i], csep[imin].to(u.arcsec))
                    
            else:
                imin = csep.argmin()
                tcat['HBA_Deep_Name'][i]  = tdeepcato['Source_Name'][imin]
                tcat['HBA_Deep_Sep'][i]  = csep[imin].to(u.arcsec).value
                tcat['HBA_Deep_ALPHA_J2000'][i]  = tdeepcato['ALPHA_J2000'][imin]
                tcat['HBA_Deep_DELTA_J2000'][i]  = tdeepcato['DELTA_J2000'][imin]
                print ('match to nearest, HBA DEEP opt source', tcat['HBA_Deep_Name'][i], csep[imin].to(u.arcsec))
                tcat['handled'][i] = True
                if np.sum(csep < tcat['Size'][i]*u.arcsec) > 1:
                    print (np.sum(csep < tcat['Size'][i]*u.arcsec), '*multiple matches*')
                    
                    imin = np.where(csep <  tcat['Size'][i]*u.arcsec)[0]
                    with open(name+'.hba_deep_ids','w') as f:
                        for ii in imin:
                            f.write('{name:s}\t{s:.2f}\n'.format(name=tdeepcato['Source_Name'][ii],s=csep[ii].to(u.arcsec).value))
        else:
            csep = C.separation(Cdeep)
            imin = csep.argmin()
            tcat['HBA_DeepWide_id'][i]  = tdeepcat['Source_id'][imin]
            tcat['HBA_DeepWide_Sep'][i]  = csep[imin].to(u.arcsec).value
            print ('match to nearest HBA DEEP Wide source', tcat['HBA_DeepWide_id'][i], csep[imin].to(u.arcsec))
            tcat['handled'][i] = True
            
        
    elif mcode == 2:
        if tcat['opt_flag'][i]:
            print('complex HBA DEEP opt TODO')
        else:
            print('complex HBA DEEP wide TODO')



    if os.path.exists(name+'.matchlist'):
        continue
    if name in mlist:
        break
    if os.path.exists(name+'.complex'):
        continue

# exclude artefacts
tcat = tcat[~tcat['Artefact']]


## add AGN classifications
hcatopt = Table.read(deepids)


tcat.add_column(Column(data=np.zeros(len(tcat),dtype=bool),name='has_hba_opt'))
tcat.add_column(Column(data=-99*np.ones(len(tcat),dtype=int),name='AGN_final'))
tcat.add_column(Column(data=-99*np.ones(len(tcat),dtype=int),name='RadioAGN_final'))
tcat.add_column(Column(data=-1*np.ones(len(tcat),dtype=float),name='Z_BEST'))
tcat.add_column(Column(data=-1*np.ones(len(tcat),dtype=float),name='SFR_conc'))
tcat.add_column(Column(data=-1*np.ones(len(tcat),dtype=float),name='Mass_conc'))
tcat.add_column(Column(data=-1*np.ones(len(tcat),dtype=float),name='HBA_Size'))
tcat.add_column(Column(data=-1*np.ones(len(tcat),dtype=float),name='HBA_RA'))
tcat.add_column(Column(data=-1*np.ones(len(tcat),dtype=float),name='HBA_DEC'))
tcat.add_column(Column(data=-1*np.ones(len(tcat),dtype=float),name='HBA_Total_flux'))
tcat.add_column(Column(data=-1*np.ones(len(tcat),dtype=float),name='HBA_E_Total_flux'))
for i in range(len(tcat)):
    if len(tcat['HBA_Deep_Name'][i]) > 0:
        hi = np.where(tcat['HBA_Deep_Name'][i]==hcatopt['Source_Name_1'])[0][0]
        tcat['has_hba_opt'][i] = 1
        tcat['AGN_final'][i] = hcatopt['AGN_final'][hi]
        tcat['RadioAGN_final'][i] = hcatopt['RadioAGN_final'][hi]
        tcat['Z_BEST'][i] = hcatopt['Z_BEST'][hi]
        tcat['SFR_conc'][i] = hcatopt['SFR_conc'][hi]
        tcat['Mass_conc'][i] = hcatopt['Mass_conc'][hi]
        tcat['HBA_Total_flux'][i] = hcatopt['Total_flux'][hi]
        tcat['HBA_E_Total_flux'][i] = hcatopt['E_Total_flux'][hi]
        tcat['HBA_RA'][i] = hcatopt['RA'][hi]
        tcat['HBA_DEC'][i] = hcatopt['DEC'][hi]
        if np.isfinite(hcatopt['LGZ_Size'][hi]):
            tcat['HBA_Size'][i] = hcatopt['LGZ_Size'][hi]
        else:
            tcat['HBA_Size'][i] = 2*hcatopt['DC_Maj'][hi]
            
            
           
tcat['flag_SFG'] = (tcat['AGN_final']==0)  &  (tcat['RadioAGN_final']==0)     #-> star-forming galaxy
tcat['flag_RQ'] = (tcat['AGN_final']==1)  &  (tcat['RadioAGN_final']==0 )    #-> 'radio-quiet' AGN
tcat['flag_LERG'] = (tcat['AGN_final']==0)  &  (tcat['RadioAGN_final']==1)     #-> 'jet-mode' radio AGN / LERG
tcat['flag_HERG'] = (tcat['AGN_final']==1)  &  (tcat['RadioAGN_final']==1 )    #-> quasar-like radio AGN / HERG
tcat['flag_AGN'] = tcat['flag_LERG'] | tcat['flag_HERG'] | tcat['flag_RQ']
tcat['flag_unc'] = (tcat['AGN_final']==-1) | (tcat['RadioAGN_final']==-1)    #-> no secure classifitcation 


# write final
tcat.write(catout, overwrite=True)
tcompcat.write(compcatout, overwrite=True)

    
# clean selection is
#~Artefact & opt_flag & ~opt_star_flag & ~ir_star_flag

