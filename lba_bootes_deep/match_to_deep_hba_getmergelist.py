import os
import sys
import numpy as np
import astropy.units as u
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import pyds9




clobber = True

dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'
imfits = 'image_full_phase_m.NS_shift.app.facetRestored.blanked.fits'

cat = 'bootes_deep_lba_hbashift.cat.fits'
tcat = Table.read(dpath+cat)


gcat = 'bootes_deep_lba_hbashift.gcat.fits'
tgcat = Table.read(dpath+gcat)

deepcat= '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.blanked.scaled.cat.fits'
deepim = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits'

mocpath = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes_merged_optical/final_moc'

mocpath+''

deepcato= '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/bootes_final_cross_match_catalogue-v1.0.fits'

tdeepcat = Table.read(deepcat)
tdeepcato = Table.read(deepcato)




def display_cat_ellipse(frameno,cat,color='white',width=1):
    print('loading catalogue into frame '+str(frameno) )
    ds9.set('frame '+str(frameno))
    for cati in cat:
        S="{csys}; ellipse {ra} {dec} {size} {size2} {pa} # color = {col} width={w}".format(csys='fk5',ra=cati['RA'], dec=cati['DEC'], size=cati['Maj'], size2=cati['Min'], pa=cati['PA']+90,col=color,w=width)
        ds9.set('regions ',S)
        
    return

def get_coords():
    ra = [] ; dec = []
    print ('select component sources in frame 1')
    print ('press \'q\' to quit')
    while True:
        crd = ds9.get('imexam key coordinate fk5')
        cc = crd.split()
        key = cc[0]
        if key.lower() == 'q':
            break
        else:
            ra.append(float(cc[1]))
            dec.append(float(cc[2]))
            print (cc[1], cc[2])
        
    c = SkyCoord(ra,dec,unit='deg')
    return c

exists = False
t = pyds9.ds9_targets()
if t is not None:
    for ti in t:
        if ti.split()[0].split(':')[1] == 'matcher':
            exists = True
            print('found ds9 session')

ds9 = pyds9.DS9('matcher',wait=300)
if not exists: 

    ds9.set('frame 1')
    ds9.set('tile yes')
    ds9.set('file '+dpath+imfits)
    ds9.set('zoom 2')
    ds9.set('scale limits -0.001 0.01')
    ds9.set('cmap Heat') 
    display_cat_ellipse(1,tcat)

    ds9.set('frame 2')
    ds9.set('file '+deepim)
    ds9.set('scale limits -0.0005 0.005')
    ds9.set('cmap Heat') 
    display_cat_ellipse(2,tdeepcat,'cyan')


ds9.set('frame lock wcs')
ds9.set('lock crosshair wcs')


Ncat = len(tcat)

if 'HBA_Source_Name' in tcat.colnames:
    restart = True
else:
    tcat.add_column(Column(name='HBA_Source_Name',data=np.zeros(Ncat, dtype='S100')))
    tcat.add_column(Column(name='HBA_deep_cover',data=np.zeros(Ncat, dtype=bool)))
    restart = False


fcodelist = 'lba_match_deephba_codes.txt'
if not os.path.isfile(fcodelist):
    print('need match codes: run match_to_deep_hba.py')
    
tmcode = Table.read(fcodelist,format='ascii')


all_merge = []

# for loop but we can step backwards if needed
i = 0
backward = False
while True:
    if i >= len(tcat): break
    
    
    #if i<550:
        #continue
    
    tcati = tcat[i]
    
    # only make the m sources - check for deblending...
    #if tcati['S_Code'] != 'S':
        #continue
    
    
    ra = tcati['RA']
    dec = tcati['DEC']
    name = tcati['Source_Name']
    sid = tcati['Source_id']
    
    #'LBABOOJ143656.19+341827.5'
    if name not in ['LBABOOJ142440.51+323605.9',
'LBABOOJ142459.44+352842.7',
'LBABOOJ142505.39+334652.5',
'LBABOOJ142526.50+344922.1',
'LBABOOJ142531.13+353343.1',
'LBABOOJ142531.79+353340.1',
'LBABOOJ142616.46+342121.1',
'LBABOOJ142639.46+350242.3',
'LBABOOJ142649.68+332931.8',
'LBABOOJ142702.71+333339.4',
'LBABOOJ142703.01+333435.2',
'LBABOOJ142706.50+333512.4',
'LBABOOJ142707.12+333428.9',
'LBABOOJ142755.95+332141.7',
'LBABOOJ142756.20+332141.5',
'LBABOOJ142759.95+324739.6',
'LBABOOJ142806.62+325935.3',
'LBABOOJ142807.39+325932.0',
'LBABOOJ142838.55+342341.1',
'LBABOOJ142842.26+342447.3',
'LBABOOJ142842.34+354325.1',
'LBABOOJ142842.68+354324.5',
'LBABOOJ142845.33+323837.4',
'LBABOOJ142848.10+323940.6',
'LBABOOJ142848.80+323822.2',
'LBABOOJ142902.81+353823.1',
'LBABOOJ142919.13+352813.0',
'LBABOOJ142919.22+324633.0',
'LBABOOJ142922.11+324626.5',
'LBABOOJ142927.98+352956.2',
'LBABOOJ142949.23+351009.6',
'LBABOOJ142949.52+353440.3',
'LBABOOJ142954.12+343517.5',
'LBABOOJ142956.47+325518.4',
'LBABOOJ143008.02+325833.5',
'LBABOOJ143008.24+331537.8',
'LBABOOJ143008.89+325840.0',
'LBABOOJ143012.28+331437.9',
'LBABOOJ143012.82+331433.1',
'LBABOOJ143014.55+345918.1',
'LBABOOJ143034.29+324321.6',
'LBABOOJ143039.15+352351.9',
'LBABOOJ143052.56+335559.9',
'LBABOOJ143053.48+330031.0',
'LBABOOJ143057.83+334501.3',
'LBABOOJ143106.60+334610.0',
'LBABOOJ143108.07+352557.5',
'LBABOOJ143108.11+352531.6',
'LBABOOJ143114.42+323225.5',
'LBABOOJ143115.19+334615.8',
'LBABOOJ143120.83+325236.7',
'LBABOOJ143122.00+353527.9',
'LBABOOJ143128.44+342726.3',
'LBABOOJ143134.71+351505.9',
'LBABOOJ143136.20+351442.0',
'LBABOOJ143143.19+353551.8',
'LBABOOJ143203.68+330744.4',
'LBABOOJ143217.30+351421.5',
'LBABOOJ143223.38+332909.7',
'LBABOOJ143227.46+354607.2',
'LBABOOJ143233.96+354535.0',
'LBABOOJ143235.24+354346.9',
'LBABOOJ143237.32+354326.1',
'LBABOOJ143239.68+323421.2',
'LBABOOJ143242.90+322047.2',
'LBABOOJ143252.09+322103.6',
'LBABOOJ143312.78+351431.0',
'LBABOOJ143329.43+343442.1',
'LBABOOJ143342.21+341135.8',
'LBABOOJ143404.64+331250.8',
'LBABOOJ143426.89+353204.8',
'LBABOOJ143429.72+342720.5',
'LBABOOJ143430.81+354214.5',
'LBABOOJ143431.46+350104.6',
'LBABOOJ143433.02+342752.0',
'LBABOOJ143433.17+352134.3',
'LBABOOJ143443.18+333007.1',
'LBABOOJ143443.26+332833.8',
'LBABOOJ143444.26+332912.9',
'LBABOOJ143445.94+332818.1',
'LBABOOJ143447.03+332823.2',
'LBABOOJ143448.35+354238.9',
'LBABOOJ143450.84+354251.2',
'LBABOOJ143602.97+334407.2',
'LBABOOJ143603.34+334350.9',
'LBABOOJ143627.51+352622.6',
'LBABOOJ143653.04+341658.5',
'LBABOOJ143657.97+355105.6',
'LBABOOJ143728.63+350749.1',
'LBABOOJ143737.93+342006.3',
'LBABOOJ143827.92+354108.1'  ]:
        i+=1
        continue
    
    name = name.replace('LBABOO ','')
    
    mcode = tmcode[tmcode['Source_Name'] == name]['match_code'][0]

    #if os.path.exists(name+'.matchlist'):
        #i += 1
        #continue
    #if os.path.exists(name+'.mergelist'):
        #i += 1
        #continue
    #if os.path.exists(name+'.complex'):
        #i += 1
        #continue
    #if mcode !=  2: 
        #i += 1
        #continue
    
    #for frame in ['1','2']:
        #ds9.set('frame '+frame)
    ds9.set('pan to {ra} {dec} wcs fk5 degrees'.format(ra=ra,dec=dec))
    ds9.set('crosshair {ra} {dec} wcs fk5 degrees'.format(ra=ra,dec=dec))
    
    
    
    tgcati = tgcat[tgcat['Source_id']==sid]
    
    c = SkyCoord(ra,dec,unit='deg')
    call = SkyCoord(tcat['RA'],tcat['DEC'],unit='deg')
    calldeep = SkyCoord(tdeepcat['RA'],tdeepcat['DEC'],unit='deg')
    
    
    sep = c.separation(call)
    snear = (sep < 0.25*u.deg) & (sep > 0*u.deg)
    tcatn = tcat[snear]
    
    dsep = c.separation(calldeep)
    sdmin = dsep.argmin()
    tdeepcatn = tdeepcat[sdmin]
    
    
    print(f'Source ({i}) {name}')
    
    print('(c) complex')
    print('(h) multiple hba matches')
    print('(l) multiple lba matches')
    print('(n) multiple lba matches - nested')
    print('(q) quit')
    
        
    s = input('input: ')
    s = s.lower()
    
    
    if s =='l':
        csel = get_coords()
        isel = []
        for cc in csel:
            csep = cc.separation(call)
            isel.append(csep.argmin())
        tsel = tcat[isel]
        print('selected sources are now in blue')
        display_cat_ellipse(1,tsel,color='blue',width=3)
        with open(name+'.mergelist','w') as f:
            f.write(' '.join(tsel['Source_Name'][:]))
        print(', '.join(tsel['Source_Name'][:]))
    elif s =='n':
        csel = get_coords()
        isel = []
        for cc in csel:
            csep = cc.separation(call)
            isel.append(csep.argsort()[:2])  # take first 2 - assuming only 2 nested srcs
        tsel = tcat[isel]
        print('selected sources are now in blue')
        display_cat_ellipse(1,tsel,color='blue',width=3)
        with open(name+'.mergelist','w') as f:
            f.write(' '.join(tsel['Source_Name'][:]))
        print(', '.join(tsel['Source_Name'][:]))
    elif s =='h':
        csel = get_coords()
        isel = []
        for cc in csel:
            csep = cc.separation(calldeep)
            isel.append(csep.argmin())  
        tsel = tdeepcat[isel]
        print('selected sources are now in blue')
        display_cat_ellipse(2,tsel,color='blue',width=3)
        with open(name+'.matchlist','w') as f:
            for ss in tsel['Source_id'][:]:
                f.write(str(ss)+' ')
        #print(', '.join(tsel['Source_id'][:]))
        
    elif s == 'c':
        print ('complex source')   
        with open(name+'.complex','w') as f:
            f.write(name)
    elif s =='q':
        sys.exit()
    #if s == '':
    else:
        print ('continuing')   
            
    i+=1
    
    


