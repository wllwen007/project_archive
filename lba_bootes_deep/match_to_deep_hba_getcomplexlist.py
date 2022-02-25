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




def display_cat_ellipse(frameno,cat,color='white',width=1,dash=False):
    print('loading catalogue into frame '+str(frameno) )
    ds9.set('frame '+str(frameno))
    if dash:
        dashed = 'dash=1'
    else:
        dashed = ''
    for cati in cat:
        S="{csys}; ellipse {ra} {dec} {size} {size2} {pa} # color = {col} width={w} {d}".format(csys='fk5',ra=cati['RA'], dec=cati['DEC'], size=cati['Maj'], size2=cati['Min'], pa=cati['PA']+90,col=color,w=width,d=dashed)
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
    display_cat_ellipse(1,tgcat,'blue',dash=True)
    display_cat_ellipse(1,tcat, 'green')

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


call = SkyCoord(tcat['RA'],tcat['DEC'],unit='deg')
callg = SkyCoord(tgcat['RA'],tgcat['DEC'],unit='deg')
calldeep = SkyCoord(tdeepcat['RA'],tdeepcat['DEC'],unit='deg')

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
    
    name = name.replace('LBABOO ','')
    
    mcode = tmcode[tmcode['Source_Name'] == name]['match_code'][0]

    if os.path.exists(name+'.matchlist'):
        i += 1
        continue
    if os.path.exists(name+'.mergelist'):
        i += 1
        continue
    if os.path.exists(name+'.complex'):
        with open(name+'.complex') as f:
            l = f.readlines()
        if len(l) > 0:
            i += 1
            continue
    if mcode !=  2: 
        i += 1
        continue
    
    #for frame in ['1','2']:
        #ds9.set('frame '+frame)
    ds9.set('pan to {ra} {dec} wcs fk5 degrees'.format(ra=ra,dec=dec))
    ds9.set('crosshair {ra} {dec} wcs fk5 degrees'.format(ra=ra,dec=dec))
    
    
    
    tgcati = tgcat[tgcat['Source_id']==sid]
    
    c = SkyCoord(ra,dec,unit='deg')
    
    
    sep = c.separation(call)
    snear = (sep < 0.25*u.deg) & (sep > 0*u.deg)
    tcatn = tcat[snear]
    
    dsep = c.separation(calldeep)
    sdmin = dsep.argmin()
    tdeepcatn = tdeepcat[sdmin]
    
    
    print(f'Source ({i}) {name}')
    for ii in np.where(snear)[0]:
        if sep[ii].to(u.arcsec) < (10.*u.arcsec):
            print('{ii} {name} {size:.1f} {sep:.1f}'.format(ii=ii,name=tcat['Source_Name'][ii], size=tcat['Maj'][ii]*3600, sep=sep[ii].to(u.arcsec)) )
    
    print('(c) complex - split lba source')
    print('(m) complex - split multiple lba source')
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
        sources = []
        addsources = True
        while addsources:
            print('(a)dd source ')
            print('(d)one')
            s2 = input('input: ')
            s2 = s2.lower()
            if s2 == 'a':
                print('mark gaus components')
                csel = get_coords()
                isel = []
                for cc in csel:
                    csep = cc.separation(callg)
                    isel.append(csep.argmin())
                tsel = tgcat[isel]
                print('selected sources are now in cyan')
                display_cat_ellipse(1,tsel,color='cyan',width=3)
                gsource = list(tsel['Gaus_id'][:])
                sources.append(gsource)
            elif s2 == 'd':
                addsources = False
        print(sources)
        with open(name+'.split','w') as f:
            for ss in sources:
                for ssi in ss:
                    f.write(f'{ssi} ')
                f.write('\n')
            
    elif s == 'm':
        print ('complex source') 
        sources = []
        lbasources = []
        addsources = True
        while addsources:
            print('(a)dd source ')
            print('(d)one')
            s2 = input('input: ')
            s2 = s2.lower()
            if s2 == 'a':
                print('mark lba source')
                csel = get_coords()
                isel = []
                for cc in csel:
                    csep = cc.separation(call)
                    isel.append(csep.argmin())
                tsel = tcat[isel]
                source = [list(tsel['Source_Name'][:])]
                lbasources.append(source)
                print('mark gaus components')
                csel = get_coords()
                isel = []
                for cc in csel:
                    csep = cc.separation(callg)
                    isel.append(csep.argmin())
                tsel = tgcat[isel]
                print('selected sources are now in cyan')
                display_cat_ellipse(1,tsel,color='cyan',width=3)
                source = [list(tsel['Gaus_id'][:])]
                sources.append(source)
            elif s2 == 'd':
                addsources = False
        print(lbasources)
        print(sources)
            
            
        
        with open(name+'.complex','w') as f:
            f.write(name)
    elif s =='q':
        sys.exit()
    #if s == '':
    else:
        print ('continuing')   
            
    i+=1
    
    


