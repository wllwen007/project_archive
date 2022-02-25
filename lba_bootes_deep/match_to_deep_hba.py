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




def display_cat_ellipse(frameno,cat,color='white'):
    print('loading catalogue into frame '+str(frameno) )
    ds9.set('frame '+str(frameno))
    for cati in cat:
        S="{csys}; ellipse {ra} {dec} {size} {size2} {pa} # color = {col}".format(csys='fk5',ra=cati['RA'], dec=cati['DEC'], size=cati['Maj'], size2=cati['Min'], pa=cati['PA']+90,col=color)
        ds9.set('regions ',S)
        
    return

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
    ds9.set('scale limits -0.001 0.01')
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


fcodelist = 'lba_match_deephba_codes.fits'
if not os.path.isfile(fcodelist):
    
    tmcode = Table((Column(tcat['Source_Name']),
                    Column(data=99*np.ones(len(tcat),dtype=int), name='match_code')))
else:
    tmcode = Table.read(fcodelist)


all_merge = []

# for loop but we can step backwards if needed
i = 0
backward = False
while True:
    print(i)
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

    if tmcode['match_code'][tmcode['Source_Name'] == name] != 99:
        if not backward:  ## don't do this if its a rewind
            i+=1
            continue
    backward = False
    
    curr_code = tmcode['match_code'][tmcode['Source_Name'] == name][0]

    
    #for frame in ['1','2']:
        #ds9.set('frame '+frame)
    ds9.set('pan to {ra} {dec} wcs fk5 degrees'.format(ra=ra,dec=dec))
    ds9.set('crosshair {ra} {dec} wcs fk5 degrees'.format(ra=ra,dec=dec))
    
    name = name.replace('LBABOO ','')
    
    
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
    
    
    print('Source ({i}) {name} - {c}'.format(i=i,name=name,c=curr_code))
    print('(0) no match')
    print('(1) simple nearest match ok - default')
    print('(2) complex')
    print('(3) lba source is artefact')
    print('(4) no match possible')
    print('(l/u) left/up go to previous')
    print('(q) quit')
    
    while True:
        s = input('input: ')
        s = s.lower()
        if s =='0':
            hasmatch = False
        if s =='q':
            tmcode.write(fcodelist, overwrite=True)
            sys.exit()
        if s == '':
            # default
            s = '1'
        if s in ['0','1','2','3','4','l','u']:
            break
        else:
            print ('enter valid value')
            
    if s in ['l','u']:
        i -= 1
        backward = True
        continue
    tmcode['match_code'][tmcode['Source_Name'] == name] = int(s)
        
    i+=1
    
tmcode.write(fcodelist, overwrite=True)
    


