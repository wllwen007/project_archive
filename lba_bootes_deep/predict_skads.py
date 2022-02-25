from astropy.table import Table
import numpy as np
import astropy.io.fits as pf
import astropy.wcs as pw

sft = ['None', 'Quiescent', 'Starburst']
agnt = ['None', 'RQ',' FRI','FRII','GPS']

Area = 20.*20   # sq deg

scat = Table.read('/net/beerze/data2/wwilliams/projects/skads/wilman_cat_all.fits')

cat = Table.read('/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/bootes_deep_lba_hbashift.cat.fits')
rmsmap = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/bootes_deep_lba_hbashift.rms.fits'
head = pf.getheader(rmsmap)
rmsdat = pf.getdata(rmsmap)
wcs = pw.WCS(head)


optmap = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/data_release/bootes/rms_starmask_optical.fits'
ohead = pf.getheader(optmap)
odat = pf.getdata(optmap)
owcs = pw.WCS(ohead)


S151 = 10**scat['itot_151']
S610 = 10**scat['itot_610']
alp = np.log10(S151/S610) / np.log10(151./610)
S54 = S151*(54./151)**alp

# 0.25 mJy flux limit
# 0.7 mJy flux limit
sel = (S54 > 0.7*5*1e-3)
#scat = scat[sel]

#sel = (scat['right_ascension'])
#scat = scat[sel]

scat['right_ascension'] = scat['right_ascension'] +217.0
scat['declination'] = scat['declination'] +35

x,y,_,_ = wcs.all_world2pix(scat['right_ascension'], scat['declination'], 0*scat['declination'], 0*scat['declination'], 0)   
_,_,Nx,Ny = rmsdat.shape

x = np.array(x, dtype=int)
y = np.array(y, dtype=int)
x[x<0] = 0
x[x>=Nx] = Nx-1
y[y<0] = 0
y[y>=Ny] = Ny-1




rmsobs = rmsdat[0,0,x,y]
detected = S151>(5*rmsobs)



x,y,_,_ = owcs.all_world2pix(scat['right_ascension'], scat['declination'], 0*scat['declination'], 0*scat['declination'], 0)   
Nx,Ny = odat.shape

x = np.array(x, dtype=int)
y = np.array(y, dtype=int)
x[x<0] = 0
x[x>=Nx] = Nx-1
y[y<0] = 0
y[y>=Ny] = Ny-1


oobs = odat[x,y]
ocover = np.isfinite(oobs)


for mask in [detected, detected&ocover]:
    scati = scat[mask]


    Area_norm = 30./ Area
    Area_norm = 1.
    print('SF :', np.sum(scati['sftype']!=0)*Area_norm)
    #i = 1
    #for sftt in sft[1:]:
        #sel = scati['sftype']==i
        #print('SF - {}: {:.0f}'.format(sftt, np.sum(sel)*Area_norm))
        ##print(np.median(alp[mask][sel]))
        #i+=1  
    print('AGN :', np.sum(scati['agntype']!=0)*Area_norm)
    #i = 1
    #for agntt in agnt[1:]:
        #sel = scati['agntype']==i
        #print('AGN - {}: {:.0f}'.format(agntt, np.sum(sel)*Area_norm))
        ##print(np.median(alp[mask][sel]))
        #i+=1
        
    
    
