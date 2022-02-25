import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import aplpy as ap
import utils.plotting as pp
import utils.make_subim as mm
import utils.cutouts as cc
import astropy.io.fits as fits
import astropy.units as u
from astropy.table import Table, Column, vstack
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky

# lotss-catalogue  # add to PYTHONPATH
#export PYTHONPATH=/home/wwilliams/scripts/git/lotss-catalogue/utils:$PYTHONPATH
from subim import extract_subim
from overlay import show_overlay
from separation import separation
from download_image_files import get_legacy,get_first,get_wise



import mocpy

def mogrify(filename):
    command='mogrify -quality 90 -trim '+filename
    p = subprocess.Popen(command,shell=True)
    return p

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



lofarfile=fits.open('/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/image_full_phase_m.NS_shift.app.facetRestored.blanked.fits')
lofarhbafile=fits.open('/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits')

spitzerfile=fits.open('/net/beerze/data2/wwilliams/projects/lofar_surveys/deepfields/Bootes_optical/SDWFS/I2_bootes.v32.fits')
ibandfile=fits.open('/net/beerze/data2/wwilliams/projects/lofar_surveys/deepfields/Bootes_merged_optical/Bootes_iband.fits')
    
    
imask = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes_merged_optical/radio_optical_overlap_masks/image_full_ampphase_di_m.NS_shift.blanked.rescaled.rms_optmask.fits'
spmask = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes_merged_optical/radio_optical_overlap_masks/image_full_ampphase_di_m.NS_shift.blanked.rescaled.rms_spmask.fits'
    
ot=Table.read('/net/beerze/data2/wwilliams/projects/lofar_surveys/deepfields/Bootes_LR/Bootes_ML_RUN_fin_overlap_srl_workflow_fixed.fits')
gals=Table.read('/net/beerze/data2/wwilliams/projects/lofar_surveys/deepfields/Bootes_merged_optical/Bootes_merged_pos.fits')


def get_pixel_values(ra,dec,image):
  import astropy.io.fits as pf
  import astropy.wcs as pw
  head = pf.getheader(image)
  data = pf.getdata(image)
  
  # Parse the WCS keywords in the primary HDU
  wcs = pw.WCS(head)

  # Some pixel coordinates of interest.
  #skycrd = np.array([ra,dec])
  if wcs.naxis ==4:
    x,y,_,_ = wcs.all_world2pix(ra,dec,0*ra,0*ra, 1)
    x=np.array(np.round(x,0), dtype=int)
    y=np.array(np.round(y,0), dtype=int)
    if len(data.shape) ==4:
        values = data[0,0,y,x]
    else:
        values = data[y,x]
  elif all_world2pix ==2 :
    x,y,_,_ = wcs.all_world2pix(ra,dec, 1)
    x=np.array(np.round(x,0), dtype=int)
    y=np.array(np.round(y,0), dtype=int)
    values = data[y,x]
      
  return values

icovered = np.isfinite(get_pixel_values(tcat['RA'],tcat['DEC'],  imask))
spcovered = np.isfinite(get_pixel_values(tcat['RA'],tcat['DEC'],  spmask))

size=180 / 3600.
scale=1
pwd = os.getcwd()
for i in range(len(tcat)):
    
    
    
    r = tcat[i]
    
    
    
    ra = r['RA']
    dec = r['DEC']
    name = r['Source_Name']
    sid = r['Source_id']
    peak = r['Peak_flux']
    sourcename = r['Source_Name']

    title = None

    iimage = 'cutouts/'+sourcename+'_I.png'
    iimaged = 'cutouts/'+sourcename+'_Id.png'
    ipimage = 'cutouts/'+sourcename+'_Ip.png'
    simage = 'cutouts/'+sourcename+'_S.png'
    spimage = 'cutouts/'+sourcename+'_Sp.png'
    
    #gseps = separation(ra,dec,gals['ra'],gals['dec'])
    #marker_ra = gals['ra'][gseps<size]
    #marker_dec = gals['dec'][gseps<size]
    # don't show optical galaxies as markers
    marker_ra = None
    marker_dec = None
    
    seps = separation(ra,dec,tcat['RA'],tcat['DEC'])
    ots = tcat[seps<size]
    
    dseps = separation(ra,dec,tdeepcat['RA'],tdeepcat['DEC'])
    dts = tdeepcat[dseps<size]
    
    if os.path.isfile(iimaged):
        print (name,'done')
        continue
    
    if spcovered[i]:
        print(spitzerfile,ra,dec,size)
        shdu=extract_subim(spitzerfile,ra,dec,size)
    else:
        print(i,'out of spitzer coverage, using wise')
        os.chdir('downloads')
        wisename=get_wise(r['RA'],r['DEC'],1)
        os.chdir(pwd)
        wisefile='./downloads/'+wisename
        print('Using wise image',wisefile)
        wisehdu=fits.open(wisefile)
        shdu = wisehdu
    if icovered[i]:
        ihdu=extract_subim(ibandfile,ra,dec,size)
    else:
        print(i,'out of iband coverage, using legacy')
        os.chdir('downloads')
        legacyname=get_legacy(r['RA'],r['DEC'],bands='zrg')
        os.chdir(pwd)
        optfile='./downloads/'+legacyname
        print('Using optical image',optfile)
        pshdu=fits.open(optfile)
        if pshdu[0].header['NAXIS']==0:
            print('*** No optical image! ***')
            logfile.write('*** No optical %s %f %f ***\n' % (sourcename,r['RA'],r['DEC']))
            continue
        # nan-blank
        pshdu[0].data=np.where(pshdu[0].data>8,np.nan,pshdu[0].data)
        ihdu = pshdu
        
        
    lhhdu=extract_subim(lofarhbafile,ra,dec,size)
    lhdu=extract_subim(lofarfile,ra,dec,size)
    
    pg=gals[(np.abs(gals['ra']-ra)<(size/np.cos(dec*np.pi/180.0))) & (np.abs(gals['dec']-dec)<size)]
    ls=[]
    ccol=[]
    for nr in ots:
        if nr['Source_Name']==r['Source_Name']:
            ls.append('solid')
        else:
            ls.append('dashed')
        ccol.append('r')


    plist=[]
    show_overlay(lhdu,ihdu,ra,dec,size,firsthdu=lhhdu,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=iimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,peak=peak,noisethresh=1)
    plist.append(mogrify(iimage))


    for nr in dts:
        ls.append('dashed')
        ccol.append('b')
    scols = ['RA','DEC','Maj','Min','PA']
    otsdts = vstack((ots[scols],dts[scols]))
    
    plist=[]
    show_overlay(lhdu,ihdu,ra,dec,size,firsthdu=lhhdu,overlay_cat=otsdts,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=iimaged,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,ellipse_color=ccol,peak=peak,noisethresh=1)
    plist.append(mogrify(iimaged))
    
    #show_overlay(lhdu,ihdu,ra,dec,size,firsthdu=lhhdu,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=ipimage,show_lofar=False,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,peak=peak,noisethresh=1,plotpos=pg,ppsize=350)
    #plist.append(mogrify(ipimage))
    
    show_overlay(lhdu,shdu,ra,dec,size,firsthdu=lhhdu,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=simage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,peak=peak)
    plist.append(mogrify(simage))
    
    #show_overlay(lhdu,shdu,ra,dec,size,firsthdu=lhhdu,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=spimage,show_lofar=False,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,peak=peak)
    #plist.append(mogrify(spimage))
       
    for p in plist:
        p.wait()     
    #sys.exit()

