import os
import matplotlib.pyplot as plt
import numpy as np
import aplpy as ap
import utils.plotting as pp
import utils.make_subim as mm
import utils.cutouts as cc
import astropy.io.fits as pf
import astropy.units as u
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky
from astropy.convolution import Gaussian2DKernel, convolve
from reproject import reproject_interp


clobber = False
#clobber=True

dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'
imfits = 'image_full_phase_m.NS_shift.int.facetRestored.blanked.fits'

cat = 'bootes_deep_lba_hbashift.cat.fits'
cat ='bootes_deep_lba.cat.matched.v0.1.fits'
tcat = Table.read(dpath+cat)


gcat = 'bootes_deep_lba_hbashift.gcat.fits'
tgcat = Table.read(dpath+gcat)

deepcat= '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.blanked.scaled.cat.fits'
deepidcat= '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/bootes_final_cross_match_catalogue-v1.0.fits'
deepim = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits'

tdeepcat = Table.read(deepcat)
tdeepidcat = Table.read(deepidcat)


def std_sigclip(x, nit=10, nclip=5.):
    
    x = x[np.isfinite(x)]
    
    for i in range(nit):
        std = np.std(x)
        mn = np.mean(x)
        
        x = x[ (np.abs (x-mn) - nclip*std < 0)]
    
    mn = np.mean(x)
    std = np.std(x)
    
    
    return mn, std


def get_scale(t):
    try:
        m, s = std_sigclip(t[0].data)
    except:
        return 0, 10, 1
    vmin = -1.*s +m
    vmax = 25.*s + m
    vmid = 10.*s + m
    return vmin, vmax, vmid


call = SkyCoord(tcat['RA'],tcat['DEC'],unit='deg')
calldeep = SkyCoord(tdeepcat['RA'],tdeepcat['DEC'],unit='deg')
calldeepid = SkyCoord(tdeepidcat['RA'],tdeepidcat['DEC'],unit='deg')
for i in range(len(tcat)-1,-1,-1):
#for i in range(len(tcat)):
    
    
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
    
    #if name not in ['LBABOOJ143152.30+323216.7']: continue
    
    
    if os.path.exists(f'cutouts/specind_{name}.png') and (not clobber):
        continue
    
    if name not in ['LBABOOJ143653.04+341658.5']: continue
    
    
    tgcati = tgcat[tgcat['Source_id']==sid]
    
    c = SkyCoord(ra,dec,unit='deg')
    
    
    sep = c.separation(call)
    snear = (sep < 0.1*u.deg) & (sep > 0*u.deg)
    tcatnear = tcat[snear]
    
    dsep = c.separation(calldeep)
    sdmin = dsep.argmin()
    sdeepnear = (dsep < 0.1*u.deg) 
    tdeepcatn = tdeepcat[sdmin]
    tdeepcatnear = tdeepcat[sdeepnear]
    
    
    didsep = c.separation(calldeepid)
    sdmin = didsep.argmin()
    sdeepidnear = (didsep < 0.1*u.deg)   # 0.05 deg = 180 arcsec
    tdeepidcatn = tdeepidcat[sdmin]
    tdeepidcatnear = tdeepidcat[sdeepidnear]
    
    
    #imsize = 0.1
    #if not (os.path.isfile(f'cutouts/{name}_I.fits') or os.path.isfile(f'cutouts/{name}_legacy_r.fits')):
        #try:
            #cc.get_NDWFS_cutout(f'cutouts/{name}_I.fits', ra, dec, imsize, band="I")
            #tndwfs = pf.open(f'cutouts/{name}_I.fits')
            
        #except:
            #cc.get_legacy(f'cutouts/{name}_legacy_r.fits',ra,dec)
            
    #if os.path.isfile(f'cutouts/{name}_I.fits'):
        #tndwfs = pf.open(f'cutouts/{name}_I.fits')
        #otitle = 'NDWFS I'
    #else:
        #tndwfs = pf.open(f'cutouts/{name}_legacy_r.fits')
        #otitle = 'LEGACY r'
            
    

    ##cc.get_SDWFS_cutout_local('cutouts/{name}_I2.fits', ra, dec, imsize)
    ##tsdwfs = pf.open('cutouts/{name}_I2.fits')
    #try: 
        #tsdwfs = mm.extract_subim('/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/Bootes_optical/SDWFS/I2_bootes.v32.fits', ra, dec, 0.1)
        #irtitle = 'SDWFS I1'
    #except:
        #outwise = cc.get_wise(ra,dec,1)
        #tsdwfs = mm.extract_subim(outwise, ra, dec, 0.1)
        #irtitle = 'WISE 1'
        

    
    t = mm.extract_subim(dpath+imfits, ra,dec, 0.15)
    tdeep = mm.extract_subim(deepim, ra,dec, 0.15)
    tdeeph = pf.open(deepim)[0].header
    th = pf.open(dpath+imfits)[0].header
    
    lbadat = t[0].data
    
    h = t[0].header
    hdeep = tdeep[0].header
    
    hdeep['RESTFRQ'] = tdeeph['RESTFRQ']
    hdeep['BUNIT'] = tdeeph['BUNIT']
    h['RESTFRQ'] = th['RESTFRQ']
    h['BUNIT'] = th['BUNIT']
    #h = pf.getheader(dpath+imfits)
    #h['WCSAXES'] = 2
    #h['NAXIS'] = 2
    #for tkey in ['CRPIX3','CDELT3','CTYPE3','CRVAL3','CUNIT3','NAXIS3','CRPIX4','CDELT4','CTYPE4','CRVAL4','CUNIT4','NAXIS4']:
        #if tkey in h.keys():
            #h.remove(tkey)
    #hdeep = pf.getheader(deepim)
    
    pf.writeto(name+'lba_deep.fits', t[0].data, header=h, overwrite=True)
    pf.writeto(name+'hba_deep.fits', tdeep[0].data, header=hdeep, overwrite=True)
    
    '''
    convolve from
    0.001666666666666667 deg 
    0.00041666666666667 deg / pixel
    to
    0.004166666666666667 
    '''
    
    gf = 2.*np.sqrt(2.*np.log(2.))
    
    beamdeep = hdeep['BMAJ'] /gf #/ abs(hdeep['CDELT2'])
    beam = h['BMAJ'] /gf #/ abs(h['CDELT2'])
    conv = np.sqrt(beam**2 - beamdeep**2.)
    conv_pix_maj = conv / abs(hdeep['CDELT2'])
    
    beamdeepm = hdeep['BMIN']/gf #/ abs(hdeep['CDELT2'])
    beamm = h['BMIN']/gf #/ abs(h['CDELT2'])
    convm = np.sqrt(beamm**2 - beamdeepm**2.)
    conv_pix_min = convm / abs(hdeep['CDELT2'])
    
    conv_pix_pa = hdeep['BPA']
    kernel = Gaussian2DKernel(x_stddev=conv_pix_maj, y_stddev=conv_pix_min, theta=conv_pix_pa)
    
    
    
    #tdeep[0].data = tdeep[0].data * 1e3
    deepdat = tdeep[0].data
    scaleF = (beamdeep*beamdeepm)/(beam*beamm)
    
    deepdatC = convolve(deepdat, kernel)
    deepdatC = deepdatC/scaleF
    hdeep['BMAJ'] = h['BMAJ']
    hdeep['BMIN'] = h['BMIN']
    hdeep['BPA'] = h['BPA']
    pf.writeto(name+'hba_deep_convolved.fits', deepdatC, header=hdeep, overwrite=True)
    # concolvce deep to 15" resolution
    #deepdatC = 
    tdeep[0].data = deepdatC
    
    deepdatCr, footprint = reproject_interp(tdeep[0], h)
    h['RESTFRQ'] = tdeeph['RESTFRQ']

    
    tdeep[0].data = deepdatCr
    
    pf.writeto(name+'hba_deep_convolved_reproj.fits', deepdatCr, header=h, overwrite=True)
    
    datdeep = tdeep[0].data
    rmsdeep = tdeepcatn['Isl_rms']*1e3
    dvmin = -1.*rmsdeep
    dvmax = 25.*rmsdeep
    dvmid = 10.*rmsdeep
    #dvmax = tdeepcatn['Peak_flux']*1e3
    #dvmid = np.min((10.*rmsdeep, 0.5*dvmax))
    
    t[0].data = t[0].data * 1e3
    rms = tcati['Isl_rms'] *1e3
    vmin = -1.*rms
    vmax = 25.*rms
    vmid = 10.*rms
    #vmax = tcati['Peak_flux']*1e3
    #vmid = np.min((10.*rms, 0.5*vmax))

    
    #nvmin, nvmax, nvmid = get_scale(tndwfs)
    #svmin, svmax, svmid = get_scale(tsdwfs)
    
    
    
    from astropy.wcs import WCS
    import matplotlib.pyplot as plt

    pp.paper_single(TW=19.91,AR=0.35)
    ax1 = plt.subplot(1,2,1, projection=WCS(h))
    ax1.imshow(t[0].data, origin='lower', vmin=-100., vmax=2000.)
    #ax1.coords['ra'].set_axislabel('Right Ascension')
    #ax1.coords['dec'].set_axislabel('Declination')
    ax1.set_title('original')

    ax2 = plt.subplot(1,2,2, projection=WCS(hdeep))
    ax2.imshow(deepdatCr, origin='lower', vmin=-2.e-4, vmax=5.e-4)
    ax2.set_title('convolved-regridded')
    
    _, rmsdeep = std_sigclip(deepdatCr)
    _, rmslba = std_sigclip(lbadat)
    specind = np.log10(lbadat/deepdatCr) / np.log10(54./144.)
    specind_mask = (lbadat > 3*rmslba) & (deepdatCr > 3*rmsdeep )
    specind[~specind_mask] = np.nan
    
    pp.paper_single(TW=19.91,AR=0.35)
    ax1 = plt.subplot(1,3,1, projection=WCS(h))
    ax1.imshow(lbadat, origin='lower', vmin=-1*rmslba, vmax=25.*rmslba)
    ax2 = plt.subplot(1,3,2, projection=WCS(h))
    ax2.imshow(deepdatCr, origin='lower', vmin=-2*rmsdeep, vmax=25.*rmsdeep)
    ax3 = plt.subplot(1,3,3, projection=WCS(h))
    ax3.imshow(specind, origin='lower')#, vmin=-2.e-4, vmax=5.e-4)


    sys.exit()

    f = plt.figure()
    
    left = 0.05
    bottom = 0.05
    top = 0.95
    right = 0.95
    #dX = 0.25
    #dY = 0.3
    wspace = 0.05
    hspace = 0.05
    nX = 3
    nY = 2
    dX = (right-left-(nX-1)*wspace)/ nX
    dY = (top-bottom-(nY-1)*hspace)/ nY
    
    
    ## bottom-left lofar
    ax1 = ap.FITSFigure(t[0],figure=f, subplot=[left,bottom,dX,dY])
    ax1.set_title('LOFAR60')
    ax1.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=vmin, vmid=vmid, vmax=vmax, stretch='arcsinh')

    ax1.show_ellipses(tcatnear['RA'], tcatnear['DEC'], tcatnear['Maj'], tcatnear['Min'], tcatnear['PA']+90, edgecolor='C2', facecolor='none') #green
    ax1.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none') #orange
    ax1.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none') # blue
    #ax1.show_ellipses(tdeepcatnear['RA'], tdeepcatnear['DEC'], tdeepcatnear['Maj'], tdeepcatnear['Min'], tcati['PA']+90, edgecolor='C3', facecolor='none')

    #ax1.add_colorbar()
    #ax1.colorbar.set_axis_label_text('Intensity (mJy/bm)') 

    ax1.recenter(ra,dec,width=0.125,height=0.125)
    #pp.fig_save_many(f, 'cutouts/{name}')
    
    
    ax1.axis_labels.hide()
    ax1.tick_labels.hide()
    
    
    ## bottom-mid lofar zoom
    ax2 = ap.FITSFigure(t[0],figure=f, subplot=[left+wspace+dX,bottom,dX,dY])
    ax2.set_title('LOFAR60 Zoom')
    
    ax2.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=vmin, vmid=vmid, vmax=vmax, stretch='arcsinh')

    ax2.show_ellipses(tcatnear['RA'], tcatnear['DEC'], tcatnear['Maj'], tcatnear['Min'], tcatnear['PA']+90, edgecolor='C2', facecolor='none')
    ax2.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
    ax2.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')
    #ax2.show_ellipses(tdeepcatnear['RA'], tdeepcatnear['DEC'], tdeepcatnear['Maj'], tdeepcatnear['Min'], tcati['PA']+90, edgecolor='C3', facecolor='none')

    #ax2.add_colorbar()
    #ax2.colorbar.set_axis_label_text('Intensity (mJy/bm)') 
    
    ax2.recenter(ra,dec,width=0.05,height=0.05)
    
    #ax2.hide_yaxis_label()
    #pp.fig_save_many(f, 'cutouts/{name}_z1')
    
    ax2.axis_labels.hide()
    ax2.tick_labels.hide()
    
    
    ## top-left hba
    ax3 = ap.FITSFigure(tdeep[0],figure=f, subplot=[left,bottom+hspace+dY,dX,dY])
    ax3.set_title('LOFAR150')
    ax3.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=dvmin, vmid=dvmid, vmax=dvmax, stretch='arcsinh')

    #ax3.show_contour(t[0], levels=rms*np.array((3.,5.,10.,30.,50.,100.)),colors='C3')
    #ax3.show_ellipses(tcatnear['RA'], tcatnear['DEC'], tcatnear['Maj'], tcatnear['Min'], tcatnear['PA']+90, edgecolor='C2', facecolor='none')
    #ax3.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
    #ax3.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')
    ax3.show_ellipses(tdeepcatnear['RA'], tdeepcatnear['DEC'], tdeepcatnear['Maj'], tdeepcatnear['Min'], tdeepcatnear['PA']+90, edgecolor='C3', facecolor='none')
    ax3.show_ellipses(tdeepidcatnear['RA'], tdeepidcatnear['DEC'], tdeepidcatnear['Maj'], tdeepidcatnear['Min'], tdeepidcatnear['PA']+90, edgecolor='C4', facecolor='none')
    ax3.show_ellipses(tdeepidcatnear['RA'], tdeepidcatnear['DEC'], tdeepidcatnear['LGZ_Size'], tdeepidcatnear['LGZ_Width'], tdeepidcatnear['LGZ_PA']+90, edgecolor='C4', facecolor='none')

    #ax3.add_colorbar()
    #ax3.colorbar.set_axis_label_text('Intensity (mJy/bm)') 

    ax3.recenter(ra,dec,width=0.125,height=0.125)
    
    
    ax3.axis_labels.hide()
    ax3.tick_labels.hide()
    
    #ax3.hide_xaxis_label()
    
    ## top-mid hba zoom
    ax4 = ap.FITSFigure(tdeep[0],figure=f, subplot=[left+wspace+dX,bottom+hspace+dY,dX,dY])
    ax4.set_title('LOFAR150 Zoom')
    ax4.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=dvmin, vmid=dvmid, vmax=dvmax, stretch='arcsinh')

    #ax4.show_contour(t[0], levels=rms*np.array((3.,5.,10.,30.,50.,100.)),colors='C3')
    #ax4.show_ellipses(tcatnear['RA'], tcatnear['DEC'], tcatnear['Maj'], tcatnear['Min'], tcatnear['PA']+90, edgecolor='C2', facecolor='none')
    #ax4.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
    #ax4.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')
    ax4.show_ellipses(tdeepcatnear['RA'], tdeepcatnear['DEC'], tdeepcatnear['Maj'], tdeepcatnear['Min'], tdeepcatnear['PA']+90, edgecolor='C3', facecolor='none')
    ax4.show_ellipses(tdeepidcatnear['RA'], tdeepidcatnear['DEC'], tdeepidcatnear['LGZ_Size'], tdeepidcatnear['LGZ_Width'], tdeepidcatnear['LGZ_PA']+90, edgecolor='C4', facecolor='none')

    #ax4.add_colorbar()
    #ax4.colorbar.set_axis_label_text('Intensity (mJy/bm)') 
    ax4.recenter(ra,dec,width=0.05,height=0.05)
    
    
    ax4.axis_labels.hide()
    ax4.tick_labels.hide()
    

    pp.fig_save_many(f, 'cutouts/{name}_HBA_LBA_specind')
    
    
    
    #plt.close('all')



    #sys.exit()

