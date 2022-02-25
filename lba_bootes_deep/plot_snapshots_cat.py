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


clobber = False
#clobber=True

dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'
imfits = 'image_full_phase_m.NS_shift.int.facetRestored.blanked.fits'

cat = 'bootes_deep_lba_hbashift.cat.fits'
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
    
    
    if os.path.exists(f'cutouts/{name}.png') and (not clobber):
        continue
    
    name = name.replace('LBABOO ','')
    
    #if name not in ['LBABOOJ143516.07+332705.6', 'LBABOOJ142250.33+325206.5', 
                    #'LBABOOJ142018.22+341447.6',
                    #'LBABOOJ144338.05+343034.8',
                    #'LBABOOJ142625.53+350326.6']: continue
    
    #if name not in ['LBABOOJ142440.51+323605.9',
#'LBABOOJ142459.44+352842.7',
#'LBABOOJ142505.39+334652.5',
#'LBABOOJ142526.50+344922.1',
#'LBABOOJ142531.13+353343.1',
#'LBABOOJ142531.79+353340.1',
#'LBABOOJ142616.46+342121.1',
#'LBABOOJ142639.46+350242.3',
#'LBABOOJ142649.68+332931.8',
#'LBABOOJ142702.71+333339.4',
#'LBABOOJ142703.01+333435.2',
#'LBABOOJ142706.50+333512.4',
#'LBABOOJ142707.12+333428.9',
#'LBABOOJ142755.95+332141.7',
#'LBABOOJ142756.20+332141.5',
#'LBABOOJ142759.95+324739.6',
#'LBABOOJ142806.62+325935.3',
#'LBABOOJ142807.39+325932.0',
#'LBABOOJ142838.55+342341.1',
#'LBABOOJ142842.26+342447.3',
#'LBABOOJ142842.34+354325.1',
#'LBABOOJ142842.68+354324.5',
#'LBABOOJ142845.33+323837.4',
#'LBABOOJ142848.10+323940.6',
#'LBABOOJ142848.80+323822.2',
#'LBABOOJ142902.81+353823.1',
#'LBABOOJ142919.13+352813.0',
#'LBABOOJ142919.22+324633.0',
#'LBABOOJ142922.11+324626.5',
#'LBABOOJ142927.98+352956.2',
#'LBABOOJ142949.23+351009.6',
#'LBABOOJ142949.52+353440.3',
#'LBABOOJ142954.12+343517.5',
#'LBABOOJ142956.47+325518.4',
#'LBABOOJ143008.02+325833.5',
#'LBABOOJ143008.24+331537.8',
#'LBABOOJ143008.89+325840.0',
#'LBABOOJ143012.28+331437.9',
#'LBABOOJ143012.82+331433.1',
#'LBABOOJ143014.55+345918.1',
#'LBABOOJ143034.29+324321.6',
#'LBABOOJ143039.15+352351.9',
#'LBABOOJ143052.56+335559.9',
#'LBABOOJ143053.48+330031.0',
#'LBABOOJ143057.83+334501.3',
#'LBABOOJ143106.60+334610.0',
#'LBABOOJ143108.07+352557.5',
#'LBABOOJ143108.11+352531.6',
#'LBABOOJ143114.42+323225.5',
#'LBABOOJ143115.19+334615.8',
#'LBABOOJ143120.83+325236.7',
#'LBABOOJ143121.19+353511.5',
#'LBABOOJ143123.18+353551.6',
#'LBABOOJ143128.44+342726.3',
#'LBABOOJ143134.71+351505.9',
#'LBABOOJ143136.20+351442.0',
#'LBABOOJ143143.19+353551.8',
#'LBABOOJ143203.68+330744.4',
#'LBABOOJ143217.30+351421.5',
#'LBABOOJ143223.38+332909.7',
#'LBABOOJ143227.46+354607.2',
#'LBABOOJ143233.96+354535.0',
#'LBABOOJ143235.24+354346.9',
#'LBABOOJ143237.32+354326.1',
#'LBABOOJ143239.68+323421.2',
#'LBABOOJ143242.90+322047.2',
#'LBABOOJ143252.09+322103.6',
#'LBABOOJ143312.78+351431.0',
#'LBABOOJ143329.43+343442.1',
#'LBABOOJ143342.21+341135.8',
#'LBABOOJ143404.64+331250.8',
#'LBABOOJ143426.89+353204.8',
#'LBABOOJ143429.72+342720.5',
#'LBABOOJ143430.81+354214.5',
#'LBABOOJ143431.46+350104.6',
#'LBABOOJ143433.02+342752.0',
#'LBABOOJ143433.17+352134.3',
#'LBABOOJ143443.18+333007.1',
#'LBABOOJ143443.26+332833.8',
#'LBABOOJ143444.26+332912.9',
#'LBABOOJ143445.94+332818.1',
#'LBABOOJ143447.03+332823.2',
#'LBABOOJ143448.35+354238.9',
#'LBABOOJ143450.84+354251.2',
#'LBABOOJ143602.97+334407.2',
#'LBABOOJ143603.34+334350.9',
#'LBABOOJ143627.51+352622.6',
#'LBABOOJ143657.97+355105.6',
#'LBABOOJ143728.63+350749.1',
#'LBABOOJ143737.93+342006.3',
#'LBABOOJ143827.92+354108.1']: continue
    
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
    
    
    imsize = 0.1
    if not (os.path.isfile(f'cutouts/{name}_I.fits') or os.path.isfile(f'cutouts/{name}_legacy_r.fits')):
        try:
            cc.get_NDWFS_cutout(f'cutouts/{name}_I.fits', ra, dec, imsize, band="I")
            tndwfs = pf.open(f'cutouts/{name}_I.fits')
            
        except:
            cc.get_legacy(f'cutouts/{name}_legacy_r.fits',ra,dec)
            
    if os.path.isfile(f'cutouts/{name}_I.fits'):
        tndwfs = pf.open(f'cutouts/{name}_I.fits')
        otitle = 'NDWFS I'
    else:
        tndwfs = pf.open(f'cutouts/{name}_legacy_r.fits')
        otitle = 'LEGACY r'
            
    

    #cc.get_SDWFS_cutout_local('cutouts/{name}_I2.fits', ra, dec, imsize)
    #tsdwfs = pf.open('cutouts/{name}_I2.fits')
    try: 
        tsdwfs = mm.extract_subim('/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/Bootes_optical/SDWFS/I2_bootes.v32.fits', ra, dec, 0.1)
        irtitle = 'SDWFS I1'
    except:
        outwise = cc.get_wise(ra,dec,1)
        tsdwfs = mm.extract_subim(outwise, ra, dec, 0.1)
        irtitle = 'WISE 1'
        

    
    t = mm.extract_subim(dpath+imfits, ra,dec, 0.25)
    tdeep = mm.extract_subim(deepim, ra,dec, 0.25)
    
    
    tdeep[0].data = tdeep[0].data * 1e3
    
    #datdeep = tdeep[0].data
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

    
    nvmin, nvmax, nvmid = get_scale(tndwfs)
    svmin, svmax, svmid = get_scale(tsdwfs)
    
    pp.paper_single(TW=19.91,AR=0.65)

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
    #pp.fig_save_many(f, 'cutouts/{name}_HBA')
    
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
    
    #
    
    ## top-right ndwfs
    try:
        ax5 = ap.FITSFigure(tndwfs[0],figure=f, subplot=[left+2*(wspace+dX),bottom+hspace+dY,dX,dY])
        ax5.set_title(otitle)
        ax5.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=nvmin, vmid=nvmid, vmax=nvmax, stretch='arcsinh')

        ax5.show_contour(t[0], levels=rms*np.array((3.,5.,10.,30.,50.,100.)),colors='C0')
        ax5.show_contour(tdeep[0], levels=rmsdeep*np.array((3.,5.,10.,30.,50.,100.)),colors='C3')
        #ax5.show_ellipses(tcatnear['RA'], tcatnear['DEC'], tcatnear['Maj'], tcatnear['Min'], tcatnear['PA']+90, edgecolor='C2', facecolor='none')
        #ax5.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
        #ax5.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')
        #ax5.show_ellipses(tdeepcatnear['RA'], tdeepcatnear['DEC'], tdeepcatnear['Maj'], tdeepcatnear['Min'], tcati['PA']+90, edgecolor='C3', facecolor='none')

        #ax5.add_colorbar()
        #ax5.colorbar.set_axis_label_text('Intensity ') 
        
        ax5.show_markers(tdeepidcatnear['ALPHA_J2000'], tdeepidcatnear['DELTA_J2000'], marker='+', edgecolor='C3', facecolor='C3',s=100)
        ax5.recenter(ra,dec,width=0.05,height=0.05)
        
        
        ax5.axis_labels.hide()
        ax5.tick_labels.hide()
        
    except:
        print ('ax5 error')
    #
    ## bottom-right sdwfs
    ax6 = ap.FITSFigure(tsdwfs[0],figure=f, subplot=[left+2*(wspace+dX),bottom,dX,dY])
    ax6.set_title(irtitle)
    ax6.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=svmin, vmid=svmid, vmax=svmax, stretch='arcsinh')

    ax6.show_contour(t[0], levels=rms*np.array((3.,5.,10.,30.,50.,100.)),colors='C0')
    ax6.show_contour(tdeep[0], levels=rmsdeep*np.array((3.,5.,10.,30.,50.,100.)),colors='C3')
    #ax6.show_ellipses(tcatnear['RA'], tcatnear['DEC'], tcatnear['Maj'], tcatnear['Min'], tcatnear['PA']+90, edgecolor='C2', facecolor='none')
    #ax6.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
    #ax6.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')
    #ax6.show_ellipses(tdeepcatnear['RA'], tdeepcatnear['DEC'], tdeepcatnear['Maj'], tdeepcatnear['Min'], tcati['PA']+90, edgecolor='C3', facecolor='none')

    #ax6.add_colorbar()
    #ax6.colorbar.set_axis_label_text('Intensity (mJy/bm)')
    ax6.show_markers(tdeepidcatnear['ALPHA_J2000'], tdeepidcatnear['DELTA_J2000'], marker='+', edgecolor='C3', facecolor='C3',s=100) 
    try:
        ax6.recenter(ra,dec,width=0.05,height=0.05)
    except:
        print('ax6 error')
    
    
    ax6.axis_labels.hide()
    ax6.tick_labels.hide()
    
    
    
    pp.fig_save_many(f, f'cutouts/{name}')
    
    
    
    #plt.close('all')



    #sys.exit()

