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
clobber=True

dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'
imfits = 'image_full_phase_m.NS_shift.int.facetRestored.blanked.fits'

cat = 'bootes_deep_lba.cat.matched.v0.1.fits'
tcat = Table.read(dpath+cat)
compcat = 'bootes_deep_lba.cat.matched.components.v0.1.fits'
tcompcat = Table.read(dpath+compcat)


#gcat = 'bootes_deep_lba_hbashift.gcat.fits'
scat = 'bootes_deep_lba_hbashift.cat.fits'
tscat = Table.read(dpath+scat)

deepcat= '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.blanked.scaled.cat.fits'
deepidcat= '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/bootes_final_cross_match_catalogue-v1.0.fits'
deepim = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits'
deepomask = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/data_release/bootes/image_full_ampphase_di_m.NS_shift.blanked.scaled.rms_spmask.fits'

tdeepcat = Table.read(deepcat)
tdeepidcat = Table.read(deepidcat)


tdeepidcat.add_column(Column(name='Size', data=np.zeros(len(tdeepidcat))))
tdeepidcat.add_column(Column(name='Width', data=np.zeros(len(tdeepidcat))))
tdeepidcat.add_column(Column(name='Angle', data=np.zeros(len(tdeepidcat))))

lgz = ~np.isnan(tdeepidcat['LGZ_Size'])
tdeepidcat['Size'][lgz] =  tdeepidcat['LGZ_Size'][lgz] /3600.
tdeepidcat['Size'][~lgz] = 2*tdeepidcat['DC_Maj'][~lgz]
tdeepidcat['Width'][lgz] =  tdeepidcat['LGZ_Width'][lgz] /3600.
tdeepidcat['Width'][~lgz] = 2*tdeepidcat['DC_Min'][~lgz]
tdeepidcat['Angle'][lgz] =  tdeepidcat['LGZ_PA'][lgz]
tdeepidcat['Angle'][~lgz] = tdeepidcat['DC_PA'][~lgz]

tdeepidcat['Size'] = np.sqrt(tdeepidcat['Size']**2. + (5./3600)**2.)   # in degrees
tdeepidcat['Width'] = np.sqrt(tdeepidcat['Width']**2. + (5./3600)**2.)

tcompcat['Maj'] = tcompcat['Maj']/3600.
tcompcat['Min'] = tcompcat['Min']/3600.

#lgz = ~np.isnan(tcat['Size'])
#tcat['Size'][~lgz] = 2*tcat['DC_Maj'][~lgz]
#tcat['Width'][~lgz] = 2*tcat['DC_Min'][~lgz]
#tcat['Angle'][~lgz] = tcat['DC_PA'][~lgz]
tcat['Size'] = np.sqrt((tcat['Size']/3600.)**2. + (15./3600)**2.)
tcat['Width'] = np.sqrt((tcat['Width']/3600.)**2. + (15./3600)**2.)

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
calls = SkyCoord(tscat['RA'],tscat['DEC'],unit='deg')
calldeep = SkyCoord(tdeepcat['RA'],tdeepcat['DEC'],unit='deg')
calldeepid = SkyCoord(tdeepidcat['RA'],tdeepidcat['DEC'],unit='deg')

tcat.sort('HBA_Deep_Sep',reverse=True)   
for i in range(len(tcat)):
#for i in range(len(tcat)):
    
    
    #if i!=959: continue
    
    tcati = tcat[i]
    
    # only make the m sources - check for deblending...
    #if tcati['S_Code'] != 'S':
        #continue
        
    if not tcati['opt_flag']: continue
    if len(tcati['HBA_Deep_Name'])==0: continue
    
    name = tcati['Source_Name']
    
    #if name not in ['LBABOOJ143009.90+352002.5']: continue
    #if name not in ['LBABOOJ142839.37+342355.7']: continue
    
    ra = tcati['RA']
    dec = tcati['DEC']
    sid = tcati['Source_id']
    size = tcati['Size']  # in deg
    size = max(size, 60./3600)
    
    
    
    if os.path.exists(f'mcutouts/{name}.png') and (not clobber):
        continue
    
    name = name.replace('LBABOO ','')
    
    
    #tgcati = tgcat[tgcat['Source_id']==sid]
    
    c = SkyCoord(ra,dec,unit='deg')
    
    
    sep = c.separation(call)
    snear = (sep < size*u.deg)# & (sep > 0*u.deg)
    tcatnear = tcat[snear]
    
    tcompcats = tcompcat[tcompcat['Source_Name']==name]
    
    sep = c.separation(calls)
    snear = (sep < size*u.deg)# & (sep > 0*u.deg)
    tscatnear = tscat[snear]
    
    dsep = c.separation(calldeep)
    sdmin = dsep.argmin()
    sdeepnear = (dsep < size*u.deg) 
    tdeepcatn = tdeepcat[sdmin]
    tdeepcatnear = tdeepcat[sdeepnear]
    
    tdeepcatmatched = tdeepidcat[tdeepidcat['Source_Name']==tcati['HBA_Deep_Name']]
    
    
    didsep = c.separation(calldeepid)
    sdmin = didsep.argmin()
    sdeepidnear = (didsep < size*u.deg)   # 0.05 deg = 180 arcsec
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
    #t = flatten(fits.open(deepomask),None,None,None,None,None)
    tdeepm = mm.extract_subim(deepomask, ra,dec, 0.25)
    
    
    tdeep[0].data = tdeep[0].data * 1e3
    
    mask = np.ones_like(tdeepm[0].data)
    mask[np.isfinite(tdeepm[0].data)] = np.nan
    tdeepm[0].data = mask
    #tdeepm[0].data[np.isfinite(tdeepm[0].data)] = 1
    #tdeepm[0].data = 1*np.isnan(tdeepm[0].data)
    
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
    
    pp.paper_single(TW=13.27,AR=0.85)

    f = plt.figure()
    
    left = 0.05
    bottom = 0.05
    top = 0.95
    right = 0.95
    #dX = 0.25
    #dY = 0.3
    wspace = 0.05
    hspace = 0.05
    nX = 2
    nY = 2
    dX = (right-left-(nX-1)*wspace)/ nX
    dY = (top-bottom-(nY-1)*hspace)/ nY
    
    
    ## bottom-left lofar
    ax1 = ap.FITSFigure(t[0],figure=f, subplot=[left,bottom,dX,dY])
    ax2 = ax1
    ax2.set_title('LOFAR60 ')
    
    ax2.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=vmin, vmid=vmid, vmax=vmax, stretch='arcsinh')

    ax2.show_ellipses(tcatnear['RA'], tcatnear['DEC'], tcatnear['Size'], tcatnear['Width'], tcatnear['Angle']+90, edgecolor='C0', facecolor='none',linestyle='dashed')
    ax2.show_ellipses(tscatnear['RA'], tscatnear['DEC'], tscatnear['Maj'], tscatnear['Min'], tscatnear['PA']+90, edgecolor='C2', facecolor='none',linestyle='dashed')
    
    ax2.show_ellipses(tcompcats['RA'], tcompcats['DEC'], tcompcats['Maj'], tcompcats['Min'], tcompcats['PA']+90, edgecolor='C0', facecolor='none',linestyle='dashed')
    ax2.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Size'], tcati['Width'], tcati['Angle']+90, edgecolor='C0', facecolor='none')
    
    ax2.show_markers(tcati['RA'], tcati['DEC'], marker='x', edgecolor='C0', facecolor='C0',s=200)
    ax2.show_markers(tdeepcatmatched['RA'], tdeepcatmatched['DEC'], marker='x', edgecolor='C4', facecolor='C4',s=200)
    ax2.show_markers(tcati['HBA_Deep_ALPHA_J2000'], tcati['HBA_Deep_DELTA_J2000'], marker='+', edgecolor='C3', facecolor='C3',s=200)

    #ax2.add_colorbar()
    #ax2.colorbar.set_axis_label_text('Intensity (mJy/bm)') 
    
    ax2.recenter(ra,dec,width=size,height=size)
    
    #ax2.hide_yaxis_label()
    #pp.fig_save_many(f, 'cutouts/{name}_z1')
    
    ax2.axis_labels.hide()
    ax2.tick_labels.hide()
    
    
    ## top-left hba
    ax3 = ap.FITSFigure(tdeep[0],figure=f, subplot=[left,bottom+hspace+dY,dX,dY])
    ax4 = ax3
    
    ## top-mid hba zoom
    ax4.set_title('LOFAR150 ')
    ax4.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=dvmin, vmid=dvmid, vmax=dvmax, stretch='arcsinh')

    ax4.show_ellipses(tdeepcatnear['RA'], tdeepcatnear['DEC'], tdeepcatnear['Maj'], tdeepcatnear['Min'], tdeepcatnear['PA']+90, edgecolor='C3', facecolor='none',linestyle='dashed')
    
    ax4.show_ellipses(tdeepidcatnear['RA'], tdeepidcatnear['DEC'], tdeepidcatnear['Size'], tdeepidcatnear['Width'], tdeepidcatnear['Angle']+90, edgecolor='C4', facecolor='none',linestyle='dashed')

    ax4.show_markers(tcati['RA'], tcati['DEC'], marker='x', edgecolor='C0', facecolor='C0',s=200)
    ax4.show_markers(tdeepcatmatched['RA'], tdeepcatmatched['DEC'], marker='x', edgecolor='C4', facecolor='C4',s=200)
    ax4.show_markers(tcati['HBA_Deep_ALPHA_J2000'], tcati['HBA_Deep_DELTA_J2000'], marker='+', edgecolor='C3', facecolor='C3',s=200)
    ax4.show_ellipses(tdeepcatmatched['RA'], tdeepcatmatched['DEC'], tdeepcatmatched['Size'], tdeepcatmatched['Width'], tdeepcatmatched['Angle']+90, edgecolor='C4', facecolor='none')
    #ax4.add_colorbar()
    #ax4.colorbar.set_axis_label_text('Intensity (mJy/bm)') 
    ax4.recenter(ra,dec,width=size,height=size)
    
    
    ax4.axis_labels.hide()
    ax4.tick_labels.hide()
    
    #
    
    ## top-right ndwfs
    try:
        ax5 = ap.FITSFigure(tndwfs[0],figure=f, subplot=[left+1*(wspace+dX),bottom+hspace+dY,dX,dY])
        ax5.set_title(otitle)
        ax5.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=nvmin, vmid=nvmid, vmax=nvmax, stretch='arcsinh')


        ax5.show_contour(tdeepm[0], levels=np.array((0,1)),colors='gray',filled=True,layer='mask',alpha=0.5)
        ax5.show_contour(t[0], levels=rms*np.array((3.,5.,10.,30.,50.,100.)),colors='C0')
        ax5.show_contour(tdeep[0], levels=rmsdeep*np.array((3.,5.,10.,30.,50.,100.)),colors='C3')
        #ax5.show_ellipses(tcatnear['RA'], tcatnear['DEC'], tcatnear['Maj'], tcatnear['Min'], tcatnear['PA']+90, edgecolor='C2', facecolor='none')
        #ax5.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
        #ax5.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')
        #ax5.show_ellipses(tdeepcatnear['RA'], tdeepcatnear['DEC'], tdeepcatnear['Maj'], tdeepcatnear['Min'], tcati['PA']+90, edgecolor='C3', facecolor='none')

        #ax5.add_colorbar()
        #ax5.colorbar.set_axis_label_text('Intensity ') 
        
        ax5.show_markers(tdeepidcatnear['ALPHA_J2000'], tdeepidcatnear['DELTA_J2000'], marker='+', edgecolor='C3', facecolor='C3',s=200)
        ax5.recenter(ra,dec,width=size,height=size)
        
        
        ax5.axis_labels.hide()
        ax5.tick_labels.hide()
        
    except:
        print ('ax5 error')
    #
    ## bottom-right sdwfs
    ax6 = ap.FITSFigure(tsdwfs[0],figure=f, subplot=[left+1*(wspace+dX),bottom,dX,dY])
    ax6.set_title(irtitle)
    ax6.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=svmin, vmid=svmid, vmax=svmax, stretch='arcsinh')

    ax6.show_contour(tdeepm[0], levels=np.array((0,1)),colors='gray',filled=True,layer='mask',alpha=0.5)
    ax6.show_contour(t[0], levels=rms*np.array((3.,5.,10.,30.,50.,100.)),colors='C0')
    ax6.show_contour(tdeep[0], levels=rmsdeep*np.array((3.,5.,10.,30.,50.,100.)),colors='C3')
    #ax6.show_ellipses(tcatnear['RA'], tcatnear['DEC'], tcatnear['Maj'], tcatnear['Min'], tcatnear['PA']+90, edgecolor='C2', facecolor='none')
    #ax6.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
    #ax6.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')
    #ax6.show_ellipses(tdeepcatnear['RA'], tdeepcatnear['DEC'], tdeepcatnear['Maj'], tdeepcatnear['Min'], tcati['PA']+90, edgecolor='C3', facecolor='none')

    #ax6.add_colorbar()
    #ax6.colorbar.set_axis_label_text('Intensity (mJy/bm)')
    ax6.show_markers(tdeepidcatnear['ALPHA_J2000'], tdeepidcatnear['DELTA_J2000'], marker='+', edgecolor='C3', facecolor='C3',s=200) 
    try:
        ax6.recenter(ra,dec,width=size,height=size)
    except:
        print('ax6 error')
    
    
    ax6.axis_labels.hide()
    ax6.tick_labels.hide()
    
    
    
    pp.fig_save_many(f, f'mcutouts/{name}')
    
    
    
    plt.close('all')



    #sys.exit()

