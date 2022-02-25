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

import mocpy

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

Ncat = len(tcat)

if 'HBA_Source_Name' in tcat.colnames:
    restart = True
else:
    tcat.add_column(Column(name='HBA_Source_Name',data=np.zeros(Ncat, dtype='S100')))
    tcat.add_column(Column(name='HBA_deep_cover',data=np.zeros(Ncat, dtype=bool)))
    restart = False


mergelist = 'done_list.txt'
if os.path.isfile(mergelist):
    
    with open(dpath+mergelist, 'r') as f:
        l = f.readlines()
    donelist = np.array([ll.strip() for ll in l])
else:
    donelist = []


all_merge = []

for i in range(len(tcat)):
    
    
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

    if name  in donelist:
        continue
    
    if name in all_merge:
        # already added to a merge...
        continue
    
    
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
    
    
    imsize = 0.1
    try:
        cc.get_NDWFS_cutout('cutouts/{name}_I.fits'.format(name=name), ra, dec, imsize, band="I")
        
        tndwfs = pf.open('cutouts/{name}_I.fits'.format(name=name))
        otitle = 'NDWFS I'
    except:
        cc.get_legacy('cutouts/{name}_legacy_r.fits'.format(name=name),ra,dec)
        
        tndwfs = pf.open('cutouts/{name}_legacy_r.fits'.format(name=name))
        otitle = 'LEGACY r'
    

    #cc.get_SDWFS_cutout_local('cutouts/{name}_I2.fits'.format(name=name), ra, dec, imsize)
    #tsdwfs = pf.open('cutouts/{name}_I2.fits'.format(name=name))
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

    ax1.show_ellipses(tcatn['RA'], tcatn['DEC'], tcatn['Maj'], tcatn['Min'], tcatn['PA']+90, edgecolor='C2', facecolor='none') #green
    ax1.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none') #orange
    ax1.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none') # blue

    #ax1.add_colorbar()
    #ax1.colorbar.set_axis_label_text('Intensity (mJy/bm)') 

    ax1.recenter(ra,dec,width=0.25,height=0.25)
    #pp.fig_save_many(f, 'cutouts/{name}'.format(name=name))
    
    
    ax1.axis_labels.hide()
    ax1.tick_labels.hide()
    
    
    ## bottom-mid lofar zoom
    ax2 = ap.FITSFigure(t[0],figure=f, subplot=[left+wspace+dX,bottom,dX,dY])
    ax2.set_title('LOFAR60 Zoom')
    
    ax2.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=vmin, vmid=vmid, vmax=vmax, stretch='arcsinh')

    ax2.show_ellipses(tcatn['RA'], tcatn['DEC'], tcatn['Maj'], tcatn['Min'], tcatn['PA']+90, edgecolor='C2', facecolor='none')
    ax2.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
    ax2.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')

    #ax2.add_colorbar()
    #ax2.colorbar.set_axis_label_text('Intensity (mJy/bm)') 
    
    ax2.recenter(ra,dec,width=0.1,height=0.1)
    
    #ax2.hide_yaxis_label()
    #pp.fig_save_many(f, 'cutouts/{name}_z1'.format(name=name))
    
    ax2.axis_labels.hide()
    ax2.tick_labels.hide()
    
    
    ## top-left hba
    ax3 = ap.FITSFigure(tdeep[0],figure=f, subplot=[left,bottom+hspace+dY,dX,dY])
    ax3.set_title('LOFAR150 deep')
    ax3.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=dvmin, vmid=dvmid, vmax=dvmax, stretch='arcsinh')

    ax3.show_contour(t[0], levels=rms*np.array((3.,5.,10.,30.,50.,100.)),colors='C3')
    ax3.show_ellipses(tcatn['RA'], tcatn['DEC'], tcatn['Maj'], tcatn['Min'], tcatn['PA']+90, edgecolor='C2', facecolor='none')
    ax3.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
    ax3.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')

    #ax3.add_colorbar()
    #ax3.colorbar.set_axis_label_text('Intensity (mJy/bm)') 

    ax3.recenter(ra,dec,width=0.25,height=0.25)
    
    
    ax3.axis_labels.hide()
    ax3.tick_labels.hide()
    
    #ax3.hide_xaxis_label()
    #pp.fig_save_many(f, 'cutouts/{name}_HBA'.format(name=name))
    
    ## top-mid hba zoom
    ax4 = ap.FITSFigure(tdeep[0],figure=f, subplot=[left+wspace+dX,bottom+hspace+dY,dX,dY])
    ax4.set_title('LOFAR150 deep Zoom')
    ax4.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=dvmin, vmid=dvmid, vmax=dvmax, stretch='arcsinh')

    ax4.show_contour(t[0], levels=rms*np.array((3.,5.,10.,30.,50.,100.)),colors='C3')
    ax4.show_ellipses(tcatn['RA'], tcatn['DEC'], tcatn['Maj'], tcatn['Min'], tcatn['PA']+90, edgecolor='C2', facecolor='none')
    ax4.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
    ax4.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')

    #ax4.add_colorbar()
    #ax4.colorbar.set_axis_label_text('Intensity (mJy/bm)') 
    ax4.recenter(ra,dec,width=0.1,height=0.1)
    
    
    ax4.axis_labels.hide()
    ax4.tick_labels.hide()
    
    #
    
    ## top-right ndwfs
    try:
        ax5 = ap.FITSFigure(tndwfs[0],figure=f, subplot=[left+2*(wspace+dX),bottom+hspace+dY,dX,dY])
        ax5.set_title(otitle)
        ax5.show_colorscale(cmap=plt.cm.cubehelix_r, vmin=nvmin, vmid=nvmid, vmax=nvmax, stretch='arcsinh')

        ax5.show_contour(t[0], levels=rms*np.array((3.,5.,10.,30.,50.,100.)),colors='C3')
        ax5.show_ellipses(tcatn['RA'], tcatn['DEC'], tcatn['Maj'], tcatn['Min'], tcatn['PA']+90, edgecolor='C2', facecolor='none')
        ax5.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
        ax5.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')

        #ax5.add_colorbar()
        #ax5.colorbar.set_axis_label_text('Intensity ') 
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

    ax6.show_contour(t[0], levels=rms*np.array((3.,5.,10.,30.,50.,100.)),colors='C3')
    ax6.show_ellipses(tcatn['RA'], tcatn['DEC'], tcatn['Maj'], tcatn['Min'], tcatn['PA']+90, edgecolor='C2', facecolor='none')
    ax6.show_ellipses(tgcati['RA'], tgcati['DEC'], tgcati['Maj'], tgcati['Min'], tgcati['PA']+90, edgecolor='C1', facecolor='none')
    ax6.show_ellipses(tcati['RA'], tcati['DEC'], tcati['Maj'], tcati['Min'], tcati['PA']+90, edgecolor='C0', facecolor='none')

    #ax6.add_colorbar()
    #ax6.colorbar.set_axis_label_text('Intensity (mJy/bm)') 
    try:
        ax6.recenter(ra,dec,width=0.05,height=0.05)
    except:
        print('ax6 error')
    
    
    ax6.axis_labels.hide()
    ax6.tick_labels.hide()
    
    def onclick(event):
        global ras, decs
        xp = event.xdata
        yp = event.ydata
        ax = event.inaxes
        ra,dec=ax.wcs.all_pix2world(xp,yp,1)
        ras.append(float(ra))
        decs.append(float(dec))
        return 
    
    f.canvas.mpl_connect('button_press_event', onclick)
    #pp.fig_save_many(f, 'cutouts/{name}'.format(name=name))
    
    
    plt.show()
    
    stop = False
    ras = []
    decs = []
    while not stop:
        print('Source ',name)
        print('Click on all sources belonging to this source')
        print('(n)ext')
        print('(p)rint list')
        print('(s)kip to next')
        print('(q)uit')
        s = input('input: ')
        if s =='q':
            stop = True
            plt.close('all')
            sys.exit()
        elif s =='n':
            stop = True
        elif s == 'p':
            print(ras)
            print(decs)
        elif s == 's':
            print('skipping to next')
            continue
        else:
            print('not handled')
    
    print(ras)
    print(decs)
    
    names = []
    C = SkyCoord(ras,decs,unit='deg')
    for ci in C:
        sep = ci.separation(call)
        argmin = sep.argmin()
        mergename = tcat['Source_Name'][argmin]
        names.append(mergename)
        
    # add this source to the list if not clicked on
    if name not in names:
        names.append(name)
        
    for mergename in names:
        all_merge.append(mergename)
        
    print(names)
    
    with open(name+'_merge.txt','w') as f:
        for n in names:
            f.write(n+'\n')
    
    plt.close('all')



    #sys.exit()

