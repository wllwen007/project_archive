import astropy.io.fits as pf
import bdsf     

'''
dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'
infile = 'image_full_phase_m.NS_shift.int.facetRestored.blanked.fits'
detectin = 'image_full_phase_m.NS_shift.app.facetRestored.blanked.fits'
#infile = 'image_full_phase_m.NS.int.restored.blanked.fits'
#detectin = 'image_full_phase_m.NS.app.restored.blanked.fits'
catprefix = 'bootes_deep_lba_hbashift'

#dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/test_apply_selfcal_screen/test_flagged/image_dirs28/'
#infile = 'image_test_gainscreen_dirs28.int.restored.fits'
#detectin = 'image_test_gainscreen_dirs28.app.restored.fits'
#catprefix = 'bootes_deep_test_flagged_dirs28'


dpath = '/data2/wwilliams/surveys/bootes_selfcal/ddf/deep'
infile = 'image_test_gainscreen_dirs12_smallf_deep.int.restored.fits'
detectin = 'image_test_gainscreen_dirs12_smallf_deep.app.restored.fits'
#infile = 'image_full_phase_m.NS.int.restored.blanked.fits'
#detectin = 'image_full_phase_m.NS.app.restored.blanked.fits'
catprefix = 'bootes_deep_dirs12'
'''

import configparser
config = configparser.ConfigParser()
config.read('image.config')
dpath = config['SETUP']['dpath']

appimfits = config['pbcor']['appimfits']
intimfits = config['pbcor']['intimfits']
pbimfits = config['pbcor']['pbimfits']

detectin = config['pbcor']['bappimfits']
infile = config['pbcor']['bintimfits']

dobands = config['sources']['dobands']

catprefix = config['sources']['catprefix']



head = pf.getheader(infile)
if 1:

    restfrq=head['RESTFRQ']

    img = bdsf.process_image(infile, detection_image=detectin, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=True, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None, frequency=restfrq)    
    img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True',clobber=True)
    img.export_image(outfile=catprefix +'.rms.fits',img_type='rms',img_format='fits',clobber=True)    
    img.export_image(outfile=catprefix +'.mean.fits',img_type='mean',img_format='fits',clobber=True)
    img.export_image(outfile=catprefix +'.resid.fits',img_type='gaus_resid',img_format='fits',clobber=True)
    img.export_image(outfile=catprefix +'.pybdsfmask.fits',img_type='island_mask',img_format='fits',clobber=True)
    img.write_catalog(outfile=catprefix +'.cat.reg',catalog_type='srl',format='ds9',correct_proj='True',clobber=True)


## for the spectral images
if dobands == 'True':
    print ('Error: not fully implemented')
    sys.exit(1)
    specpath = config['sources']['specpath']
    for i in range(2):
        dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/specim'
        dpath = specpath
        infile = 'image_full_ampphase_di_m.NS_Band{i}_shift.int.facetRestored.fits'.format(i=i)
        detectin = 'image_full_ampphase_di_m.NS_Band{i}_shift.app.facetRestored.fits'.format(i=i)
        catprefix = 'bootes_deep_lba_band{i}'.format(i=i)
        head = pf.getheader(infile)

        restfrq=head['RESTFRQ']

        if 0:
            img = bdsf.process_image(infile, detection_image=detectin, thresh_isl=4.0, thresh_pix=5.0, rms_box=(150,15), rms_map=True, mean_map='zero', ini_method='intensity', adaptive_rms_box=True, adaptive_thresh=150, rms_box_bright=(60,15), group_by_isl=False, group_tol=10.0, output_opts=True, output_all=True, atrous_do=True, atrous_jmax=4, flagging_opts=True, flag_maxsize_fwhm=0.5,advanced_opts=True, blank_limit=None, frequency=restfrq)    
            img.write_catalog(outfile=catprefix +'.cat.fits',catalog_type='srl',format='fits',correct_proj='True')
            img.export_image(outfile=catprefix +'.rms.fits',img_type='rms',img_format='fits',clobber=True)
            img.export_image(outfile=catprefix +'.resid.fits',img_type='gaus_resid',img_format='fits',clobber=True)
            img.export_image(outfile=catprefix +'.pybdsfmask.fits',img_type='island_mask',img_format='fits',clobber=True)
            img.write_catalog(outfile=catprefix +'.cat.reg',catalog_type='srl',format='ds9',correct_proj='True')
