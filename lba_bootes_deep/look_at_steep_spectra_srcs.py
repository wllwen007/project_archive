import os
import utils.plotting as pp
from utils.make_subim import extract_subim
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table, Column, MaskedColumn, join
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyregion
from astropy.wcs import WCS
from astropy.io import fits

#clobber = True
clobber = 1
stilts = '/net/beerze/data2/wwilliams/software/bin/stilts '

#lofar_cat_name = 'bootes_deep_lba.cat.fits'
lofar_cat_name = 'bootes_deep_lba_hbashift.cat.fits'

# TODO
#nudeep = 148.72 ??  - in header
nudeep = 144.

facet_regfile= 'facets_fin_finer.reg'

## match to deep opt

lofar_code_name = 'lba_match_deephba_codes.fits'

m_lofar_cat_name = 'bootes_deep_lba.cat.match_{n}.fits'.format(n='deephba')
m_lofar_cat_name = 'bootes_deep_lba_hbashift.cat.match_{n}.fits'.format(n='deephba')
cat_name = '/net/beerze/data2/wwilliams/projects/lofar_surveys/deep/Bootes/image_full_ampphase_di_m.NS_shift.blanked.scaled.cat.fits'
fmatch_cat = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/Bootes/bootes_final_cross_match_catalogue-v1.0-agnclass.fits'
# this is these 2 matched in topcat
#fmatch_cat = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/Bootes/bootes_final_cross_match_catalogue-v1.0.fits'
#fagn_class_cat = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/SEDfits_v1.0/AGNclasses_Bootes_v1.fits'
m_lofar_deephba_cat_name = '/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/Bootes/radio_bootes_final_cross_match_catalogue-v1.0.fits'

ocat = Table.read('/net/beerze//data2/wwilliams/projects/lofar_surveys/deep/science_ready_catalogs/Bootes_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits')


C = SkyCoord( '14h35m16.07s +33d27m05.6s')
