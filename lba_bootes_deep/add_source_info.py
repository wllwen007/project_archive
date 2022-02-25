import os
import numpy as np
from astropy.table import Table, Column, join
from utils.catalogues import name_from_coords

clobber = True



dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'

cat = 'bootes_deep_lba_hbashift.cat.fits'
tcat = Table.read(dpath+cat)



gcat = 'bootes_deep_lba_hbashift.gcat.fits'
tgcat = Table.read(dpath+gcat)

if clobber:
    if 'Ngaus' in tcat.colnames:
        tcat.remove_column('Ngaus')
        
        

if (not 'Ngaus' in tcat.colnames):
    
    ngaus = np.ones(len(tcat), dtype=int)
    for i,tcati in enumerate(tcat):
        if tcati['S_Code'] == 'S':
            continue
        
        ngaus[i] = np.sum(tcati['Source_id'] == tgcat['Source_id'])
        
        
    if not clobber:
        tcat.add_column(Column(name='Ngaus', data=ngaus))
    else:
        tcat['Ngaus']= ngaus

    tcat.write(dpath+cat, overwrite=True)

