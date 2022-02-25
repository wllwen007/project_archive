import os
from astropy.table import Table, Column, join
from utils.catalogues import name_from_coords

clobber = True



dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'

cat = 'bootes_deep_lba_hbashift.cat.fits'
tcat = Table.read(dpath+cat)



gcat = 'bootes_deep_lba_hbashift.gcat.fits'
tgcat = Table.read(dpath+gcat)

if clobber:
    if 'Source_Name' in tcat.colnames:
        tcat.remove_column('Source_Name')
    if 'Source_Name' in tgcat.colnames:
        tgcat.remove_column('Source_Name')
    if 'Gaussian_Name' in tgcat.colnames:
        tgcat.remove_column('Gaussian_Name')

if (not 'Source_Name' in tcat.colnames) or clobber:
    
    names = name_from_coords(tcat['RA'], tcat['DEC'], prefix='LBABOOJ')
    
    if not clobber:
        tcat.add_column(Column(name='Source_Name', data=names))
    else:
        tcat['Source_Name']= names

    tcat.write(dpath+cat, overwrite=True)


tj = join(tgcat, tcat, keys='Source_id', join_type='left')


if (not 'Gaussian_Name' in tgcat.colnames) or clobber:
 
    names = name_from_coords(tgcat['RA'], tgcat['DEC'], prefix='LBABOOJ')
    
    if not clobber:
        tgcat.add_column(Column(name='Gaussian_Name', data=names))
    else:
        tgcat['Gaussian_Name']= names

    if not clobber:
        tgcat.add_column(Column(name='Source_Name', data=tj['Source_Name']))
    else:
        tgcat['Source_Name'] = tj['Source_Name']

    tgcat.write(dpath+gcat, overwrite=True)
