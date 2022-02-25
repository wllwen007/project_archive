import os
import numpy as np
from astropy.table import Table, Column
from utils.catalogues import name_from_coords

clobber = False



dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7/fin_im/'

cat = 'bootes_deep_lba.cat.fits'
tcat = Table.read(dpath+cat)

artlist = 'artefact_list.txt'
with open(dpath+artlist, 'r') as f:
    l = f.readlines()
artefacts = np.array([ll.strip() for ll in l])

if (not 'Artefact' in tcat.colnames) or clobber:
    
    art = np.zeros(len(tcat), dtype=bool)
    for artefact in artefacts:
        if artefact not in tcat['Source_Name']:
            print(artefact, 'not in list')
            continue
        ind = np.where(artefact==tcat['Source_Name'])
        art[ind] = 1
    
    if 'Artefact' in tcat.colnames:
        tcat['Artefact'] = art
    else:
        tcat.add_column(Column(name='Artefact',data=art))
    tcat.write(dpath+cat, overwrite=True)
