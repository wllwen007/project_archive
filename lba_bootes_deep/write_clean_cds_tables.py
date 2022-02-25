import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column



dpath = '/net/lofar2/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7-v2/fin_im/'
lofar_cat_name = 'bootes_deep_lba_hbashift.cat.fits'
mcat_name = 'bootes_deep_lba.cat.matched.v1.0.fits'

cat = Table.read(dpath+lofar_cat_name)
mcat = Table.read(dpath+mcat_name)
cat.sort('Total_flux')
mcat.sort('Total_flux')

#mcat.rename_column('','HBA_Deep_flag')

Sf = 1000.   # scale flux in cat to print

#iters = np.hstack((np.arange(-1,-10,-1),[999], np.arange(0,10,1)) )

#for i in iters:
    #if i == 999:
        #s = '\\vdots & '
        #print(s+'\\\\[+0.5em]')
        #continue
    #cc=cat[i]
    #sline = " {sname:20s} & ${RA:9.5f}$ & ${sRA:6.1f}$ & ${DEC:9.5f}$  & ${sDEC:6.1f}$ & ${Si:6.1f}$ & ${eSi:6.1f}$ & ${Sp:6.1f}$ & ${eSp:6.1f}$  ".format(sname=cc['Source_Name'].replace('BOOJ','BOO\\,J'), RA=cc['RA'], sRA=3600*cc['E_RA'], DEC=cc['DEC'], sDEC=3600*cc['E_DEC'], Si=Sf*cc['Total_flux'], eSi=Sf*cc['E_Total_flux'] , Sp=Sf*cc['Peak_flux'], eSp=Sf*cc['E_Peak_flux']  )
    
    
    #sline += ' & ${rms:6.2f}$  & ${Ng:3d}$ '.format(rms=Sf*cc['Isl_rms'], Ng=cc['Ngaus'])
    
    #if cc['DC_Maj'] != 0:
        #sline += '& ${a:6.1f}$ & ${sa:6.1f}$ & ${b:6.1f}$ & ${sb:6.1f}$ & ${pa:6.1f}$ & ${spa:6.1f}$'.format(a=3600*cc['DC_Maj'], sa=3600*cc['E_DC_Maj'], b=3600*cc['DC_Min'], sb=3600*cc['E_DC_Min'], pa=cc['DC_PA'], spa=cc['E_DC_PA'])
        
    #sline +='\\\\' 
    
    #print(sline)


keepcols = ['Source_Name', 'RA', 'E_RA', 'DEC', 'E_DEC', 'Total_flux', 'E_Total_flux', 'Peak_flux', 'E_Peak_flux', 'Isl_rms', 'Ngaus', 'DC_Maj', 'E_DC_Maj', 'DC_Min', 'E_DC_Min', 'DC_PA', 'E_DC_PA']
cat.keep_columns(keepcols)
cat = cat[keepcols]
mcat.rename_column('opt_flag','flag_in_optical')
keepcols = ['Source_Name', 'RA', 'E_RA', 'DEC', 'E_DEC', 'Total_flux',  'E_Total_flux', 'Peak_flux', 'E_Peak_flux', 'Isl_rms', 'Ngaus', 'Size', 'Width', 'Angle', 'flag_in_optical', 'HBA_Deep_Name']
mcat.keep_columns(keepcols)
mcat = mcat[keepcols]
mcat['flag_in_optical'] = 1* mcat['flag_in_optical']   
cat.write('cds_table2.fits',overwrite=True)
mcat.write('cds_table_matched.fits',overwrite=True)
