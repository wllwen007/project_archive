from astropy.table import Table      

#bootes_deep_lba_hbashift.cat.fits
mcatin = Table.read('bootes_deep_lba.cat.matched.v0.1.fits')


mcatin = mcatin[mcatin['opt_flag']==True]
keep_cols = [
 'Source_Name',
 'RA',
 'E_RA',
 'DEC',
 'E_DEC',
 'Total_flux',
 'E_Total_flux',
 'Peak_flux',
 'E_Peak_flux',
 'Maj',
 'E_Maj',
 'Min',
 'E_Min',
 'PA',
 'E_PA',
 'DC_Maj',
 'E_DC_Maj',
 'DC_Min',
 'E_DC_Min',
 'DC_PA',
 'E_DC_PA',
 'Isl_rms',
 'S_Code',
 'Ngaus',
 'Size',
 'Width',
 'Angle',
 'FLAG_MASKED',
 'confused',
 'Deblended_From',
 'HBA_Deep_Name',
 'HBA_Deep_Sep']





descriptions = {
     'Source_Name': "Source Name in format LBABOOJhhmmss.ss+ddmmss.s",
'RA': "Right ascension in J2000",
'E_RA': "Error on the right ascension",
'DEC': "Declination in J2000",
'E_DEC': "Error on the declination",
'Total_flux': "Integrated flux density",
'E_Total_flux': "Error on the integrated flux density",
'Peak_flux': "Peak brightnes (affected by smearing)",
'E_Peak_flux': "Error on the peak bringhtness",
'Maj': "Source size, major axis ",
'E_Maj': "Error on the major axis size",
'Min': "Source size, minor axis",
'E_Min': "Error on the minor axis",
'PA': "Position angle",
'E_PA': "Error on the position angle",
'DC_Maj': "Deconvolved source size, major axis",
'E_DC_Maj': "Error on the deconvolved major axis size",
'DC_Min': "Deconvolved source size, minor axis",
'E_DC_Min': "Error on the deconvolved minor axis size",
'DC_PA': "Deconvolved position angle",
'E_DC_PA': "Error on the deconvolved position angle",
'Isl_rms': "Local rms noise",
'S_Code': "`S' indicates an isolated source which is fit with a single Gaussian; `C' represents sources that are fit by a single Gaussian but are within an island of emission that also contains other sources; and `M' is used for sources which are extended and fitted with multiple Gaussians; `Z' is used for sources that have been (re)-grouped after visual inspection",
'Ngaus': "Number of Gaussian components",
'Size': " Final source size",
'Width': " Final source width",
'Angle': " Final source position angle",
'FLAG_MASKED': "Flag indicating whether source was within the deconvolution mask",
'confused': "Flag indicating whether the source matches multiple HBA sources",
'Deblended_From': "Name of PyBDSF source from which the source was deblended",
'HBA_Deep_Name': "Name of matched HBA Deep source",
'HBA_Deep_Sep': "Distance between LBA and HBA positions"
}


mcatin.keep_columns(keep_cols)
tt = mcatin[keep_cols]

outname = 'bootes_deep_lba.cat.hba_matched.v1.0.fits'
tt.write(outname,overwrite=True)

tt['HBA_Deep_Sep'].unit = 'arcsec'

for tk in tt.colnames:
    tt[tk].description = descriptions[tk]

with open(outname.replace('.fits','.README'),'w') as f:
    f.write(f'# Description of columns in {outname}\n')
    f.write('# \n')
    f.write('# Column Name   | Unit      | Dtype | Description\n')
    f.write('# \n')
    for tk in tt.colnames:
        if tt[tk].unit is not None:
            f.write('{name:15s} | {unit:9s} | {dt:5s} | {desc}'.format(name=tk, dt=str(tt[tk].dtype).replace('|',''), desc= tt[tk].description, unit=tt[tk].unit))
        else:
            f.write('{name:15s} | {unit:9s} | {dt:5s} | {desc}'.format(name=tk, dt=str(tt[tk].dtype).replace('|',''), desc= tt[tk].description, unit=''))
        f.write('\n')
