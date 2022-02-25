import os
import glob
import numpy as np
from astropy.table import Table, Column
from utils.catalogues import Make_Shape, name_from_coords

clobber = True



dpath = '/data1/wwilliams/surveys/postcal/ddf/outless5C/DEEP-obs7/fin_im/'

cat = 'bootes_deep_lba.cat.fits'
gcat = 'bootes_deep_lba.gcat.fits'
catout = 'bootes_deep_lba.cat.mergedeblend.fits'

tcat = Table.read(dpath+cat)
tgcat = Table.read(dpath+gcat)

#deblend_files = glob.glob('deblend/*.txt')

deblendlist = 'to_deblend_list.txt'
with open(dpath+deblendlist, 'r') as f:
    l = f.readlines()
to_deblend = np.array([ll.strip() for ll in l])

mergefiles = glob.glob('LBABO*merge.txt')

tcat.add_column(Column(data=np.zeros(len(tcat),dtype=bool),name='MergeDup'))
tcat.add_column(Column(data=np.zeros(len(tcat),dtype=bool),name='Merge_Flag'))
tcat.add_column(Column(data=np.nan*np.zeros(len(tcat),dtype=bool),name='Merge_Size'))
tcat.add_column(Column(data=np.nan*np.zeros(len(tcat),dtype=bool),name='Merge_Width'))
tcat.add_column(Column(data=np.nan*np.zeros(len(tcat),dtype=bool),name='Merge_PA'))

if clobber:
    
    mergelist = []
    for mergei, mergef in enumerate(mergefiles):
        with open(mergef, 'r') as f:
            l = []
            for line in f:
                li = line.strip()
                l.append(li)
            l.sort()   # sort so we can check for duplicate lists - still need to check for duplicates in lists
            mergelist.append(l)

    mergelist = np.array(mergelist)
    umergelist = np.unique(mergelist)
    print(len(mergelist),'groups for merging')
    print(len(umergelist),'unique groups for merging')

    allmergelist = []
    for mergegroup in umergelist:
        for mergei in mergegroup:
            allmergelist.append(mergei)
    allmergelist = np.array(allmergelist)
    uallmergelist, ucnt = np.unique(allmergelist, return_counts=True)
    print(len(allmergelist),'sources in unique merge groups')
    print(len(uallmergelist),'unique sources in unique merge groups')

    if len(allmergelist) != len(uallmergelist):
        print('!!!!! there is a problem - at least one source appears in multiple mergelists, please resolve before rerunning !!!!!')
        for probsrc in uallmergelist[ucnt>1]:
            print(probsrc)
            for mergegroup in mergelist:
                if probsrc in mergegroup:
                    print(mergegroup)
        sys.exit()
    else:
        print('all sources in mergelist are unique')
        #if len(l) == 1:
        #    print(l,'nothing to merge, continuing')
        #    continue
    #allmergelist = np.array(allmergelist)
    #allmergelisti = np.array(allmergelisti)
    #uallmergelist, ucount = np.unique(allmergelist, return_counts=True)
    #print(len(allmergelist),'sources in full merge list')
    #print(len(uallmergelist),'unique sources in full merge list')

    #dupmergelist = uallmergelist[ucount>1]
    #for di in dupmergelist:
    #    print(allmergelist[allmergelist==di])
    #    print(allmergelisti[allmergelist==di])

    mergelist = []
    for mergegroup in umergelist:
        l = mergegroup
        if len(l) == 1:
            print(l[0], 'size of 1, dont merge')
            continue
        okmerge = True
        for ll in mergegroup:
            if ll not in tcat['Source_Name']:
                print(ll,'not in catalogue')
                print(ll,ll in mergelist)
                if ll not in mergelist:
                    #its not in the catalogue and its not in the already mearged list
                    pass
                okmerge = False

        if not okmerge:
            print('problem with',mergef)
            print(l)
            continue
        inds = np.array([np.where(tcat['Source_Name']==ll)[0][0] for ll in l ])
        print(l, inds)

        for ll in l:
            mergelist.append(ll)

        sh = Make_Shape(tcat[inds])

        tcat['MergeDup'][inds[1:]] = 1 ## set duplicate flag to 1 for all except first one
        tcat['Merge_Flag'][inds] = 1 ## source has been merged

        tcat['Source_Name'][inds] = name_from_coords(sh.ra, sh.dec, prefix='LBABOOJ')
        tcat['RA'][inds] = sh.ra
        tcat['DEC'][inds] = sh.dec
        tcat['E_RA'][inds] = np.nan
        tcat['E_DEC'][inds] = np.nan
        tcat['Isl_rms'][inds] = np.mean(tcat['Isl_rms'][inds])
        tcat['Total_flux'][inds] = np.sum(tcat['Total_flux'][inds])
        # total flux error is error on the sum
        tcat['E_Total_flux'][inds] = np.sqrt(np.sum(tcat['Total_flux'][inds]**2.0))
        # peak flux and error from brightest component
        maxpk=np.argmax(tcat['Peak_flux'][inds])
        tcat['Peak_flux'][inds] = tcat['Peak_flux'][inds[maxpk]]
        tcat['E_Peak_flux'][inds] = tcat['E_Peak_flux'][inds[maxpk]]
            
        tcat['S_Code'][inds] = 'M'

        tcat['Isl_id'][inds] = -99
            
        for t in ['Maj', 'Min', 'PA']:
            tcat[t][inds] = np.nan
            tcat['E_'+t][inds] = np.nan
            tcat['DC_'+t][inds] = np.nan
            tcat['E_DC_'+t][inds] = np.nan
        for t in tcat.colnames:
            if ('_img_plane' in t) or ('_max' in t):
                tcat[t][inds] = np.nan
                        
        tcat['Merge_Size'][inds] = sh.length()
        tcat['Merge_Width'][inds] = sh.width()
        tcat['Merge_PA'][inds] = sh.pa()
        

        'Source_id,Isl_id,RA,E_RA,DEC,E_DEC,Total_flux,E_Total_flux,Peak_flux,E_Peak_flux,RA_max,E_RA_max,DEC_max,E_DEC_max,Maj,E_Maj,Min,E_Min,PA,E_PA,Maj_img_plane,E_Maj_img_plane,Min_img_plane,E_Min_img_plane,PA_img_plane,E_PA_img_plane,DC_Maj,E_DC_Maj,DC_Min,E_DC_Min,DC_PA,E_DC_PA,DC_Maj_img_plane,E_DC_Maj_img_plane,DC_Min_img_plane,E_DC_Min_img_plane,DC_PA_img_plane,E_DC_PA_img_plane,Isl_Total_flux,E_Isl_Total_flux,Isl_rms,Isl_mean,Resid_Isl_rms,Resid_Isl_mean,S_Code,Artefact,Source_Name'

print(len(tcat))
tcat = tcat[tcat['MergeDup']==0]
tcat.remove_column('MergeDup')

print(len(tcat))

tcat.write(catout, overwrite=True)
    
sys.exit()
   
if clobber: 
    for deblend in to_deblend:
        
        gname = deblend
            
        # for a given source to be made from a list of Gaussian components
        compi = np.where(tgcat['Source_Name']==gname)[0]
        
        print(deblend, len(compi))
        #sh = Make_Shape(gcat[compi])
    
    
    #tcat.write(dpath+catout, overwrite=True)

'''
## to merge
LBABOOJ144301.14+350842.2 ... with??
LBABOOJ143426.30+342812.4
LBABOOJ143436.66+342814.6
LBABOOJ143433.14+352135.4
LBABOOJ143652.29+341637.3 - 4 sources
LBABOOJ143654.38+341728.5
LBABOOJ143525.17+325127.3
LBABOOJ143433.14+352135.4
LBABOOJ143229.82+354559.8
LBABOOJ143124.49+342554.9
LBABOOJ143103.75+334548.4
LBABOOJ143013.52+351858.6
LBABOOJ142946.69+335624.4
LBABOOJ142939.80+335718.1
LBABOOJ142921.22+335425.6
LBABOOJ142913.84+335619.6
LBABOOJ142848.74+323824.7
LBABOOJ142651.54+322227.2
LBABOOJ142636.88+362507.4
LBABOOJ142628.21+332012.9
LBABOOJ143322.24+344956.5
LBABOOJ144301.14+350842.2
LBABOOJ144159.84+340241.1
LBABOOJ143858.70+324954.1
LBABOOJ143756.28+351942.6 - merge artefact?
LBABOOJ143652.29+341637.3 - DDRG
LBABOOJ143656.04+341827.5
LBABOOJ143648.94+341531.5
LBABOOJ144439.79+344438.0 ?
LBABOOJ144337.87+343035.2
LBABOOJ144307.52+334218.6
LBABOOJ144305.49+334146.0
LBABOOJ144249.41+330533.2
LBABOOJ144155.65+340047.0?
LBABOOJ144151.36+335934.4
LBABOOJ144147.82+340035.8
LBABOOJ143911.55+315306.2
LBABOOJ143815.09+351000.8 - no second source to merge
LBABOOJ143751.33+331157.4?
LBABOOJ143747.02+331221.3
LBABOOJ143735.23+343154.4
LBABOOJ143705.27+362615.9
LBABOOJ143704.99+362723.5
LBABOOJ143628.26+343022.6??
LBABOOJ143625.44+343041.8??
LBABOOJ143606.30+314134.5??
LBABOOJ143549.62+331421.4
LBABOOJ143546.93+350520.1
LBABOOJ143543.82+350435.2
LBABOOJ143537.47+350353.1
LBABOOJ143546.91+354818.8
LBABOOJ143542.90+354720.9
LBABOOJ143524.16+325013.7 - art?
LBABOOJ143455.60+315902.8
LBABOOJ143451.68+315944.1
LBABOOJ143451.63+345214.3 ?
LBABOOJ143445.96+332819.5
LBABOOJ143319.02+313046.3
LBABOOJ143313.56+313218.5
LBABOOJ143237.63+354322.6
LBABOOJ143200.94+351751.8
LBABOOJ143158.89+351648.7
LBABOOJ143159.68+321051.8
LBABOOJ143158.76+321203.8
LBABOOJ143138.73+321704.2 ?
LBABOOJ143006.21+352109.6
LBABOOJ143008.00+353601.4
LBABOOJ142936.99+321751.2
LBABOOJ142934.42+321745.1
LBABOOJ142923.50+352852.7
LBABOOJ142735.83+341933.8 ?
LBABOOJ142712.95+362402.4 ?
LBABOOJ142717.59+331226.2 ?
LBABOOJ142629.59+323534.3
LBABOOJ142625.10+333149.9 ?
LBABOOJ142619.84+332007.9
LBABOOJ142558.07+351257.9
LBABOOJ142512.58+350542.9
LBABOOJ142458.43+322320.1?
LBABOOJ142456.13+322231.4?
LBABOOJ142408.62+360859.3
LBABOOJ142402.55+360618.4
LBABOOJ142407.63+345131.9
LBABOOJ142410.28+341706.8
LBABOOJ142410.09+331302.3
LBABOOJ142406.95+331255.8
LBABOOJ142350.36+360510.5
LBABOOJ142307.81+360036.0?
LBABOOJ142248.01+325549.3
LBABOOJ142242.39+325633.8
LBABOOJ142202.02+352111.5
LBABOOJ142159.60+352323.3
LBABOOJ142151.87+352133.6?
LBABOOJ142135.43+351223.1
LBABOOJ142133.08+351544.9
LBABOOJ142131.18+342443.6
LBABOOJ142129.12+342511.3
LBABOOJ142123.48+331625.1
LBABOOJ142002.70+350832.4
LBABOOJ144252.17+330625.7
LBABOOJ144251.74+330607.0
LBABOOJ144259.08+350923.8
LBABOOJ143738.05+343249.7
LBABOOJ143556.95+364732.6
LBABOOJ143433.07+342750.8
LBABOOJ143429.13+342726.5
LBABOOJ143108.41+352531.8
LBABOOJ143108.12+352555.4
LBABOOJ143010.89+353523.3?
LBABOOJ142919.10+352815.4
LBABOOJ142510.90+350733.6
LBABOOJ142127.91+351516.4 - big diffuse
LBABOOJ144441.82+344409.1
LBABOOJ144340.02+342957.2
LBABOOJ144157.11+340200.4
LBABOOJ144002.38+344257.4
LBABOOJ143828.01+354122.1
LBABOOJ143827.67+354048.6
LBABOOJ143553.62+364719.4
LBABOOJ143526.31+325121.2
LBABOOJ143500.57+334419.6
LBABOOJ143453.37+343625.1
LBABOOJ143447.37+332826.3
LBABOOJ143233.55+354416.2
LBABOOJ143227.61+332756.6
LBABOOJ143146.60+323736.0
LBABOOJ143055.15+334448.3
LBABOOJ142906.95+354238.4
LBABOOJ142845.46+323844.0
LBABOOJ142734.86+341857.6
LBABOOJ142737.64+330817.8
LBABOOJ142624.53+323650.5
LBABOOJ142411.00+345150.9
LBABOOJ142352.17+353307.0 - this is a GRG?
LBABOOJ142350.05+353224.8
LBABOOJ142334.10+352736.5
LBABOOJ142135.40+341627.6
LBABOOJ142045.40+353912.3
LBABOOJ142001.20+350844.7
LBABOOJ143748.39+345438.6
LBABOOJ143553.14+331504.3
LBABOOJ143443.19+341220.9
LBABOOJ143359.75+360104.8
LBABOOJ142927.89+352958.5
LBABOOJ142910.65+362000.4
LBABOOJ142627.20+323611.4
LBABOOJ142556.53+351312.5
LBABOOJ142408.56+341713.8
LBABOOJ142310.43+355919.6
LBABOOJ143757.10+352007.6
LBABOOJ143443.02+332924.7
LBABOOJ143434.67+352114.3
LBABOOJ142650.67+322245.5
## to deblend
LBABOOJ142842.74+342535.0 - merge & deblend
LBABOOJ143959.23+325033.6 - m&d - art?
LBABOOJ143449.25+354214.4 m&d
LBABOOJ142838.65+342340.5 m&d
LBABOOJ143426.90+353206.0 - into 2 sources
LBABOOJ143450.47+354252.6
LBABOOJ143641.01+315900.3 - into 2
LBABOOJ143603.23+334352.9
LBABOOJ143528.35+331144.3
LBABOOJ143450.47+354252.6
LBABOOJ143426.90+353206.0
LBABOOJ143225.74+332839.4
LBABOOJ143208.92+335716.1
LBABOOJ143134.76+351506.4
LBABOOJ143130.13+342819.8
LBABOOJ143039.96+314947.7
LBABOOJ142949.25+353438.0
LBABOOJ142948.47+351004.3
LBABOOJ142918.61+343725.4
LBABOOJ142916.69+321726.8
LBABOOJ142842.73+342429.5
LBABOOJ142806.09+363119.6
LBABOOJ142759.86+345458.8
LBABOOJ142759.41+324734.6
LBABOOJ142713.15+350024.4
LBABOOJ142538.61+344858.0
LBABOOJ141942.84+342046.0
LBABOOJ143114.97+365917.2 - one artefact
LBABOOJ144242.78+324253.7 ??
LBABOOJ144229.08+323449.4 - can't actually deblend
LBABOOJ144142.51+334719.7
LBABOOJ144010.05+341332.5
LBABOOJ143958.84+325009.2
LBABOOJ143912.06+321622.6
LBABOOJ143848.76+323529.5 - can't
LBABOOJ143749.76+345458.6
LBABOOJ143641.01+315900.3 - 2 doubles
LBABOOJ143126.45+344111.9?
LBABOOJ143122.90+344024.4?
LBABOOJ143121.94+343931.2?
LBABOOJ142903.06+353823.7 - can't
LBABOOJ142738.33+355058.7 - can't
## interesting
LBABOOJ143445.96+332819.5
LBABOOJ143433.14+352135.4
LBABOOJ143652.29+341637.3
LBABOOJ143639.26+361331.1 - remnant?
LBABOOJ143623.23+352714.8 - possible halo around source?
LBABOOJ143352.27+355957.0
LBABOOJ143239.69+361808.5
LBABOOJ143229.82+354559.8
LBABOOJ143244.73+340737.6 - no hba??
LBABOOJ144159.84+340241.1 - diffuse
LBABOOJ144155.65+340047.0
LBABOOJ143516.05+332706.1 - no hba?
LBABOOJ143229.15+320450.9 - no hba? but star/qso?
LBABOOJ142827.60+363554.1 - artefact?
LBABOOJ142702.74+362143.7 - artefact?
LBABOOJ142658.13+351850.2 - diffuse in hba too
LBABOOJ142522.05+345242.3 - no hba ... but near gal...
LBABOOJ143659.00+333924.9 no hba?
LBABOOJ143455.10+364015.9 - no hba
LBABOOJ143252.34+322107.9 - diffuse blob
LBABOOJ143038.06+320049.9 - no hba
LBABOOJ142935.18+332738.7 - no hba
LBABOOJ143447.37+332826.3 - diffuse halo around source...
LBABOOJ142518.86+355725.1 - int the lobes aren't detected at lba
LBABOOJ142352.17+353307.0 - large source! end of lobe!
LBABOOJ142346.78+360227.0 - int lobes aren't seen at lba
LBABOOJ142625.49+350328.3 diffuse
LBABOOJ142333.11+340147.3 - only one lobe detected...
'''
