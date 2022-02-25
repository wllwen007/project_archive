import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import utils.plotting as pp
import pyrap.quanta as qa
import pyrap.tables as pt
import pyrap.measures as pm
import sys
import os
import numpy as np
import glob


targets = [ {'name' : 'CasA', 'ra' : 6.123487680622104,  'dec' : 1.0265153995604648},
            {'name' : 'CygA', 'ra' : 5.233686575770755,  'dec' : 0.7109409582180791},
            {'name' : 'TauA', 'ra' : 1.4596748493730913, 'dec' : 0.38422502335921294},
            {'name' : 'HerA', 'ra' : 4.4119087330382163, 'dec' : 0.087135562905816893},
            {'name' : 'VirA', 'ra' : 3.276086511413598,  'dec' : 0.21626589533567378},
            {'name' : 'Sun'},
            {'name' : 'Jupiter'},
            {'name' : 'Moon'}]

min_separation = 30
outputimage = None

msname = '/data1/wwilliams/surveys/postcal/mss/L667876_BAND20.ms'

# Create a measures object
me = pm.measures()

# Open the measurement set and the antenna and pointing table
ms = pt.table(msname)  


# Get the position of the first antenna and set it as reference frame
ant_table = pt.table(msname + '::ANTENNA')  
ant_no = 0
pos = ant_table.getcol('POSITION')
x = qa.quantity( pos[ant_no,0], 'm' )
y = qa.quantity( pos[ant_no,1], 'm' )
z = qa.quantity( pos[ant_no,2], 'm' )
position =  me.position( 'wgs84', x, y, z )
me.doframe( position )
ant_table.close()

# Get the first pointing of the first antenna
field_table = pt.table(msname + '::FIELD')
field_no = 0
direction = field_table.getcol('PHASE_DIR')
ra = direction[ ant_no, field_no, 0 ]
dec = direction[ ant_no, field_no, 1 ]
targets.insert(0, {'name' : 'Pointing', 'ra' : ra, 'dec' : dec})
field_table.close()

# Get a ordered list of unique time stamps from the measurement set
time_table = pt.taql('select TIME from $1 orderby distinct TIME', tables = [ms])
time = time_table.getcol('TIME')
time1 = time/3600.0
time1 = time1 - np.floor(time1[0]/24)*24

ra_qa  = qa.quantity( targets[0]['ra'], 'rad' )
dec_qa = qa.quantity( targets[0]['dec'], 'rad' )
pointing =  me.direction('j2000', ra_qa, dec_qa)

separations = []

f,ax = pp.paper_single_ax()

print('SEPARATION from A-Team sources')
print('------------------------------')
print('The minimal accepted distance to an A-Team source is: ' + str(min_separation) + ' deg.')
for target in targets:

    t = qa.quantity(time[0], 's')
    t1 = me.epoch('utc', t)
    me.doframe(t1)

    if 'ra' in target.keys():
        ra_qa  = qa.quantity( target['ra'], 'rad' )
        dec_qa = qa.quantity( target['dec'], 'rad' )
        direction =  me.direction('j2000', ra_qa, dec_qa)
    else :
        direction =  me.direction(target['name'])
    
    separations.append(me.separation(pointing, direction))

    # Loop through all time stamps and calculate the elevation of the pointing
    el = []
    for t in time:
        t_qa = qa.quantity(t, 's')
        t1 = me.epoch('utc', t_qa)
        me.doframe(t1)
        a = me.measure(direction, 'azel')
        elevation = a['m1']
        el.append(elevation['value']/np.pi*180)
    
    el = np.array(el)
    ax.plot(time1, el, label=target['name'])
    
plt.legend()
pp.set_attrib(ax, ylabel='Elevation (deg)', xlabel='Time (UTC)')
pp.fig_save_many(f,'beerze_plots/elevation1')



msnames = glob.glob('/data1/wwilliams/surveys/postcal/mss/L*_BAND15.ms')
msnames.sort()
f,ax = pp.paper_single_ax()
for msname in msnames :
    name = msname.split('/')[-1].split('_')[0]
    

    # Create a measures object
    me = pm.measures()
    # Open the measurement set and the antenna and pointing table
    ms = pt.table(msname)  


    # Get the position of the first antenna and set it as reference frame
    ant_table = pt.table(msname + '::ANTENNA')  
    ant_no = 0
    pos = ant_table.getcol('POSITION')
    x = qa.quantity( pos[ant_no,0], 'm' )
    y = qa.quantity( pos[ant_no,1], 'm' )
    z = qa.quantity( pos[ant_no,2], 'm' )
    position =  me.position( 'wgs84', x, y, z )
    me.doframe( position )
    ant_table.close()

    # Get the first pointing of the first antenna
    field_table = pt.table(msname + '::FIELD')
    field_no = 0
    direction = field_table.getcol('PHASE_DIR')
    ra = direction[ ant_no, field_no, 0 ]
    dec = direction[ ant_no, field_no, 1 ]
    target = {'name' : 'Pointing', 'ra' : ra, 'dec' : dec}
    field_table.close()

    
    # Get a ordered list of unique time stamps from the measurement set
    time_table = pt.taql('select TIME from $1 orderby distinct TIME', tables = [ms])
    time = time_table.getcol('TIME')
    time1 = time/3600.0
    time1 = time1 - np.floor(time1[0]/24)*24
    
    
    ra_qa  = qa.quantity( target['ra'], 'rad' )
    dec_qa = qa.quantity( target['dec'], 'rad' )
    pointing =  me.direction('j2000', ra_qa, dec_qa)

    direction =  me.direction('j2000', ra_qa, dec_qa)

    t = qa.quantity(time[0], 's')
    t1 = me.epoch('utc', t)
    me.doframe(t1)
    

    # Loop through all time stamps and calculate the elevation of the pointing
    el = []
    for t in time:
        t_qa = qa.quantity(t, 's')
        t1 = me.epoch('utc', t_qa)
        me.doframe(t1)
        a = me.measure(direction, 'azel')
        elevation = a['m1']
        el.append(elevation['value']/np.pi*180)
    
    el = np.array(el)
    
    
    if name == 'L667882':
        ls = ':'
    else:
        ls = '-'
    ax.plot(time1, el, label=name, linestyle=ls)

plt.legend()
pp.set_attrib(ax, ylabel='Elevation (deg)', xlabel='Time (UTC)')
pp.fig_save_many(f,'beerze_plots/elevation_all')
