import numpy as np
import matplotlib.pyplot as plt
import utils.plotting as pp
import shapely
from shapely.geometry import Polygon, Point
from shapely.ops import cascaded_union

regfilein = '../image_full_phase.tessel.reg'
regfileout = 'facets_fin.reg'

with open(regfilein, 'r') as f:
    flines = f.readlines()
    
plt.figure()
points = []
pimfacets = []
pfacets = []
lines = []
limfacets = []
lfacets = []
cals = []
f = ''
s = ''
fi = 0
centres = []
for l in flines:
    if 'point' in l:
        fi+= 1
        #if fi==15: break
        C = l[l.find('(')+1:l.find(')')]
        C = C.split(',')
        ra = float(C[0])
        dec = float(C[1])
        label = l[l.find('{[')+2:l.find(']}')]
        C = label.split('_')
        f = C[0]
        s = C[1]
        points.append([ra,dec])
        pfacets.append(s)
        pimfacets.append(f)
        plt.text(ra,dec, f+'_'+s)
        plt.plot(ra,dec, 'kx')
        centres.append(Point(ra,dec))
    elif 'line' in l:
        C = l[l.find('(')+1:l.find(')')]
        C = C.split(',')
        ra1 = float(C[0])
        dec1 = float(C[1])
        ra2 = float(C[2])
        dec2 = float(C[3])
        lines.append([ra1, dec1])
        lines.append([ra2, dec2])
        lfacets.append(s)
        lfacets.append(s)
        limfacets.append(f)
        limfacets.append(f)
        print(ra1,dec1,ra2,dec2,s,f)
        plt.plot([ra1,ra2],[dec1,dec2])
        plt.plot([ra1,ra2],[dec1,dec2],'.')
        
lines = np.array(lines)
lfacets = np.array(lfacets)
limfacets = np.array(limfacets)

lofacets = lfacets.copy()
loimfacets = limfacets.copy()


## correct the facet names.... somehow this is wrong in the region file...
facets = np.unique(lofacets)
for fi,f in enumerate(facets):
    thisfacet_ims = loimfacets[(lofacets == f)]
    thisfacet_ims_u = np.unique(thisfacet_ims)
    ps = []
    for si, s in enumerate(thisfacet_ims_u):
        ind = (lofacets == f) & (loimfacets == s)
        linesf = lines[ind]
        p = Polygon(linesf)
        for ci, c in enumerate(centres):
            if p.contains(c):
                print(f,s, pfacets[ci], pimfacets[ci])
                lfacets[ind] = pfacets[ci]
                limfacets[ind] = pimfacets[ci]




facets = np.unique(lfacets)

#sys.exit()
plt.figure()

sout = []
for fi,f in enumerate(facets):
    thisfacet_ims = limfacets[(lfacets == f)]
    thisfacet_ims_u = np.unique(thisfacet_ims)
    ps = []
    for si, s in enumerate(thisfacet_ims_u):
        linesf = lines[(lfacets == f) & (limfacets == s)]
        plt.scatter(linesf.T[0], linesf.T[1],marker='.',c='C'+str(fi))
        plt.plot(linesf.T[0], linesf.T[1],c='C'+str(fi))
        p = Polygon(linesf)
        ps.append(p)
    pu = cascaded_union(ps)
    po = pu.convex_hull
    px,py = po.boundary.coords.xy
    s = 'polygon('
    for pxi,pyi in zip(px,py):
        s+= str(pxi)+','+str(pyi)+','
    s = s[:-1]
    s+= ')'
    print(s)
    sout.append(s)
    
    
    
headout = '''# Region file format: DS9 version 4.1
global color=black width=1
fk5
'''
with open(regfileout,'w') as f:
    f.write(headout)
    for s in sout:
        f.write(s+'\n')
    
    
