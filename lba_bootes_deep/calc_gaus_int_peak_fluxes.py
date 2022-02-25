import numpy as np
import numexpr as ne
import re

import utils.plotting as pp

import collections
import functools

class memoized(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)
      
@memoized
def gaussian(xsize,ysize,x0,y0,sx,sy,pa):
    X,Y = np.meshgrid(np.arange(0,xsize,1.0), np.arange(0,ysize,1.0))
    pa*=np.pi/180.0
    a=0.5*((np.cos(pa)/sx)**2.0+(np.sin(pa)/sy)**2.0)
    b=0.25*((-np.sin(2*pa)/sx**2.0)+(np.sin(2*pa)/sy**2.0))
    c=0.5*((np.sin(pa)/sx)**2.0+(np.cos(pa)/sy)**2.0)
    
    return ne.evaluate('exp(-(a*(X-x0)**2.0+2*b*(X-x0)*(Y-y0)+c*(Y-y0)**2.0))')


def restore_gaussian(image,norm,x,y,bmaj,bmin,bpa,guard,verbose=False):
    # deal with real image co-ords
    yd,xd=np.shape(image)
    if x is int:
        xp=x
        x=x+0.5
        yp=y
        y=y+0.5
    else:
        xp=int(np.trunc(x))
        yp=int(np.trunc(y))
    if verbose:
       print(x,y,xp,yp,norm)
    if xp<0 or yp<0 or xp>xd-1 or yp>yd-1:
       raise Exception('position out of range')
    xmin=xp-guard
    xmax=xp+guard
    ymin=yp-guard
    ymax=yp+guard
    if xmin<0:
        xmin=0
    if ymin<0:
        ymin=0
    if xmax>=xd:
        xmax=xd-1
    if ymax>=yd:
        ymax=yd-1
    x0=x-xmin
    y0=y-ymin
    image[ymin:ymax,xmin:xmax]+=norm*gaussian(xmax-xmin,ymax-ymin,x0,y0,bmaj,bmin,bpa)



ps = 3.  # in arcsec
BMAJ = 15.  # in arcsec
BMIN = BMAJ    #in arcsec

BMAJ = BMAJ/ps
BMIN = BMIN/ps


sizes = np.linspace(15.,50.,10)/ps

peak = 1.
totals = np.zeros_like(sizes)
totals2 = np.zeros_like(sizes)
for si,s in enumerate(sizes):
    image = np.zeros((501,501))
    restore_gaussian(image,peak,250,250,s,s,0,10,verbose=False)
    totals[si] = np.sum(image)
    
    image2 = np.zeros((501,501))
    restore_gaussian(image2,peak,250,250,s,BMAJ,0,10,verbose=False)
    totals2[si] = np.sum(image2)
    #print(s*3600., peak, total)

peaks = 1./totals
peaks2 = 1./totals2

f,ax = pp.paper_single_ax()
ax.plot(sizes*ps, peaks, label='circle')
ax.plot(sizes*ps, peaks2, label='ellipse')
pp.set_attrib(ax, xlabel='size (arcsec)', ylabel='peak intensity')
