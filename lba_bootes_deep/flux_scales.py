import numpy as np
import matplotlib.pyplot as plt

'''# S&H: scale to RCB below 325MHz (avoids issues with CasA dec at low flux densities in usual Baars77)
# ratio of Baars 1977 to KPW @ 1400 MHz is 1.029

nu > 325: RCB and KPW consistent
 Roger, Costain & Bridle (1973)  - RCB
 Kellerman, Pauliny-Toth & Williams (1969) - KPW 
nu < 325: need to scale from KPW to RCB 
'''
 #B77 data needs to be scaled to RCB: using polynomial fit to ratios of Baars to KPW
ff = [38,  
81.5,
178,
400,
750,
1400,
2700,
5000,
10600,
15000]
frat = [0.981, 1.020, 1.051, 1.065, 1.046, 1.029, 1.011, 0.993, 0.974, 0.966]
# wenss : average scaling to B77, then to KPW
# 6C already on RCB
# VLA-P on SH


#RACS Perley & Butler (2017) and Reynolds (1994).

plt.figure()
plt.scatter(ff, frat)


p = np.polyfit(np.log10(ff),frat,5) 
ffrange = np.geomspace(ff[0],ff[-1],100)
fratrange = np.polyval(p, np.log10(ffrange)) 
plt.plot(ffrange, fratrange)
np.polyval(p, np.log10(887)) 
np.interp(np.log10(887),np.log10(ff),frat )

# 1.041

# S&H 3C295 best model @ 325MHz
A = [97.763, -0.582, -0.298, 0.583, -0.363]
nu = 325.  
S = A[0] * np.prod(np.array([ 10**(A[i]*(np.log10(nu/150.)**i)) for i in range(1,4) ]))

# 60.7038275160885

# in WENSS is 61.840   -> WENSS flux should be x 0.9816272237401116


'''
From Risely+ 2016:
It is important to note that in all our comparisons, the flux density measurements from the literature have been adjusted to bring
them in line with the Scaife & Heald (2012) flux scale adopted in
this work (hereafter SH12). The flux density scale of WENSS is
complex, having been set using observations of the sources 3C 48,
3C 147, 3C 286 and 3C 295. The correction factor to scale the
WENSS survey to the SH12 flux scale is an average correction
across the discrete set of WENSS calibrators. It would be incorrect to apply this factor to any localised area of the sky as the difference between the flux densities of the WENSS calibrators and
the SH12 source models range from ∼ 1 per cent to ∼ 18 per cent
(see Rengelink et al. 1997 and SH12). We note that, given the small
difference between the native WENSS flux densities in this region
and the GMRT values measured here (see §4.1.3) as well as the RA
range, this field was likely calibrated using either 3C 147 or 3C 48.
These calibrators both have ratios to the SH12 scale within the measurement uncertainty. Hence, we do not apply a correction factor to
the flux densities from WENSS during catalogue verification
'''
