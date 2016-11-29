import numpy as np
from astropy import units as u

HELIUM_CORR = 1.366  # factor by which to correct H masses for existence of He, as used in Peeples et al. (2014)

## atomic masses
O_AMU = 15.999
FE_AMU = 55.845

####### Solar metallicities, Caffau et al. (2011), used in Peeples et al. (2014)
Z_SUN = 0.0153
## number abundances, 12+log(X/H):
T_LOG_OH_SUN = 8.76
T_LOG_FEH_SUN = 7.52
## mass abundances
Z_SUN_O = (O_AMU / HELIUM_CORR) * np.power(10.0,T_LOG_OH_SUN - 12)
Z_SUN_FE = (FE_AMU / HELIUM_CORR) * np.power(10.0,T_LOG_FEH_SUN - 12)

## yields, as used in Peeples et al. (2014)
Y_Z_II = 0.03
Y_O_II = 0.015

DUST_OXY_FRAC = 0.27  ## as adopted in Peeples et al. (2014)

Zejmax = Y_O_II / 0.2 ## approx. metallicity of Type II SNe ejecta, assumes 20% of Mstar in high mass stars


## CHABRIER (2003) IMF cumulative fraction of mass lass by stellar particle
## Leitner & Kravtsov (2011) and Jungweirt et al. (2001), using buggy SPS
FML_C0_OLD = 0.046
FML_LAMBDA_OLD = 2.76e5*u.yr
## New values from Behroozi, Wechsler, and Conroy (2013):
FML_C0 = 0.05
FML_LAMBDA = 1.4e6*u.yr
CHI_RECY = 0.40  ## approx assymptotic limit for fml(t~5Gyr)
