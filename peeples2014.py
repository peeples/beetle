import numpy as np
from constants import *

###### this matches what is in metals.c from the Peeples et al. (2014) metal census paper #####
###### in each case, the input variable "x" is log Mstar in solar masses at redshift = 0
###### the returned value is in units of log Msun
###### assumes Chabrier (2003) IMF everywhere


## kefits() from Kewley & Ellison (2008) on a Chabrier IMF
## average of all of the coefficients except two lowest ones
## as calculated for Peeples et al. (2014)
## log Moxy/Mgas = a + b*x + c*x^2 + d*x^3, with x=log Mstar
KE_a = 27.8612
KE_b = -7.05649
KE_c = 0.818368
KE_d = -0.0302926
DUST_CORR = 0.1 ### assumed assumed dust depletion correction in O/H calibrations

def Fg(x):
    """F_g = Mgas / Mstar"""
    return 0.5*(np.power(10.0,(-0.43*x + 3.75 + np.log10(HE_CORR))) + np.power(10.0,(-0.4814*x + 4.3676)))

def metalslost(x):
    """metals lost = metals made - metals in stars - metals in ISM """
    value = np.power(10.0,metalsmade(x)) - np.power(10.0,starz(x)) - np.power(10.0,ismz(x)) - np.power(10.0,dustz(x))
    return np.log10(value)

def oxylost(x):
    """oxygen lost = oxygen made - oxygen in stars - oxygen in ISM """
    value = np.power(10.0,oxymade(x)) - np.power(10.0,staroxy(x)) - np.power(10.0,ismoxy(x)) - np.power(10.0,dustoxy(x))
    return np.log10(value)

def metalsmade(x):
    return np.log10(np.power(10.0,zii(x)) + np.power(10.0,zia(x)) + np.power(10.0,zagb(x)))

def zii(x):
    return (1.0146*x + np.log10(Y_Z_II) + 0.109109)

def zia(x):
    return (1.043*x - 2.678) ## t^-1 DTD

def zagb(x):
    m = x - 10.5
    return (7.5+ 0.611 + 0.914*m - 0.1797*np.sin(m + np.cos(0.7188*m + 0.336*np.sin(0.611 + 1.797*m) + 0.218*m*np.sin(0.611 + 1.797*m) - 0.218*m*m)))

def starz(x):
    if (x >= 9.05):
        lzstar = (x-10.72)/(4.005+4.275*(x-10.72)+2.051*(x-10.72)*(x-10.72)) + 0.04
    else:
        lzstar = -0.11 + 0.40*(x-6) + np.log10(0.019)
    lmz = 1.08*lzstar - 0.16

    return (x + lmz + np.log10(Z_SUN))

def ismz(x):
    return (ismoxy(x) + np.log10(Z_SUN/Z_SUN_O))

def oxymade(x):
    return np.log10(np.power(10.0,oxyii(x)) + np.power(10.0,oxyia(x)) + np.power(10.0,oxyagb(x)))

def oxyii(x):
    return (1.0146*x + np.log10(Y_O_II) + 0.109109)

def oxyia(x):
    return (1.0216*x - 3.37046)

def oxyagb(x):
    m = x - 10.5
    if(x >= 9.79):
        lagb = 6.52532793646306 + m + np.sin(0.558326101347158*m/np.cos(m*m - m) + np.cos(0.707421425956075*m*m/np.cos(0.101415124817729 + -0.558326101347158*m/np.cos(m*m - m)) - m));
        return (-1.0*np.power(10.0,lagb))  ## OXYGEN DESTRUCTION!
    else:
        return (np.pow(10.0, 7.198 + 1.255*m + 0.01446/np.cos(2.212*m) + (0.9399 + m)/(0.5036 + m)))
    return -10

def staroxy(x):
    alphafe =  -0.495 + 0.28*(x - 0.63)/4.52
    return (alphafe + starz(x) + np.log10(Z_SUN_O/Z_SUN))

def ismoxy(x):
    print "--> not checking lMstar range validity for ismoxy !!! <---"
    lMg = np.log10(Fg(x)) + x
    oh = KE_a + KE_b*x + KE_c*x*x + KE_d*x*x*x
    oh = oh - DUST_CORR
    zg = (O_AMU / HELIUM_CORR) * np.power(10.0, oh-12)
    return (np.log10(zg) + lMg)

def dustz(x):
    return (0.864*x - 1.3065)

def dustoxy(x):
    return (dustz(x) + np.log10(DUST_OXY_FRAC))

def ovioxy(x):
    if ((x >= 9.33) and (x < 10.83)):
        return (np.log10(2e7))
    else:
        return -20

def oviz(x):
    return (ovioxy(x) + np.log10(Z_SUN/Z_SUN_O))

def cgmz(x):
    if((x >= 9.33) and (x <= 10.83)):
        return 7.365
    else:
        return -20

def cgmoxy(x):
    return cgmz(x) + np.log10(Z_SUN_O/Z_SUN)

def hotz(x):
    if((x >= 10.37) and (x <= 11.3)):
        return (0.98*x - 2.89)
    else:
        return -20

def hotoxy(x):
    return(hotz(x) + np.log10(Z_SUN_O/Z_SUN))

def igdustz(x):
    if((x >= 9.929) and (x <= 10.429)):
        return (np.log10(5.e7))
    else:
        return -20
def igdustoxy(x):
    return (igdustz(x) + np.log10(DUST_OXY_FRAC))


#-----------------------------------------------------------------------------------------------------

def fmr(lmstar, lsfr):
    ## copied from code for Peeples et al. (2014), probably from Mannucci et al. (2010)?
    ## takes log stellar mass and log star formation rate and returns 12+log(O/H)
    m = lmstar - 10
    return (8.90 + 0.37*m - 0.14*lsfr - 0.19*m*m + 0.12*m*lsfr - 0.054*lsfr*lsfr)
