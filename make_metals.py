import numpy as np
from math import *
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import QTable, Table, Column
from astropy import units as u
from scipy import interpolate
from scipy.interpolate import UnivariateSpline

from constants import *
import peeples2014 as P

def make_metals(galaxy, **kwargs):
    print "about to make metals. whereby metals, I mean oxygen."

    Nspline = kwargs.get("Nspline",10)
    given_age = galaxy['age'].to(u.yr).value
    indices = np.append(np.arange(0,len(galaxy),len(galaxy)/(Nspline-1)), len(galaxy)-1)
    short_age = given_age[indices]

    ## first, how many metals has the galaxy produced?
    dMoxyiidt = Y_O_II * galaxy['sfr']
    Moxymade = np.cumsum(dMoxyiidt * galaxy['dt'])
    ## Need to account for metals made by stars existing in 0th timestep; assume from eq(2) of peeples14
    Moxyinit = np.power(10.0,P.oxyii(np.log10(galaxy['Ms'][0].value))) * u.Msun
    galaxy['Moxymade'] = Moxymade + Moxyinit
    galaxy['dMoxymadedt'] = dMoxyiidt

    ## how much oxygen is in the ISM?
    if "tlogoh" not in galaxy.colnames:
        print "--->you don't have any ISM oxygen! that's OK! I'll make some assuming the FMR!"
        tlogoh = P.fmr(np.log10(np.array(galaxy['Ms'])), np.log10(np.array(galaxy['sfr'])))
        galaxy['tlogoh'] = tlogoh
    Zism = (O_AMU / HELIUM_CORR) * np.power(10.0,galaxy['tlogoh'] - 12)
    Moxyism = Zism * galaxy['Mg']
    galaxy['Zism'] = Zism
    galaxy['Moxyism'] = Moxyism

    ## what is rate at which ISM oxygen mass is changing?
    dMoxyism = UnivariateSpline(short_age,galaxy['Moxyism'][indices].value/1.e9,k=3)
    dMoxyismdt_func = dMoxyism.derivative(n=1)
    dMoxyismdt = dMoxyismdt_func(galaxy['age'].value) * u.Msun / u.yr
    galaxy['dMoxyismdt'] = dMoxyismdt

    ## gotta put some metals into the stars
    print """----make_metals assumes that the galaxy started with stellar metallicity
        equal to whatever the initial gas phase metallicity was;
        at the earliest time you've got a stellar mass of""",galaxy['Ms'][0],"""with Zism=""",galaxy['Zism'][0]
    alpha = np.zeros(len(galaxy))
    Mzstar = np.zeros(len(galaxy)) * u.Msun
    ## currently going to mass-weight the metals
    for i in range(len(galaxy)):
        weight = np.zeros(len(galaxy)) * u.Msun
        for j in range(i+1):
            time_elapsed = galaxy['age'][i] - galaxy['age'][j]
            ## this actually needs to be the Mrecy for consistency right?
            weight[j] = galaxy['sfr'][j]*galaxy['dt'][j] * (1 - fml(time_elapsed))
        alpha[i], sum_weight = np.average(galaxy['Zism'],weights=weight,returned=True)
        Mzstar[i] = np.sum(galaxy['Zism']*weight) ## / np.cumsum(weight)
    galaxy['Zstar'] = alpha
    galaxy['Mzstar'] = Mzstar

    ## CGM?
    if "Mcgm" in galaxy.columns:
        print "assuming that all metals ever made remain in halo; if not in galaxy, then in CGM!"
        print "this breaks if Moxymade isn't tracked by the stellar mass  buildup"
        assert(np.greater_equal(galaxy['Moxymade'],galaxy['Moxyism'] + galaxy['Mzstar']).any()), ">>>>>>> galaxy has made less oxygen than is found in stars and ISM ..."
        galaxy['Mzcgm'] = galaxy['Moxymade'] - galaxy['Moxyism'] - galaxy['Mzstar']
        galaxy['Zcgm'] = galaxy['Mzcgm'] / galaxy['Mcgm']
    else:
        print "no CGM :-( setting Zcgm = 0"
        galaxy['Zcgm'] = 0.0 * galaxy['Zism']

    return galaxy


#-----------------------------------------------------------------------------------------------------


def fml(t):
    return (FML_C0 * np.log(t/FML_LAMBDA + 1))
