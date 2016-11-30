from astropy.table import QTable, Table, Column
from astropy.io import ascii
from astropy import units as u
import numpy as np
from astropy.cosmology import WMAP9 as cosmo

from constants import *
from calculate_rates import fml, dfmldt

def create_galaxies(**kwargs):
    """
    Creates a galaxy

    Parameters
    ----------
    lifetime : `float`, optional
    A positive float for how long the galaxy lives, in years

    Returns
    -------
    galaxy : `~astropy.QTable`
        An astropy QTable with galaxy properties.
    """


    print "creating a galaxy!"

    lifetime = kwargs.get("lifetime",3.0e9 * u.yr)
    dt = kwargs.get("dt",1.e7 *u.yr)
    const_sfr = kwargs.get("sfr",1.0 *u.Msun / u.yr)
    init_Mstar = kwargs.get("init_Mstar", 3.0e10 * u.Msun)
    print "these are hardcoded and should be passed in"
    print "lifetime, dt = ", lifetime, dt
    Nstep = int((lifetime / dt).value) + 1
    const_Mg = init_Mstar
    const_Mh = 1.0e12 * u.Msun
    const_tlogoh = T_LOG_OH_SUN
    print "dt = ",dt," for ",lifetime," with SFR = ",const_sfr
    # initialize
    galaxy = initialize_galaxy(Nstep)

    for i in range(Nstep):
        galaxy['age'][i] = i * dt
        galaxy['lookback'][i] = (Nstep - i - 1)*dt
        galaxy['sfr'][i] = const_sfr
        galaxy['Mg'][i] = const_Mg
        galaxy['Mh'][i] = const_Mh
        galaxy['tlogoh'][i] = const_tlogoh
        galaxy['dt'][i] = dt
        Mrecy = CHI_RECY * const_sfr
        galaxy['dMrdt'][i] = Mrecy
        if i == 0:
            galaxy['Ms'][i] = init_Mstar
        else:
            # figure out how much recylced material
            # assume that previous stars also recyling at rate of CHI_RECY*const_sfr
            # Mrecy = 0 * u.Msun / u.yr
            #for j in range(i):
            #    # recycle rate = cumulative sum of (mass formed in previous timestep) * (dfml(time elapsed)/dt)
            #    Mrecy += (const_sfr * dt) * dfmldt(galaxy['age'][i] - galaxy['age'][j])
            ## Mrecy += const_sfr * fml(galaxy['age'][i])
            galaxy['Ms'][i] = galaxy['Ms'][i-1] + const_sfr*dt - galaxy['dMrdt'][i]*dt

        galaxy['Mcgm'][i] = galaxy['Mh'][i]*cosmo.Ob0/cosmo.Om0 - galaxy['Mg'][i] - galaxy['Ms'][i]
    ## check to make sure sorted so earlier is earlier
    galaxy.sort('age')

    return galaxy


def initialize_galaxy(Nstep):
    galaxy = QTable()
    dummy = np.zeros(Nstep)
    col_dummy = Column(name='age',data=dummy,unit=u.yr)
    galaxy.add_column(col_dummy)
    col_dummy = Column(name='Ms',data=dummy,unit=u.Msun)
    galaxy.add_column(col_dummy)
    col_dummy = Column(name='Mg',data=dummy,unit=u.Msun)
    galaxy.add_column(col_dummy)
    col_dummy = Column(name='Mh',data=dummy,unit=u.Msun)
    galaxy.add_column(col_dummy)
    col_dummy = Column(name='Mcgm',data=dummy,unit=u.Msun)
    galaxy.add_column(col_dummy)
    col_dummy = Column(name='dMrdt',data=dummy,unit=u.Msun/u.yr)
    galaxy.add_column(col_dummy)
    col_dummy = Column(name='sfr',data=dummy,unit=u.Msun/u.yr)
    galaxy.add_column(col_dummy)
    col_dummy = Column(name='tlogoh',data=dummy,unit=u.m/u.m)
    galaxy.add_column(col_dummy)
    col_dummy = Column(name='lookback',data=dummy,unit=u.yr)
    galaxy.add_column(col_dummy)
    col_dummy = Column(name='dt',data=dummy,unit=u.yr)
    galaxy.add_column(col_dummy)

    return galaxy
