import numpy as np
from math import *
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import QTable, Table, Column
from astropy import units as u

from constants import *

def calculate_rates(galaxy, **kwargs):
    print "calculating rates...."

    print """****calculate_rates assumes that galaxy masses are in Msun.
****actual units are: galaxy['Ms'].unit = """ + str(galaxy['Ms'].unit) + """ and galaxy['Mg'].unit = """ + str(galaxy['Mg'].unit)

    Nspline = kwargs.get("Nspline",10)

    given_age = galaxy['age'].to(u.yr).value
    indices = np.append(np.arange(0,len(galaxy),len(galaxy)/(Nspline-1)), len(galaxy)-1)
    short_age = given_age[indices]

    assert ("Mg" in galaxy.colnames), "galaxy column 'Mg' not found, galaxy lacks gas :-("
    print galaxy['Mg'][indices].value
    dMg = UnivariateSpline(short_age,galaxy['Mg'][indices].value,k=3)
    dMgdt_func = dMg.derivative(n=1)
    dMgdt = dMgdt_func(given_age) * u.Msun / u.yr
    col_dMgdt = Column(data=dMgdt,name='dMgdt')
    galaxy.add_column(col_dMgdt)


    if "Mhi" in galaxy.colnames:
        dMhi = UnivariateSpline(short_age,galaxy['Mhi'][indices].value,k=3)
        dMhidt_func = dMhi.derivative(n=1)
        dMhidt = dMhidt_func(given_age) * u.Msun / u.yr
        col_dMhidt = Column(data=dMhidt,name='dMhidt')
        galaxy.add_column(col_dMhidt)

    if "Mh2" in galaxy.colnames:
        dMh2 = UnivariateSpline(short_age,galaxy['Mh2'][indices].value,k=3)
        dMh2dt_func = dMh2.derivative(n=1)
        dMh2dt = dMh2dt_func(given_age) * u.Msun / u.yr
        col_dMh2dt = Column(data=dMh2dt,name='dMh2dt')
        galaxy.add_column(col_dMh2dt)

    assert ("Ms" in galaxy.colnames), "galaxy column 'Ms' not found, galaxy that lacks stars is not a galaxy :-("
    dMs = UnivariateSpline(short_age,galaxy['Ms'][indices].value,k=3)
    dMsdt_func = dMs.derivative(n=1)
    dMsdt = dMsdt_func(given_age) * u.Msun / u.yr
    col_dMsdt = Column(data=dMsdt,name='dMsdt')
    galaxy.add_column(col_dMsdt)

    if "Mh" in galaxy.colnames and "dMhdt" not in galaxy.colnames:
        dMh = UnivariateSpline(short_age,galaxy['Mh'][indices].value,k=3)
        dMhdt_func = dMh.derivative(n=1)
        dMhdt = dMhdt_func(given_age) * u.Msun / u.yr
        col_dMhdt = Column(data=dMhdt,name='dMhdt')
        galaxy.add_column(col_dMhdt)
    if "dMhdt" in galaxy.colnames:
        galaxy['dMcgmdt'] = galaxy['dMhdt'] - galaxy['dMsdt'] - galaxy['dMgdt']

    if "dt" not in galaxy.colnames:
        dt = np.zeros(len(galaxy))*u.yr
        for i in range(len(galaxy)):
            dt[i] = calculate_timestep(galaxy['age'].to(u.yr),i,len(galaxy))
        col_dt = Column(data=dt, name='dt')
        galaxy.add_column(col_dt)

    return galaxy

#-----------------------------------------------------------------------------------------------------

def calculate_return(galaxy):
    if "dMrdt" in galaxy.columns:
        print "dMrdt is already defined; not calculating"
        return galaxy

    ## if dMrdt is not in galaxy.columns then .....:
    print "calculating recycled masses...."
    dMrdt = np.zeros(len(galaxy))*u.Msun / u.yr
    dMrdt_expl = np.zeros(len(galaxy))*u.Msun / u.yr
    Mr = np.zeros(len(galaxy))*u.Msun  ## cumulative mass recyclced
    mass_made = 0
    mass_made_z0 = 0
    recy_prev = 0*u.Msun
    for i in range(len(galaxy)):
        ## how much mass is made in this timestep
        mass_made += galaxy['sfr'][i]*galaxy['dt'][i]
        mass_made_z0 += galaxy['sfr'][i]*galaxy['dt'][i] * (1 - fml(galaxy['lookback'][i]))

        ## time needs to be until this timestep, not z=0 and rates need to be cumulative
        for j in range(i):
            DT = galaxy['dt'][j]
            time_elapsed = galaxy['age'][i] - galaxy['age'][j]
            dMrdt[i] += galaxy['sfr'][j]*DT * dfmldt(time_elapsed)
            Mr[i] += galaxy['sfr'][j]*DT * (1 - fml(time_elapsed))
        dMrdt_expl[i] = (Mr[i] - recy_prev) / galaxy['dt'][i]
        recy_prev = Mr[i]

    col_dMrdt = Column(data=dMrdt, name='dMrdt')
    galaxy.add_column(col_dMrdt)
    col_dMrdt_expl = Column(data=dMrdt_expl, name='dMrdt_expl')
    galaxy.add_column(col_dMrdt_expl)
    col_Mr = Column(data=Mr, name='Mr')
    galaxy.add_column(col_Mr)
    return galaxy


#-----------------------------------------------------------------------------------------------------

def calculate_etas(galaxy):
    print "calculating etas...."

    # eta_aw = eta_acc - eta_w
    if "dMrdt" in galaxy.columns:
        print galaxy
        eta_aw = (galaxy['dMgdt'] - galaxy['sfr'] + galaxy['dMrdt']) / galaxy['sfr']
    else:
        print "assuming dMs/dt = SFR - dMrecy/dt, which if mergers, it doesn't"
        eta_aw = (galaxy['dMgdt'] - galaxy['dMsdt']) / galaxy['sfr']

    # figure out metal rates

    ## first figure out metallicities
    print "assuming that ISM is well mixed before outlow; setting fe = 1"
    fe = 1
    Zw = (1-fe)*Zejmax + fe*galaxy['Zism']
    galaxy['Zw'] = Zw
    if "Zcgm" in galaxy.columns:
        Za = galaxy['Zcgm']
    else:
        print "No CGM metallicity? Assuming no metal accretion."
        Za = 0

    # zeta_aw = zeta_a - zeta_w
    print "assuming metals from stellar mass loss is zero"
    dMroxydt = galaxy['dMoxymadedt']
    zeta_aw = (galaxy['dMoxyismdt'] + dMroxydt)/galaxy['sfr'] - galaxy['Zism']
    etaa = (zeta_aw + Zw*eta_aw) / (galaxy['Zcgm'] + Zw)
    etaw = etaa - eta_aw
    zetaa = (galaxy['Zcgm']/galaxy['Zism']) * etaa
    zetaw = (galaxy['Zw']/galaxy['Zism']) * etaw

    # dMoutoxydt = (Y_O_II - galaxy['Zoxy'])*galaxy['sfr']*u.Msun/u.yr - galaxy['dMoxyismdt'] - dMoxyaccdt
    # dMwdt = dMoutoxydt / Zw
    # etaw = dMwdt / galaxy['sfr']
    # etaa = eta_est + etaw

    galaxy['etaw'] = etaw
    galaxy['etaa'] = etaa

    return galaxy


#-----------------------------------------------------------------------------------------------------

def fml(t):
    "returns cumulative fraction of mass loss by population of stars formed time t ago"
    return (FML_C0 * np.log(t/FML_LAMBDA + 1))



#-----------------------------------------------------------------------------------------------------

def dfmldt(t):
    return (FML_C0 / (t + FML_LAMBDA))



#-----------------------------------------------------------------------------------------------------

def calculate_timestep(age, i, length):
    if i == 0:
        DT = 0.5*age[i+1]
    elif i == length-1:
        DT = 0.5*(age[i] - age[i-1])
    else:
        DT = 0.5*(age[i+1] - age[i-1])
    return DT
