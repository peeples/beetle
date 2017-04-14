import numpy as np
from math import pi
from constants import *
import peeples2014 as P
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import QTable, Table, Column
from astropy import units as u
import sys

from calculate_rates import fml, dfmldt
from create_galaxies import initialize_galaxy

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['font.family'] = 'stixgeneral'
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors


def get_alpha(lMstar, mzr, Fgyear):
    '''
    alpha = (1 - frecy) * (dlogMg/dlogMstar + dlogZg/dlogMstar)
    '''
    dlogZg = get_dlogZg(lMstar, mzr)
    dlogMg = get_dlogMg(lMstar, Fgyear)
    frecy = 0.2
    return (1-frecy)*(dlogMg - dlogZg)

def get_dlogMg(lMstar, Fgyear):
    '''
    Fg = Mg/Mstar = Kf*(Mstar)^-gamma
    log Mg = log Kf + (1-gamma)*log Mstar
    dlogMg/dlogMstar = 1-gamma
    '''
    Kf, gamma = get_Fg_params(Fgyear)
    return (1-gamma)

def get_Fg(lMstar, Fgyear):
    '''
    Fg = Mg/Mstar = Kf*(Mstar)^-gamma
    '''
    Kf, gamma = get_Fg_params(Fgyear)
    Mstar = np.power(10.0,lMstar)
    return Kf*np.power(Mstar,-1.0*gamma)

def get_Fg_params(Fgyear):
    '''
    Fg = Mg/Mstar = Kf*(Mstar)^-gamma
    '''
    if Fgyear == "P14":
        ## P14
        return np.power(10,4.39), 0.48
    elif Fgyear == "PS11":
        ## ps11
        return 316228, 0.57

def get_dlogZg(lMstar, mzr):
    '''
    log Moxy/Mgas = a + b*x + c*x^2 + d*x^3, with x=log Mstar
    dlogZg/dlogMstar = b + 2*c*x + 3*d^2
    '''
    if mzr == "t04":
        a, b, c, d = get_Zg_coeffs(mzr)
        return (b + 2*c*lMstar + 3*d*lMstar*lMstar)
    else:
        tlogoh_asm = 8.798
        mto = np.power(10.0, 8.901)
        gamma = 0.640
        mstar = np.power(10.0, lMstar)
        return -1*np.power(mstar/mto, gamma-1) * \
                (gamma/mto) / (1+np.power(mstar/mto,gamma))

def get_tlogoh(lMstar, mzr):
    '''
    12 + log(O/H) = a + b*x + c*x^2 + d*x^3, with x=log Mstar
    '''
    ## first let's use the Tremonti values
    if mzr == "t04" or mzr == "p14":
        a, b, c, d = get_Zg_coeffs(mzr)
        return (a + b*lMstar + c*lMstar*lMstar + d*lMstar*lMstar*lMstar)
    else:
        # Andrews & Martini 2013 table 4
        tlogoh_asm = 8.798
        mto = np.power(10.0, 8.901)
        gamma = 0.640
        mstar = np.power(10.0, lMstar)
        tlogoh = tlogoh_asm - np.log10(1+np.power(mto / mstar, gamma))
        return tlogoh

def get_logZg(lMstar, mzr):
    tlogoh = get_tlogoh(lMstar, mzr)
    # print tlogoh
    logZg = tlogoh - 12 + np.log10((15.999/1.0079)/(1.0079*0.75+4.0026*0.25))
    return logZg

def get_Zg_coeffs(mzr):
    # return KE_a, KE_b, KE_c, KE_d
    ### Tremonti+2004 from KE08,PS11
    if mzr == "t04":
        return -0.759210, 1.30177, 0.003261, -0.00364112
    if mzr == "p14":
        return P.KE_a, P.KE_b, P.KE_c, P.KE_d

def make_galaxy():
    ## integrate Zg(mstar)
    #### take mstar as an array
    ## Zstar = int(Zg, mstar) ?

    ## can be done numerically? or analytically?

    ## let's start numerically because I'm lazy
    sfr = 15.*u.Msun/u.yr # msun per year
    nsteps = 2000
    time_start = 2*u.Gyr
    init_Mstar = 1.e6 * u.Msun
    dt = (cosmo.age(z=0) - time_start) / nsteps
    galaxy = initialize_galaxy(nsteps)
    weight = np.zeros(nsteps) * u.Msun
    zstar = np.zeros(len(galaxy))
    Mzstar = np.zeros(len(galaxy)) * u.Msun
    for i in range(nsteps):
        galaxy['age'][i] = i * dt
        galaxy['dt'][i] = dt
        galaxy['lookback'][i] = (nsteps - i - 1)*dt
        galaxy['sfr'][i] = sfr
        if i == 0:
            galaxy['Ms'][i] = init_Mstar
        else:
            galaxy['Ms'][i] = galaxy['Ms'][i-1] + sfr*dt - sfr*(1-fml(galaxy['age'][i]))*dt

        weight[i] = sfr*(1-fml(galaxy['age'][i]))*dt
        galaxy['tlogoh'][i] = get_tlogoh(np.log10(galaxy['Ms'][i].value), "am13")
        Zism = (O_AMU / HELIUM_CORR) * np.power(10.0,galaxy['tlogoh'] - 12)
        galaxy['Zism'] = Zism
        zstar[i], sum_weight = np.average(galaxy['Zism'],weights=weight,returned=True)
        Mzstar[i] = np.sum(galaxy['Zism']*weight) ## / np.cumsum(weight)

    ## how many metals made?
    dMoxyiidt = Y_O_II * galaxy['sfr']
    Moxymade = np.cumsum(dMoxyiidt * galaxy['dt'])
    ## Need to account for metals made by stars existing in 0th timestep; assume from eq(2) of peeples14
    Moxyinit = np.power(10.0,P.oxyii(np.log10(galaxy['Ms'][0].value))) * u.Msun
    galaxy['Moxymade'] = Moxymade + Moxyinit
    galaxy['dMoxymadedt'] = dMoxyiidt

    galaxy['Zstar'] = zstar
    galaxy['Mzstar'] = Mzstar

    return galaxy

def mstar_mhalo(z, mhalo):
    ## Moster+09, eqns 16,17 in PS11
    mhalo0 = np.power(10, 11.88 + 0.018*(1+z))
    mMz = 0.0282*np.power(1+z, -0.72)
    gammaz = 0.556*np.power(1+z, -0.26)
    betaz = 0.17*z + 1.06
    denom = np.power(mhalo/mhalo0, -1.0*betaz) + np.power(mhalo/mhalo0,gammaz)
    mstar = mhalo * 2*mMz/denom

    ## vvir, eqn 19 PS11
    d = cosmo.Om(z) - 1
    Delta = 18*pi*pi + 82*d - 39*d*d
    factor = np.power((cosmo.Om0 / 0.25) * (1/cosmo.Om(z)) * (Delta / (18*pi*pi)), 1/6.)
    vvir = 112.6 * np.power(mhalo / 1e12, 1/3.) * factor * np.sqrt(1+z)
    # +0.05 for Chabrier
    return (np.log10(mstar) + 0.05), vvir

def ps11():
    # OK, first we want some galaxies
    lMstar = np.arange(8,11.5,0.05)
    ## start with lmhalo
    lMhalo = np.arange(10.5, 13.5, 0.025)
    lMstar, vvir = mstar_mhalo(0.0, np.power(10.0,lMhalo))
    mzr = "t04"
    logZg = get_logZg(lMstar, mzr)
    Zg = np.power(10.0, logZg)
    Fgyear = "PS11"
    Fg = get_Fg(lMstar, Fgyear)
    alpha = get_alpha(lMstar, mzr, Fgyear)
    ## (y/Zg - 1) = zetaw - zetaa + alpha*Fg
    ## zetaw - zetaa = y/Zg - 1 - alpha*Fg
    zetaw_old = Y_O_II/Zg - 1 - alpha*Fg

    Fgyear = "P14"
    Fg = get_Fg(lMstar, Fgyear)
    alpha = get_alpha(lMstar, mzr, Fgyear)
    ## (y/Zg - 1) = zetaw - zetaa + alpha*Fg
    ## zetaw - zetaa = y/Zg - 1 - alpha*Fg
    zetaw_tp = Y_O_II/Zg - 1 - alpha*Fg

    mzr = "p14"
    logZg = get_logZg(lMstar, mzr)
    Zg = np.power(10.0, logZg)
    Fg = get_Fg(lMstar, Fgyear)
    alpha = get_alpha(lMstar, mzr, Fgyear)
    zetaw_p = Y_O_II/Zg - 1 - alpha*Fg

    mzr = "am13"
    logZg = get_logZg(lMstar, mzr)
    Zg = np.power(10.0, logZg)
    Fg = get_Fg(lMstar, Fgyear)
    alpha = get_alpha(lMstar, mzr, Fgyear)
    zetaw_new = Y_O_II/Zg - 1 - alpha*Fg

    Fgyear = "PS11"
    Fg = get_Fg(lMstar, Fgyear)
    alpha = get_alpha(lMstar, mzr, Fgyear)
    ## (y/Zg - 1) = zetaw - zetaa + alpha*Fg
    ## zetaw - zetaa = y/Zg - 1 - alpha*Fg
    zetaw_amps = Y_O_II/Zg - 1 - alpha*Fg

    ## let's plot stuff!
    fig = plt.figure(figsize=(11,8),dpi=300)
    ax = fig.add_subplot(111)
    ax.plot(lMstar,np.log10(zetaw_old),color='blue',lw=3, label="T04, PS11")
    ax.plot(lMstar,np.log10(zetaw_amps),color='red',lw=3, label="AM13, PS11")
    ax.plot(lMstar,np.log10(zetaw_new),color='orange',lw=3, label="AM13, P14")
    ax.plot(lMstar,np.log10(zetaw_tp),color='purple',lw=3, label="T04, P14")
    ax.plot(lMstar,np.log10(zetaw_p),color='green',lw=3, label="P14, P14")
    # ax.plot(vvir,np.log10(zetaw),color='blue',lw=3)
    # ax.plot(lMstar,np.log10(Zg/Y_O_II) ,color='blue',lw=3)

    # ax.set_xscale('log')
    # plt.xlim(40,300)

    plt.xlabel(r"log $M_{\star}$",fontsize=20)
    # plt.xlabel(r"v$_{\rm vir}$",fontsize=20)
    plt.ylabel(r"log $\zeta_{\rm w}$",fontsize=20)
    # plt.ylabel(r"log $Z_{\rm g}/y$",fontsize=20)
    lg = plt.legend(loc='best')
    lg.draw_frame(False)
    plt.tight_layout()

    plt.savefig("zetaw_lMstar.png")
    # plt.savefig("zetaw_vvir.png")
    plt.close(fig)

def plot_metals_retained(galaxy):
    fig = plt.figure(figsize=(8,8),dpi=400)
    # plt.style.use('seaborn-white')
    ax = fig.add_subplot(111)

    lstars = np.log10(galaxy['Ms'].value)
    starsz = galaxy['Mzstar']/galaxy['Moxymade']
    ismz = galaxy['Mzism']/galaxy['Moxymade']

    z = np.zeros( lstars.size ) ###0.0*lstars
    starbar = ax.fill_between(lstars,z,starsz,facecolor='#d73027',edgecolor='#d73027', zorder=10, label="stars")   # stars, red
    ismbar = ax.fill_between(lstars,starsz,starsz+ismz,facecolor='#4575b4',edgecolor='#4575b4', zorder=10, label="ISM gas")  ## ISM, blue

    plt.plot([8.,11.5],[1,1], linestyle='dashed', color='black', linewidth=2)

    #lg = ax.legend(loc='best')
    #lg.draw_frame(False)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    # plt.xlim(8,11.5)
    plt.xlim(8.5,11.5)
    # plt.ylim(0,0.05)
    plt.ylim(0,1.18)
    plt.xlabel(r'log M$_{\star}$/M$_{\odot}$',fontsize=25)
    plt.ylabel(r'fraction of available metals, < 150kpc',fontsize=25)

    plt.tight_layout()
    plt.savefig('metal-census.png')


def metals_retained():
    galaxy = make_galaxy()
    lMstar = np.log10(galaxy['Ms'].value)
    Fgyear = "PS11"
    Fg = get_Fg(lMstar, Fgyear)
    galaxy['Mg'] = Fg * galaxy['Ms']
    galaxy['Mzism'] = galaxy['Zism'] * galaxy['Mg']

    ## now plot stuff :-)
    plot_metals_retained(galaxy)

if __name__ == "__main__":
# def main():
    ## anything that's defined in here is a global variable (!!!)

#    args = parse_args()
    # ps11()  ## have to pass args into beetle to not break
    metals_retained()
    sys.exit("~~~*~*~*~*~*~all done!!!! galaxies are fun!")
