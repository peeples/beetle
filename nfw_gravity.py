from astropy.cosmology import WMAP9 as cosmo
from astropy.constants import G
import astropy.units as u

import numpy as np

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['font.family'] = 'stixgeneral'
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors


c = 10  ## NFW concentration parameter
z = 0.0 ## redshift of interest

def Mencl(r):
    ## NFW!

    ## start with point source
    M = 1.e10 * u.Msun
    return M

def gravity(r):
    g = -1.0*G*Mencl(r) / (r*r)

    return g

def radius(r0,v0,g,t):
    return (r0 + v0*t + g*t*t).to(u.kpc)

def plot_nfw():
    ## want radius versus time for different launch velocities
    ## eventually: include pressure and hydrostatic drag 

    ## actually first start as if g is constant for all things

    ## r = np.arange(10,400,1) * u.kpc
    ## g = gravity(r)
    ## for now gravity is constant
    r0 = 10 * u.kpc
    g = gravity(r0)
    t = np.arange(0,1000,1) * u.Myr

    fig = plt.figure(figsize=(11,8),dpi=300)
    ax = fig.add_subplot(111)

    v0 = 100 * u.km / u.s
    ax.plot(t,radius(r0,v0,g,t),color='red',lw=3)
    v0 = 200 * u.km / u.s
    ax.plot(t,radius(r0,v0,g,t),color='orange',lw=3)
    v0 = 300 * u.km / u.s
    ax.plot(t,radius(r0,v0,g,t),color='green',lw=3)
    v0 = 400 * u.km / u.s
    ax.plot(t,radius(r0,v0,g,t),color='blue',lw=3)
    v0 = 500 * u.km / u.s
    ax.plot(t,radius(r0,v0,g,t),color='purple',lw=3)

    plt.xlim(0,1000)
    plt.ylim(0,300)
    plt.xlabel("time [Myr]",fontsize=20)
    plt.ylabel("radius [kpc]",fontsize=20)
    plt.tight_layout()
    plt.savefig("test.png")


if __name__ == "__main__":
    ## anything that's defined in here is a global variable (!!!)
    plot_nfw()
