import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table, Column
from astropy import units as u

import sys
import argparse

from constants import *

def test_fml(galaxy):
    mass_made = 0
    mass_made_z0 = 0
    for i in range(len(galaxy)):
        ## calculate timestep
        DT = calculate_timestep(galaxy['age'].quantity, i, len(galaxy))

        ## how much mass is made in this timestep
        mass_made += galaxy['sfr'].quantity[i]*DT
        mass_made_z0 += galaxy['sfr'].quantity[i]*DT * (1 - fml(galaxy['lookback'].quantity[i]))
        
    print "mass_made = ",mass_made.to(u.solMass)," vs. mass_made_z0 = ",mass_made_z0.to(u.solMass)," vs. actual Mstar = ",galaxy['Ms'].quantity[len(galaxy)-1]


def test_table():
    g1 = Table()
    dummy = np.zeros(17)
    col_dummy = Column(name='age',data=dummy,unit=u.yr)
    g1.add_column(col_dummy)
    col_dummy = Column(name='bar',data=dummy)
    g1.add_column(col_dummy)
    col_dummy = Column(name='foo',data=dummy)
    g1.add_column(col_dummy)
    g1['blah'] = "this is blah"
    galaxies = {}
    for i in range(3):
        galaxies[i] = g1
    print galaxies[2]['age'][12]
    galaxy = galaxies[1]
    print galaxy
    print galaxies[1]['foo']
    print galaxies[0]['blah']

def test_return(animal):
    if animal == "elephant":
        print "animal is elephant!"
        return animal
    print "setting animal to giraffe"
    animal = "giraffe"
    return animal


class Galaxy:
    def __init__(self,sfr,mstar):
        self.sfr = sfr * u.Msun / u.yr
        self.mstar = mstar * u.Msun
    animal = "Galaxy does not have an animal"
    ID = "I don't have an ID"
    def ssfr(self):
        return self.sfr/self.mstar
    def fmr(self):
        m = np.log10(self.mstar.value) - 10
        lsfr = np.log10(self.sfr.value)
        return (8.90 + 0.37*m - 0.14*lsfr - 0.19*m*m + 0.12*m*lsfr - 0.054*lsfr*lsfr)
    elements = ["H", "C", "O", "Fe"]
        
def blah():
    g = Galaxy(1,3.e10)
    print g.fmr()
    g.animal = "elephant"
    print g.animal
    gg = Galaxy(2, 1.e10)

    galaxies = {1 : g, 2 : gg}  # this is a dictionary
    print galaxies
    print galaxies[1].sfr

def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--time_trend', dest='trend_with_time', action='store_true',
                        default=False,
                        help='Run time trending')

    parser.add_argument('--no_check', dest='check_for_new', action='store_false',
                        default=True,
                        help='Run in serial')

    parser.add_argument('--processors', dest='n_processors', action='store',
                        default=1,
                        help='Number of processors for making gainmaps')

    parser.add_argument('--regress', dest='regress', action='store_true',
                        default=False,
                        help='Run Regression set')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    print "hello world"

    for arg in sys.argv:
        print arg

    if "fart" in sys.argv:
        print "ptttth"
    

    args = parse_args()

    print args
    
