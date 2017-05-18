#! /usr/bin/env python

'''
AUTHOR: MOLLY PEEPLES
DATE: 11/29/2016
NAME: beetle.py
DESCRIPTION:
'''

import dill
import sys
import argparse

from read_in_galaxies import *
from create_galaxies import *
from calculate_rates import *
from make_metals import *
from make_plots import *

def parse_args():
    '''
    Parse command line arguments.  Returns args object.
    '''
    parser = argparse.ArgumentParser(description="a chemical / galaxy evolution code")
    parser.add_argument('mode', metavar='mode', type=str, action='store',
                        help='What mode to run in. options: create, read')

    parser.add_argument('--Nspline', dest='Nspline', action='store',
                        default=10,
                        help='Number points to fit spline over when calculating rates')

    parser.add_argument('--index', dest='index', type=int, action='store',
                        default=5,
                        help='Index for when galaxies file has more than one galaxy in it')

    parser.add_argument('--lifetime', dest='lifetime', type=float, action='store',
                        default=3.e9,
                        help='Lifetime for galaxy being created, in years')

    # parser.error("blah")

    args = parser.parse_args()
    return args


def beetle():
    """
    This is the main thing that beetle does
    """

    # do we create or read in galaxies?
    if args.mode == "create":
        print "let's create a galaxy!"
        galaxy = create_galaxies(lifetime=args.lifetime)
        outfile_base = "equilibrium"
    elif args.mode == "read":
        print "let's read in a galaxy!"
        filename = "sfh/test_sfh_gas_history_from_peter_and_gergo.dat"
        # filename = "galaxiesq.pkl"
        index = args.index
        galaxy = read_in_galaxy(filename, index)
        outfile_base = "galaxiesq_"+str(index)
    elif args.mode == "scatter":
        print "let's make some galaxies with scattered SFHs!"
        filename = "test_sfh_format.dat"
        galaxy = make_scattered_galaxy(filename)
        outfile_base = "scattered"
    elif args.mode == "compile":
        print "I'm going to compile your galaxies!"
        index = args.index  # this is going to be a maximum
        galaxy = read_in_galaxy(filename, index)
        outfile_base = "galaxiesq_"+str(index)
    elif args.mode == "plots":
        print "I'm just going to make some plots!"
        print "but I don't know how to read in a file!"
    else:
        print "****I don't know what to do here."
        print "****mode set to", args.mode
        exit

    outfile = outfile_base+"_basic.pkl"
    dill.dump(galaxy, open(outfile, "w"))
    print "basic galaxy printed to", outfile

    gasfile = "../Popping/SHAM+gas/Molly_sfh_gas_history.dat.gz"
    recyfile = "from_peter_with_recy.dat.gz"
    outfile = "galaxiesq.pkl"
    Nm = 16
    Niter = 1000
    # combine_galaxy_files(gasfile, recyfile, outfile, Nm=Nm, Niter=Niter)

    ## actually just want one redshift

    # calculate derivatives
    galaxy = calculate_rates(galaxy, Nspline=args.Nspline)
    galaxy = calculate_return(galaxy)
    print galaxy

    # metals
    galaxy = make_metals(galaxy, Nspline=args.Nspline)

    ## calculate ism metallicity, calculate rate of oxy production, calculate rate oxy leaving to get Zism
    ## if no metal loading that gives etaw and can back out etaa
    galaxy = calculate_etas(galaxy)

    # output galaxy?
    outfile = outfile_base+"_with_rates.pkl"
    dill.dump(galaxy, open(outfile, "w"))
    print "galaxy with everything printed to", outfile

    # plot some stuff!
    make_plots(galaxy, outfile_base)

if __name__ == "__main__":
# def main():
    ## anything that's defined in here is a global variable (!!!)

    args = parse_args()
    beetle()  ## have to pass args into beetle to not break
    sys.exit("~~~*~*~*~*~*~all done!!!! galaxies are fun!")
