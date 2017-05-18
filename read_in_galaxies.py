#! /usr/bin/env python

'''
AUTHOR: MOLLY PEEPLES
DATE: 09/02/2015
NAME: read_in_galaxies.py
DESCRIPTION:
'''

import pickle
from astropy.table import QTable, Table, Column
from astropy.io import ascii
from astropy import units as u
import numpy as np
from astropy.cosmology import WMAP9 as cosmo

from constants import *


def read_in_galaxy(filename, index):
    print "attempting to read in",filename

    if ".pkl" in filename:
        print "assuming that",filename,"is a pickle file with only a dict of galaxies in it"
        galaxies = pickle.load(open(filename,"r"))
        print "for now we're only dealing with one galaxy at a time; taking galaxy #",index
        assert(index in galaxies.keys()), "index %i not in file :-(" % index
        galaxy = galaxies[index]
        col_z = Column(name='z',data=galaxies['z'])
        galaxy.add_column(col_z,0)
    elif "peter" in filename:
        print "assuming that",filename,"is an ascii file from peter behroozi and gergo popping with Mh,Ms,dMhdt,sfr,Mg,Mhi,Mh2"
        print "also assuming that",filename,"only has one iteration"
        z,Mh,Ms,dMhdt,sfr,Mg,Mhi,Mh2 = np.loadtxt(filename,usecols=(0,7*index+1,7*index+2,7*index+3,7*index+4,7*index+5,7*index+6,7*index+7),unpack=True)
        galaxy = QTable([z,Mh,Ms,dMhdt,sfr,Mg,Mhi,Mh2],names=('z','Mh','Ms','dMhdt','sfr','Mg','Mhi','Mh2'))

    ## need to add units
    galaxy['Mh'].unit = u.Msun
    galaxy['Ms'].unit = u.Msun
    galaxy['Mg'].unit = u.Msun
    galaxy['Mhi'].unit = u.Msun
    galaxy['Mh2'].unit = u.Msun
    galaxy['dMhdt'].unit = u.Msun / u.yr
    if 'dMrdt' in galaxy.colnames:
        galaxy['dMrdt'].unit = u.Msun / u.yr
    galaxy['sfr'].unit = u.Msun / u.yr

    ### Helium is not included in the gas masses!!
    galaxy['Mg'] = galaxy['Mg'] * HELIUM_CORR

    ### Let's make a CGM! ###
    galaxy['Mcgm'] = galaxy['Mh']*cosmo.Ob0/cosmo.Om0 - galaxy['Ms'] - galaxy['Mg']

    age = cosmo.age(galaxy['z'])
    col_age = Column(name='age', data=age)
    galaxy.add_column(col_age,1)

    tlb = cosmo.lookback_time(galaxy['z'])
    col_tlb = Column(name='lookback', data=tlb)
    galaxy.add_column(col_tlb,2)

    ## check to make sure sorted so earlier is earlier
    galaxy.sort('age')

    return galaxy


#-----------------------------------------------------------------------------------------------------


def read_in_galaxies(filename, **kwargs):
    print "attempting to read in", filename
    Nm = kwargs.get("Nm", 16)
    Niter = kwargs.get("Niter", 1000)
    index = kwargs.get("index", 0) ## index = which mass

    if ".pkl" in filename:
        print "assuming that", filename,"is a pickle file with only a dict of galaxies in it"
        galaxies = pickle.load(open(filename, "r"))
    # else:
    #     print "file isn't a pickle file I don't know what to do :-("
    #     continue

    iters = np.arange(Niter) + 1
    Nz = len(galaxies['z'])

    ## get the values
    Mh = [galaxies[k][zz]['Mh'] for zz in range(len(galaxies['z'])) for k in iters]  ## this is just one mass
    Mh.reshape((Niter, Nz))


    ## need to add units
    galaxy['Mh'].unit = u.Msun
    galaxy['Ms'].unit = u.Msun
    galaxy['Mg'].unit = u.Msun
    galaxy['Mhi'].unit = u.Msun
    galaxy['Mh2'].unit = u.Msun
    galaxy['dMhdt'].unit = u.Msun / u.yr
    if 'dMrdt' in galaxy.colnames:
        galaxy['dMrdt'].unit = u.Msun / u.yr
    galaxy['sfr'].unit = u.Msun / u.yr

    ### Helium is not included in the gas masses!!
    galaxy['Mg'] = galaxy['Mg'] * HELIUM_CORR



#-----------------------------------------------------------------------------------------------------


def combine_galaxy_files(gasfile, recyfile, outfile, **kwargs):
    # need to open both files and read in galaxy by galaxy
    # adopted from read_galaxies.py code snippet Gergo Popping sent me
    Nm = kwargs.get("Nm", 16)
    Niter = kwargs.get("Niter", 1000)


    #read in data, takes a while because np.genfromtxt is sloooowwwww
    print "reading in ",gasfile
    gasgal = np.genfromtxt(gasfile)
    Ngasparam = 7 #Number of parameters per galaxy for gas file
    len_realizations_gas = len(gasgal)/Niter #the length of a realization. i.e., the number of output redshifts


    print "reading in ",recyfile
    recygal = np.genfromtxt(recyfile)
    Nrecyparam = 5 #Number of parameters per galaxy for recy file
    len_realizations_recy = len(recygal)/Niter #the length of a realization. i.e., the number of output redshifts
    Nextra = len_realizations_recy - len_realizations_gas

    print "done reading data!"

    redshifts = gasgal[:,0]  ## gas file has fewer redshift outputs so these are the ones we want

    galaxies = {"z" : redshifts[0:len_realizations_gas]}
    for i in range(Niter):
        #Going through the 1000 realizations one by one
        for j in range(Nm):
            #Going through the 16 galaxies for each realization one by one
            Mh     = gasgal[i*len_realizations_gas:(i+1)*len_realizations_gas, j*Ngasparam+1] #Halo mass
            Sm     = gasgal[i*len_realizations_gas:(i+1)*len_realizations_gas, j*Ngasparam+2] #Steller mass
            dMh_dt = gasgal[i*len_realizations_gas:(i+1)*len_realizations_gas, j*Ngasparam+3] #Hallo accretion rate
            SFR    = gasgal[i*len_realizations_gas:(i+1)*len_realizations_gas, j*Ngasparam+4] #SFR
            Mgas   = gasgal[i*len_realizations_gas:(i+1)*len_realizations_gas, j*Ngasparam+5] #Cold gas mass (no Helium included)
            MHI    = gasgal[i*len_realizations_gas:(i+1)*len_realizations_gas, j*Ngasparam+6] #HI mass (no Helium included)
            MH2    = gasgal[i*len_realizations_gas:(i+1)*len_realizations_gas, j*Ngasparam+7] #H2 mass (no Helium included)

            Mh_z0 = Mh[-1] #z=0.0 halo mass
            """
            Time evolution of the above parameters is read in for one galaxy with z=0.0 halo mass Mh_z0. This can be stored in a new file.
            """
            dMrdt  = recygal[i*len_realizations_recy+Nextra-1:(i+1)*len_realizations_recy-1, j*Nrecyparam+5] # dMrecy/dt
            galaxy = QTable([Mh, Sm, dMh_dt, SFR, dMrdt, Mgas, MHI, MH2], names=('Mh','Ms','dMhdt','sfr','dMrdt','Mg','Mhi','Mh2'))
            galaxies[(i+1)*(j+1)] = galaxy

    ## save galaxies as a pickle file!
    pickle.dump(galaxies,open(outfile,"w"))

    ## output as an ascii file...
    print "I don't know how to output as an ascii file"




#-----------------------------------------------------------------------------------------------------

def make_scattered_galaxy(filename, **kwargs):
    print "attempting to read in ", filename

    galaxy = ascii.read(filename)
    galaxy.remove_column('i')
    galaxy.rename_column('SFR', 'sfr')

    galaxy['sfr'] = np.power(10.0,galaxy['sfr'])

    galaxy['Ms'].unit = u.Msun
    galaxy['sfr'].unit = u.Msun / u.yr


    return galaxy
