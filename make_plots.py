from astropy import units as u
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['font.family'] = 'stixgeneral'
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors


from astropy.cosmology import WMAP9 as cosmo

def plot_masses(galaxy,output_name):
    fig = plt.figure(figsize=(11,8),dpi=300)
    ax = fig.add_subplot(111)
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['Mg'],color='blue',lw=3,label="all gas, inc. He")
    if "Mhi" in galaxy.colnames:
        ax.plot(galaxy['age'].to(u.Gyr),galaxy['Mhi'],color='green',lw=3,label="HI")
    if "Mh2" in galaxy.colnames:
        ax.plot(galaxy['age'].to(u.Gyr),galaxy['Mh2'],color='cyan',lw=3,label=r"H$_2$")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['Ms'],color='red',lw=3,label=r"stars")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['Mh']*cosmo.Ob0/cosmo.Om0,color='black',lw=3,label=r"halo baryons")
    if "Mcgm" in galaxy.colnames:
        ax.plot(galaxy['age'].to(u.Gyr),galaxy['Mcgm'],color='purple',lw=3,label=r"CGM")

    ax.set_yscale('log')
    plt.xlabel("time [Gyr]",fontsize=16)
    plt.ylabel(r"mass [M$_{\odot}$]",fontsize=16)
    lg = plt.legend(loc='best')
    lg.draw_frame(False)
    plt.tight_layout()

    plt.savefig(output_name+"_masses.png")
    plt.close(fig)




#-----------------------------------------------------------------------------------------------------

def plot_rates(galaxy,output_name):
    fig = plt.figure(figsize=(11,8),dpi=300)
    ax = fig.add_subplot(111)
    if "dMhdt" in galaxy.colnames:
        ax.plot(galaxy['age'].to(u.Gyr),galaxy['dMhdt']*cosmo.Ob0,color='black',lw=3,label=r"halo baryons")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['sfr'],color='red',lw=3,label=r"star formation")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['dMgdt'],color='blue',lw=3,label=r"all gas, inc. He")
    if "dMhidt" in galaxy.colnames:
        ax.plot(galaxy['age'].to(u.Gyr),galaxy['dMhidt'],color='green',lw=3,label=r"HI")
    if "dMh2dt" in galaxy.colnames:
        ax.plot(galaxy['age'].to(u.Gyr),galaxy['dMh2dt'],color='purple',lw=3,label=r"H$_2$")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['dMsdt'],color='orange',lw=3,label=r"stellar mass")

    plt.xlabel("time [Gyr]",fontsize=16)
    plt.ylabel(r"dM/dt [M$_{\odot}$ yr$^{-1}$]",fontsize=16)
    # plt.ylim(-6,6)
    lg = plt.legend(loc='best')
    lg.draw_frame(False)
    plt.tight_layout()

    plt.savefig(output_name+"_rates.png")
    plt.close(fig)




#-----------------------------------------------------------------------------------------------------

def plot_compare_recycle(galaxy,output_name):
    fig = plt.figure(figsize=(11,8),dpi=300)
    ax = fig.add_subplot(111)
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['sfr'] - galaxy['dMsdt'],color='red',lw=3,label=r"$\dot{M}_{\star}$-SFR")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['dMrdt'],color='green',lw=3,label=r"recycle from dfml/dt")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['dMrdt_expl'],color='limegreen',lw=3,label=r"recycle from fml, explicit")

    plt.xlabel("time [Gyr]",fontsize=16)
    plt.ylabel(r"dM/dt [M$_{\odot}$ yr$^{-1}$]",fontsize=16)
    # plt.ylim(-6,6)
    lg = plt.legend(loc='best')
    lg.draw_frame(False)
    plt.tight_layout()

    plt.savefig(output_name+"_compare_recycle.png")
    plt.close(fig)



#-----------------------------------------------------------------------------------------------------

def plot_recycle(galaxy,output_name):
    fig = plt.figure(figsize=(11,8),dpi=300)
    ax = fig.add_subplot(111)
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['Mr'],color='green',lw=3,label=r"recycle from fml")

    plt.xlabel("time [Gyr]",fontsize=16)
    plt.ylabel(r"dM/dt [M$_{\odot}$ yr$^{-1}$]",fontsize=16)
    # plt.ylim(-6,6)
    lg = plt.legend(loc='best')
    lg.draw_frame(False)
    plt.tight_layout()

    plt.savefig(output_name+"_recycle.png")
    plt.close(fig)




#-----------------------------------------------------------------------------------------------------

def plot_metal_mass(galaxy,output_name):
    fig = plt.figure(figsize=(11,8),dpi=300)
    ax = fig.add_subplot(111)
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['Moxymade'],color='black',lw=3,label=r"M$_{oxy}$ produced")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['Moxyism'],color='blue',lw=3,label=r"M$_{oxy}$ in the ISM")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['Mzstar'],color='red',lw=3,label=r"M$_{oxy}$ in stars")
    if 'Mzcgm' in galaxy.columns:
        ax.plot(galaxy['age'].to(u.Gyr),galaxy['Mzcgm'],color='green',lw=3,label=r"M$_{oxy}$ in CGM")
    # ax.plot(galaxy['Mh']),galaxy['Moxymade']/galaxy['Moxyism']), color='black', lw=3, label=r'M$_{oxy}$ made / ISM')
    # ax.plot(galaxy['age'].to(u.Gyr),galaxy['Moxyism']/galaxy['Moxymade'], color='black', lw=3, label=r'M$_{oxy}$ ISM / made')

    ax.set_yscale('log')
    plt.xlabel("time [Gyr]",fontsize=16)
    plt.ylabel(r"mass [M$_{\odot}$]",fontsize=16)
    # plt.ylim(-6,6)
    lg = plt.legend(loc='best')
    lg.draw_frame(False)
    plt.tight_layout()

    plt.savefig(output_name+"_metal_mass.png")
    plt.close(fig)


#-----------------------------------------------------------------------------------------------------

def plot_metallicity(galaxy,output_name):
    fig = plt.figure(figsize=(11,8),dpi=300)
    ax = fig.add_subplot(111)
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['Zcgm'],color='purple',lw=3,label=r"Z$_{cgm}$")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['Zism'],color='blue',lw=3,label=r"Z$_{ism}$")
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['Zstar'],color='red',lw=3,label=r"Z$_{stars}$")

    # ax.set_yscale('log')
    plt.xlabel("time [Gyr]",fontsize=16)
    plt.ylabel(r"metallicity",fontsize=16)
    # plt.ylim(0.005,0.01)
    lg = plt.legend(loc='best')
    lg.draw_frame(False)
    plt.tight_layout()

    plt.savefig(output_name+"_metallicity.png")
    plt.close(fig)



#-----------------------------------------------------------------------------------------------------

def plot_etas(galaxy,output_name):
    fig = plt.figure(figsize=(11,8),dpi=300)
    ax = fig.add_subplot(111)
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['etaa'], color='green', lw=3, label=r'$\eta_{\rm a} = \dot{M}_{\rm acc}$/SFR')
    ax.plot(galaxy['age'].to(u.Gyr),galaxy['etaw'], color='orange', lw=3, label=r'$\eta_{\rm w} = \dot{M}_{\rm wind}$/SFR')

    plt.xlabel("time [Gyr]",fontsize=16)
    plt.ylabel(r"$\eta\equiv\dot{\rm M}$/SFR",fontsize=16)
    plt.ylim(-3,10)
    lg = plt.legend(loc='best')
    lg.draw_frame(False)
    plt.tight_layout()

    plt.savefig(output_name+"_etas.png")
    plt.close(fig)



#-----------------------------------------------------------------------------------------------------

def make_plots(galaxy, output_name):
    print "making some plots..."
    print """----make_plots assumes that galaxy masses are in Msun.
----needs to be updated to explicitly check or convert"""

    plot_masses(galaxy, output_name)
    plot_rates(galaxy, output_name)
    # plot_compare_recycle(galaxy, output_name)
    # plot_recycle(galaxy, output_name)
    plot_metal_mass(galaxy, output_name)
    plot_metallicity(galaxy, output_name)
    plot_etas(galaxy, output_name)
