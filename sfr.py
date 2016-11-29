import numpy as np
from math import *
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from matplotlib import colors

t, z, mstar, frac, sfr = np.loadtxt('/Users/molly/Dropbox/fmr_stars/leitner/out/avgsfh-M10.70N1kfpdt0.00IMFc03.dat', usecols=(0,1,2,3,4), unpack=True)

lmstar = np.log10(mstar)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(z, sfr, lw=4, label='Leitner', color='black')
ax.plot(z, np.power(10.0,0.26+0.92*z-0.23*z*z)-1, lw=3, label='van Dokkum', color='blue')
ax.legend(loc='upper right')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(0,3)
plt.ylim(0,20)
plt.xlabel(r'z',fontsize=20)
plt.ylabel(r'SFR',fontsize=20)
#plt.tight_layout()

plt.savefig('sfr.png')
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t, sfr, lw=4, label='Leitner', color='black')
ax.legend(loc='upper right')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
# plt.xlim(0,3)
plt.ylim(0,20)
plt.xlabel(r't',fontsize=20)
plt.ylabel(r'SFR',fontsize=20)
#plt.tight_layout()

plt.savefig('sfr_t.png')

