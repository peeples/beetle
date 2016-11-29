import numpy as np
from math import *
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from matplotlib import colors

t, z, mstar, frac, sfr = np.loadtxt('/Users/molly/research/metals/fmr/leitner/out/avgsfh-M10.70N1kfpdt0.00IMFc03.dat', usecols=(0,1,2,3,4), unpack=True)

lmstar = np.log10(mstar)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(z, lmstar, lw=4, label='Leitner', color='black')
ax.plot(z, 10.7-0.045*z-0.13*z*z, lw=3, label='van Dokkum', color='blue')
ax.legend(loc='lower left')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(0,1.5)
plt.ylim(9,11)
plt.xlabel(r'z',fontsize=20)
plt.ylabel(r'log Mstar',fontsize=20)
#plt.tight_layout()

plt.savefig('lmstar.eps')

