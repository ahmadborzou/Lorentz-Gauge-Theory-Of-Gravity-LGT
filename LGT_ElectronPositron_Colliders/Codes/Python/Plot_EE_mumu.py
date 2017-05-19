# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 11:34:14 2017

@author: Ahmad Borzou

@email: ahmad_borzou@baylor.edu

@affiliation: Baylor University and Isfahan University of Technology

@more info: https://github.com/ahmadborzou/Lorentz-Gauge-Theory-Of-Gravity-LGT/wiki

This code needs numpy and matplotlib libraries in addition to python.
The whole package can be installed through ANAconda at:
https://www.continuum.io/downloads

"""

import numpy as np
from scipy.integrate import quad
from matplotlib import pyplot as plt

# angular part of LGT Cross Section
def LGTXS(theta):
    return np.cos(theta/2)**4

# angular part of QED Cross Section
def QEDXS(theta):
    return 1 + np.cos(theta)**2


# Total Cross Section
def TotalXS(theta,Ecm,r):
    #input in GeV
    #output in femtobarn(fb)
    # r=g/e ratio of coupling constants
    #Ecm is center of mass energy in GeV
    #0.39e12 is to return in fb
    a=1/137.04 # alpha in QED
    b=a*r*np.sqrt(12.-18.*r**2) # beta defined in the paper
    return 0.39e12/(4*Ecm**2)*(
                   a**2*QEDXS(theta)-b**2*LGTXS(theta)
                 )
# Asymmetry in number of events in Left-Right hemispheres of detector
def LREventAsym(Ecm,r,Lumin):
    def integrand(theta,Ecm,r):
        return TotalXS(theta,Ecm,r)*np.sin(theta)
    # numerical integration 
    DSigTot=quad(integrand,np.pi/2.,np.pi,args=(Ecm,r))[0] - quad(integrand,0.,np.pi/2.,args=(Ecm,r))[0]
    return Lumin*2*np.pi*DSigTot
                  

# length delta in the paper
def ParticleLength(r):
    # return delta, particle size, in units of Planck length
    return np.sqrt(40.*np.pi/(0.3*r))
    



# break the theta into 1000 pieces
theta = np.atleast_2d(np.linspace(0, np.pi, 1000)).T


# plot asymmetry in # events between left and right hemisphere
#set luminosity in fb^-1
Lumin=2.#/fb
Ecm=50.#GeV

# start plotting 
fig3=plt.figure()

# different values for r
r=np.linspace(0.,0.1,1000,retstep=False)
# compute left right asymmetry for different r
# and record it in a container
LRNAs=np.empty(np.size(r))
for ir in range(0,np.size(r)):
    LRNAs[ir]=LREventAsym(Ecm,r[ir],Lumin)

# plot    
plt.plot(r,LRNAs,color='b', ls='-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\frac{g}{e}$',fontsize=20)
plt.ylabel(r'$N(\theta:\frac{\pi}{2}-\pi)-N(\theta:0-\frac{\pi}{2})$',fontsize=15)
plt.text(0.0002, 5, r'$\sqrt{s} = 50$ GeV', fontsize=10)
plt.text(0.0002, 50, r'Lum.$= 2$ fb$^{-1}$', fontsize=10)

x1=fig3.add_subplot(111)
x2 = x1.twiny()
x2.set_xscale('log')
x2.set_xlim(x1.get_xlim())

x2.set_xticks(x1.get_xticks())
#x2.set_xticklabels(['1','e','d','f'])
x2.set_xbound(x1.get_xbound())
UpXValue=[]
for x in x1.get_xticks():
    UpXValue.append(int(ParticleLength(x)))
    
UpXValue[UpXValue==np.inf]=1.e10

x2.set_xticklabels(UpXValue)

x2.set_xlabel(r'$\frac{\delta}{l_{\mathrm{Pl}}}$',fontsize=20)


plt.plot()
plt.tight_layout()
plt.show()
fig3.savefig('N_AntiSymmetry.pdf')
