# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 16:29:14 2017

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



# Total Cross Section
def TotalXS(theta,Ecm,r):
    #input in GeV
    #output in femtobarn(fb)
    # r=g/e
    #Ecm is center of mass energy in GeV
    #0.39e12 is to return in fb
    a=1/137.04 # QED alpha. This cross section has nothing to do with QED. We only define g=r*e where e is the unit electric charge
    return 0.39e12*r**4*a**2/(16.*Ecm**2)*(
                   (18.*(1+np.sin(2.*theta))-5.*np.cos(2.*theta)-4.*np.cos(4.*theta)-np.cos(6.*theta)+np.sin(4.*theta))/(4.*np.sin(theta)**2)
                 )

# number of events:
# Lumin*\int \frac{d\sigma}{d|Omega}d\Omega
def NEvent(Ecm,r,Lumin):
    def integrand(theta,Ecm,r):
        return TotalXS(theta,Ecm,r)*np.sin(theta)
    DSigTot=quad(integrand,np.pi/18.,np.pi/2.,args=(Ecm,r))[0] 
    return Lumin*2.*np.pi*DSigTot
                  

# length delta in the paper
def ParticleLength(r):
    # return delta, particle size, in units of Planck length
    return np.sqrt(40.*np.pi/(0.3*r))
    


# break the theta into 1000 pieces
theta = np.atleast_2d(np.linspace(0, np.pi, 1000)).T





#set luminosity in fb^-1
Lumin=2.#/fb
Ecm=50.#GeV

# prepare figure
fig3=plt.figure()

# different values for r
r=np.linspace(0.,0.01,100,retstep=False)

# a container to hold number of events for each r
NEvs=np.empty(np.size(r))

for ir in range(0,np.size(r)):
    NEvs[ir]=NEvent(Ecm,r[ir],Lumin)

# plot    
plt.plot(r,NEvs,color='b', ls='-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\frac{g}{e}$',fontsize=20)
plt.ylabel(r'$N$',fontsize=15)
plt.text(1.5e-4, 1e-5, r'$\sqrt{s} = 50$ GeV', fontsize=10)
plt.text(1.5e-4, 1e-4, r'Lum.$= 2$ fb$^{-1}$', fontsize=10)


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

x2.set_xlabel(r'$\frac{\delta}{l_{\mathrm{Pl}}}$',fontsize=20,labelpad=15)

plt.ylim(1e-13, 1e-2)
plt.plot()
plt.tight_layout()
plt.show()
fig3.savefig('NAnnihilatedLeptons.pdf')