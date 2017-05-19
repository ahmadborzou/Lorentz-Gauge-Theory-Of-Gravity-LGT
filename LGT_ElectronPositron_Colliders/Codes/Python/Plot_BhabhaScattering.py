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
from matplotlib import pyplot as plt



# the cross section of e+e->e+e
def TotalXS(theta,Ecm,r):
    #input in GeV
    #output in nano barn(nb)
    # r=g/e
    #Ecm is center of mass energy in GeV
    #0.39e6 is to return in nb
    a=1/137.04 # QED coupling
    return 0.39e6*a**2/(4.*Ecm**2)*(
                  2.*( (1.+np.cos(theta/2.)**4)/(np.sin(theta/2.)**4) + (1.+np.cos(theta)**2)/2. - 2.*(np.cos(theta/2.)**4)/(np.sin(theta/2.)**2)  )
                  -3.*r**2/64.*np.sin(theta)**8/np.sin(theta/2.)**12
                  +18.*r**4*np.cos(theta/2.)**4/(np.tan(theta/2.)**4)
                 )



#set luminosity in fb^-1
#center of mass energy in GeV
Lumin=2.#/fb
Ecm=50.#GeV
fig3=plt.figure()


# break the theta into 1000 pieces
theta=np.linspace(0,np.pi,20,retstep=False)

# a container for cross section
XS=np.empty(np.size(theta))

# 3 values for r
r=[0.,0.01,1.]
# colors, line style,mark,legend
col=['k','b','g']
linest=['-','none','-.']
mark=['','^','']
lb=[r'$r=0$ (QED)',r'$r=0.01$',r'$r=1$']

# plot
for ir in range(0,np.size(r)):
    for it in range(0,np.size(theta)):
        XS[it]=TotalXS(theta[it],Ecm,r[ir])
    plt.plot(np.cos(theta),XS,color=col[ir], linestyle=linest[ir],marker=mark[ir],label=lb[ir])
    #plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'$\cos(\theta)$',fontsize=15)
    plt.ylabel(r'$\frac{d\sigma}{d\Omega}(\mathrm{nb})$',fontsize=20)
    plt.text(0.1, 50, r'$\sqrt{s} = 50$ GeV', fontsize=10)
    plt.text(0.1, 20, r'Lum.$= 2$ fb$^{-1}$', fontsize=10)


plt.legend(loc='best',fontsize = 'large',shadow=True,fancybox=True)
plt.plot()
plt.tight_layout()
plt.show()
fig3.savefig('XSBhabha.pdf')