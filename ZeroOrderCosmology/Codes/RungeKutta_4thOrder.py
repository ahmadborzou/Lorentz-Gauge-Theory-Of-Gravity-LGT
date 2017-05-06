# -*- coding: utf-8 -*-

"""
Created on Sat Nov 12 18:49:40 2016

@author: Ahmad Borzou

@email: ahmad_borzou@baylor.edu

@affiliation: Baylor University and Isfahan University of Technology

@more info: https://github.com/ahmadborzou/Lorentz-Gauge-Theory-Of-Gravity-LGT/wiki

This code needs numpy and matplotlib libraries in addition to python.
The whole package can be installed through ANAconda at:
https://www.continuum.io/downloads

"""

import numpy as np
import matplotlib.pyplot as plt

## fixing some parameters
# the scale factor at present time
a=1
# the hubble constant at present time
H=71.
# the deceleration parameter at present time
q=-1.
# Newton's gravitational constant. Since the critic density will be 
# fixed below, this G will cancel another one in the field equations 
# so its value is not important
G=1.
# radiation energy density to go into GR calculations
rho0_r=0.000083
# critic energy density
rho_critic_GR = 3.*71.*71./(8*3.14*G)
rho_critic = 3.*H*H/(8*3.14*G)
# matter density
Mdensity_GR = 0.3
Mdensity=0.3
rho0 = Mdensity * rho_critic

# the source in LGT
def Source(a,y):
    return 4.*3.14*G*rho0/(y*(a**3.))

# second derivative of y in LGT which is known from field equations
def F(a,y,yp):
    return -yp**2./y+2.*y/(a**2.)-yp/a+Source(a,y)


# find the dy and dyp from RungeKutta method
def RK4(delta,a,y,yp):
    h=delta
    ##F=ypp
    dyp1=h*F(a,y,yp)
    dy1=h*yp
    dyp2=h*F(a-h/2.,y-dy1/2.,yp-dyp1/2.)
    dy2=h*(yp-dyp1/2.)
    dyp3=h*F(a-h/2.,y-dy2/2.,yp-dyp2/2.)
    dy3=h*(yp-dyp2/2.)
    dyp4=h*F(a-h,y-dy3,yp-dyp3)
    dy4=h*(yp-dyp3)
    return (dy1+2.*dy2+2.*dy3+dy4)/6.,(dyp1+2.*dyp2+2.*dyp3+dyp4)/6.


# hubble parameter in GR directly from the known solution
def H_GR(a):
    H= np.sqrt(8.*3.14*G*rho_critic_GR/3.*(
                        Mdensity_GR/(a**3.)+
                        (1-Mdensity_GR)+
                        rho0_r/(a**4.)
                      )
                    )
    return H 

# \frac{d\dot{a}}{da}
def yp_GR(a):
    return (   G * rho_critic_GR* 1.45* (2.*a**4.* (1-Mdensity_GR) - a *Mdensity_GR - 2.* rho0_r)     )/(
                            a**4.*(( G*rho_critic_GR* ( a**4. *(1-Mdensity_GR) + a *Mdensity_GR + rho0_r) )/a**4.)**.5
                            )

# Rate of interaction between electron, neutrino, neutron, proton
def WeakRate(a):
    return 1.e-73/a**5.

# define the main function 
def main(amax,yin,ypin):
    
    # This is a recursive function and needs to be terminated as some point
    # This if condition will do the job
    if amax < 1.e-11:
        return True

    # some initial values
    Ninterval =1000000
    #Ninterval =10
    age=0.
    age_GR=0.
    amin = amax*0.0001 
    delta = (amax - amin)/Ninterval
    
    
    # arrays or containers
    a=np.linspace(amin, amax, num=Ninterval, retstep=False)
    a=a[::-1]
    y=np.empty(np.size(a))
    yp=np.empty(np.size(a))
    ypp=np.empty(np.size(a))
    age=np.empty(np.size(a))
    age_GR=np.empty(np.size(a))
    source=np.empty(np.size(a))
    
    
    
    print(" amax: %g amin: %g Ninterval: %d delta: %g" % (amax,amin,Ninterval,delta))
    print "a size: ",np.size(a)
    
    # a loop that goes back in time
    for ia in range(0,np.size(a)):
        # setting the initial values
        if ia==0:
            y[ia]=yin
            yp[ia]=ypin
            age[ia]=0.
            age_GR[ia]=0.
        else:
            # Runge-Kutta method
            dy,dyp=RK4(delta,a[ia-1],y[ia-1],yp[ia-1])
            #print dy,dyp
            y[ia]=y[ia-1]-dy
            yp[ia]=yp[ia-1]-dyp
            # needed for finding the age
            age[ia]=age[ia-1]+1000.*delta/y[ia]
            age_GR[ia]=age_GR[ia-1]+1000.*delta/(H_GR(a[ia])*a[ia])
        
        
        
        # the source in LGT
        source[ia] = Source(a[ia],y[ia])

        # second derivative of y in LGT which is known from field equations
        ypp[ia]=F(a[ia],y[ia],yp[ia])

        
        # print out some information
        if ia % 100000==0:
            print " a: " , a[ia], " y: " ,y[ia], " yp: " , yp[ia], " ypp: ",ypp[ia], " HRatio: ",y[ia]/(a[ia]*H_GR(a[ia])/4.65e43)
            
    
    
    ## Do the plotting here
    # \frac{\dot{a}}{a}
    fig0=plt.figure()
    if amin < 1e-11:
        plt.plot(a,y/(4.65e43*a),color='b', ls='-', label=r'$\mathrm{H}_{\mathrm{LGT}}$') # the number is to convert to GeV
        plt.plot(a,H_GR(a)/4.65e43,color='g', ls='-.',linewidth=3.0, markersize=5, label=r'$\mathrm{H}_{\mathrm{GR}}$')
        plt.plot(a,WeakRate(a),color='r', ls='--',linewidth=3.0, markersize=5, label=r'$\Gamma$')
        plt.ylabel(r'(GeV)',fontsize=12)
    else:
        plt.plot(a,y/(4.65e43*a),color='b', ls='-', label=r'LGT') # the number is to convert to GeV
        plt.plot(a,H_GR(a)/4.65e43,color='g', ls='-.',linewidth=3.0, markersize=5, label=r'GR')  
        plt.ylabel(r'H (GeV)',fontsize=12)
        
    plt.yscale('log')
    plt.xscale('log')
    plt.legend(loc='best',fontsize = 'large',shadow=True,fancybox=True, borderpad=2)
    plt.xlabel(r'$a$',fontsize=15)
    if amax > 0.1:
        plt.xlim(plt.xlim()[0],0.1)
    plt.plot()
    plt.tight_layout()
    fig0.savefig('adot_a_%d_%d.pdf'%(-np.log10(amin),-np.log10(amax)))
    
    # only in the first iteration
    if amax >0.1:
        fig1=plt.figure()
        plt.plot(a,y/(4.65e43*a),color='b', ls='-', label=r'LGT') # the number is to convert to GeV
        plt.plot(a,H_GR(a)/4.65e43,color='g', ls='-.',linewidth=3.0, markersize=5, label=r'GR')
        plt.yscale('log')
        plt.legend(loc='best',fontsize = 'large',shadow=True,fancybox=True, borderpad=2)
        plt.ylabel(r'H (GeV)',fontsize=12)
        plt.xlabel(r'$a$',fontsize=15)
        plt.xlim(0.01,1.)
        plt.ylim(1e-42,1e-36)
        plt.plot()
        plt.tight_layout()
        fig1.savefig('adot_a_001_1.pdf')
    
        # \frac{\dot{da}}{dt}
        fig2=plt.figure()
        plt.plot(a,y/4.65e43,color='b', ls='-', label=r'LGT')
        plt.plot(a,H_GR(a)*a/4.65e43,color='g', ls='-.',linewidth=3.0, markersize=5, label=r'GR')
        plt.yscale('log')
        #plt.xscale('log')
        plt.legend(loc='best',fontsize = 'large',shadow=True,fancybox=True, borderpad=2)
        plt.ylabel(r'$\dot{a}$ (GeV)',fontsize=15)
        plt.xlabel(r'$a$',fontsize=15)
        plt.plot()
        plt.tight_layout()
        fig2.savefig('da_dt_%d_%d.pdf'%(-np.log10(amin),-np.log10(amax)))
            
    # age
    age=age[-1]-age ## this means, in what we report in paper, the age of the universe from a=0 to 10^-4 is neglected
    age_GR=age_GR[-1]-age_GR
    fig3=plt.figure()
    plt.plot(a,age,color='b', ls='-', label=r'LGT')
    plt.plot(a,age_GR,color='g', ls='-.',linewidth=3.0, markersize=5, label=r'GR')
    plt.legend(loc='best',fontsize = 'large',shadow=True,fancybox=True, borderpad=2)
    plt.ylabel(r'age ($10^{9}$ year)',fontsize=10)
    plt.xlabel(r'$a$',fontsize=15)
    plt.plot()
    plt.tight_layout()
    fig3.savefig('age_%d_%d.pdf'%(-np.log10(amin),-np.log10(amax)))
    
    
    # continue going back in time with a smaller step-size 
    main(a[-1],y[-1],yp[-1])


# Finally run the code
main(a,H,-H*q)
