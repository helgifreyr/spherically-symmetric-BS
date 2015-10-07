from scipy import *
from scipy.integrate import odeint
from scipy.optimize import fsolve
from scipy.integrate import trapz
from pylab import *

from functions import *
from equations import *

def read_parameters():
    pars = genfromtxt('parameters')
    alpha = pars[0]
    mu = pars[1]
    sigma0 = pars[2]
    w = pars[3]
    a_0 = pars[4]
    r1 = pars[5]
    rmax = pars[6]
    return alpha,mu,sigma0,w,a_0,r1,rmax

def solve(pars):
    alpha,mu,sigma0,w,a_0,r1,rmax = pars
    # changing how I move through the r2s matters depending on where I am in the
    # spiral. at the center of the spiral, 5,10,11,... worked fine. in the weak
    # field regime, the starting point has to be further out it seems
    # use test.py to play around with it
    # r2s = arange(28,rmax+1)
    # r2s = arange(1,rmax+1)
    r2s = concatenate((array([5.0,10.0]),arange(11,rmax+1)))

    for r2 in r2s:
        a_0 = fsolve(odeShooting,a_0,args=(pars,r2))
        # print r2,a_0[0]
    print 'Shooting parameter: ',a_0[0]

    rs,u1,u2,u3,u4 = odeShooting(a_0,pars,rmax,printStuff=1)
    funcsOut = open('functions.dat','w')
    for r,m,sigma,phi,psi in zip(rs,u1,u2,u3,u4):
        r_out = '%8.6f'%r
        m_out = '%14.8f'%m
        sigma_out = '%14.8f'%sigma
        phi_out = '%14.8f'%phi
        psi_out = '%14.8f'%psi
        funcsOut.write(r_out+m_out+sigma_out+phi_out+psi_out+'\n')

    parsOut = open('pars.dat','w')
    parsOut.write(str(alpha)+'  ')
    parsOut.write(str(mu)+'  ')
    parsOut.write(str(sigma0)+'  ')
    parsOut.write('%10.9f'%(w/u2[-1])+'  ')
    parsOut.write('%6.4f'%w+'  ')
    parsOut.write('%10.9f'%a_0+'  ')
    parsOut.write('%10.9f'%u1[-1])
    print 'Initial frequency: ',w
    print 'Real frequency: ', w/u2[-1]
    print 'Mass from asymptotics: ', u1[-1]
    
    return rs,u1,u2/u2[-1],u3,u4

def plot_functions(rs,m,sigma,phi,psi):
    plot(rs,phi)
    savefig('phi.pdf')
    clf()
    plot(rs,psi)
    savefig('psi.pdf')
