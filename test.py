from scipy import *
from scipy.integrate import odeint
from scipy.optimize import fsolve
from scipy.integrate import trapz
from pylab import *

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

def function(U,r):
    u1,u2,u3,u4 = U
    # u1 = m, u2 = sigma, u3 = phi, u4 = psi
    H = 1.0-2.0*u1/r
    du1dr = r*r * alpha*alpha * ( mu*mu * u3*u3 + w*w * u3*u3/( H * u2*u2 ) + H * u4*u4 )
    dH = -2.0*du1dr/r + 2.0*u1/(r*r)
    du2dr = 2.0*r*alpha*alpha*u2 * ( u4*u4 + w*w * u3*u3/( H*H * u2*u2 ) )
    du3dr = u4
    du4dr = u3 * ( mu*mu/H - w*w/( H*H * u2*u2 ) ) - u4 * (2.0/r + dH/H + du2dr/u2 )
    return [du1dr, du2dr, du3dr, du4dr]

def odeShooting(a,printStuff=0):
    # u1 = m, u2 = sigma, u3 = phi, u4 = psi
    a = a[0]

    phi2 = a*( mu*mu - w*w/( sigma0*sigma0 ) )/6.0
    m3 = a*a * alpha*alpha * ( mu*mu + w*w/( sigma0*sigma0 ) )/3.0
    sigma2 = a*a * w*w * alpha*alpha/sigma0

    # u1 = m, u2 = sigma, u3 = phi, u4 = psi
    u1_0 = m3 * r1**3
    u2_0 = sigma0 + sigma2*r1*r1
    u3_0 = a + phi2*r1*r1
    u3_1 = 0.0
    u4_0 = 2.0*phi2*r1
    # print u1_0, u2_0, u3_0, u4_0

    rs = arange(r1,r2,0.001)

    U = \
    odeint(function,[u1_0,u2_0,u3_0,u4_0],rs,atol=1e-12,rtol=1e-12,h0=1e-3)
    # print U
    u1 = U[:,0]
    u2 = U[:,1]
    u3 = U[:,2]
    u4 = U[:,3]
    if printStuff==0:
        return u3[-1]-u3_1
    else:
        return rs,u1,u2,u3,u4

alpha,mu,sigma0,w,a_0,r1,rmax = read_parameters()
r2s = concatenate((array([5.0,10.0]),arange(11,rmax+1)))

r2 = 20.0
print a_0
a=fsolve(odeShooting,a_0)
print a
print odeShooting(a)
r2 = 21.0
a=fsolve(odeShooting,a)
print a

# for r2 in r2s:
#     a_0 = fsolve(odeShooting,a_0)
#     print r2,a_0[0]
#
# rs,u1,u2,u3,u4 = odeShooting(a_0,printStuff=1)
