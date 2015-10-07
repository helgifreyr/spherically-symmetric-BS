from scipy import *

def max_phi(rs,phi):
    maxPos = rs[argmax(phi)]
    maxVal = max(phi)
    return maxPos, maxVal

def energy(rs,m,sigma,phi,psi,w):
    H = 1.0 - 2.0*m/rs
    T44=-(psi**2*H) + (-1 - w**2/(H*sigma**2))*phi**2
    Ttot=2.*phi**2 - (4.*w**2*phi**2)/(H*sigma**2)

    intE = trapz(-rs**2 * sigma * Ttot,x=rs)
    return intE

def charge(rs,m,sigma,phi,psi,w):
    intQ = trapz(rs**2 * ( 2*phi**2/sigma ),x=rs)
    Q = w * intQ
    return Q

def light_rings(rs,m,sigma):
    Veff = ( rs - 2*m ) * sigma**2 / rs**3
    maxPos = argmax(Veff)
    minPos = argmin(Veff)
    if maxPos == 0:
        maxOut = 'No maxima'
    else:
        maxOut = [rs[maxPos],max(Veff)]
    if minPos == len(rs)-1:
        minOut = 'No minima'
    else:
        minOut = [rs[minPos],min(Veff)]
    return minOut,maxOut

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

def plot_functions(rs,m,sigma,phi,psi):
    plot(rs,phi)
    savefig('phi.pdf')
    clf()
    plot(rs,psi)
    savefig('psi.pdf')
