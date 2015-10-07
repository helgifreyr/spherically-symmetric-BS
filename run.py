from functions import *
from equations import *
from shooting import *

pars = read_parameters()
rs, m, sigma, phi, psi = solve(pars)

sigma=sigma/sigma[-1]
old_w = pars[3]
w = old_w/sigma[-1]

intE = energy(rs,m,sigma,phi,psi,w)
Q = charge(rs,m,sigma,phi,psi,w)

# print max_phi(rs,phi)
minLR = light_rings(rs,m,sigma)[0]
maxLR = light_rings(rs,m,sigma)[1]
print 'Light rings: ', minLR,', ', maxLR


intData = open('integral-data.dat','w')
intData.write('%4.2f'%old_w+'%12.9f'%w+'%12.9f'%Q+'%12.9f'%intE)

plot_functions(rs,m,sigma,phi,psi)
