from scipy import *
from scipy.integrate import odeint
from scipy.optimize import fsolve
from pylab import *

# w''(t) = 3/2 * w^2
# w(0) = 4, w(1) = 1

def function(U,y):
    # u1 = w', u2 = w'
    u1,u2 = U
    du1dy = u2
    du2dy = 3.0/2.0 * u1**2
    return [du1dy, du2dy]

u1_0 = 4
u1_1 = 1
dspan = linspace(0,1)

def objective(u2_0):
    U = odeint(function, [u1_0,u2_0],dspan)
    u1 = U[:,0]
    return u1[-1]-u1_1

u2_x = linspace(-100,0)
plot(u2_x,[objective(u2) for u2 in u2_x])
plot((-100, 0), (0, 0), 'k-')
show()
clf()

u2_0 = fsolve(objective,-5.0)

U = odeint(function, [u1_0, u2_0],dspan)

print u2_0

plot(dspan, U[:,0])
plot([1],[1],'ro')
show()
