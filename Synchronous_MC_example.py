# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 10:40:32 2023

@author: Lenovo
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 06:16:22 2023

@author: Lenovo
"""
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from sympy.plotting import plot3d

# Define your symbol
x = sp.Symbol('x')
y = sp.Symbol('y')
s = sp.Symbol('s') 


# Values of all the parameters used
dlta = 0.512
b = 2.2143
a = 1.2143
s = sp.Symbol('s')
lmda = 2
c = 47.6190
n = 1
d = 0.072
ep = 0.1

# The Function V(x)
y1 = sp.asin((sp.sin(dlta) / (1 + x)) - dlta)
N = sp.cos(y1 + dlta) - sp.cos(dlta)
V = -a * x + b * N
Vx = sp.integrate(V, (x, 0, x))
print("Simplified expression for V(x): ", Vx)

## The function Wx
h1 = y1
M = sp.integrate((sp.sin(s + dlta) - sp.sin(h1 + dlta)), (s,0,x))
W = (lmda / 2) * ((y - h1) ** 2) + (n * c / lmda) * (1 + x) * M
print("\nSimplified expression for W(x,y):  ", W)

## The Function V2
v_dot = sp.diff(h1, x)
w_dot1 = sp.diff(W, y)
w_dot2 = sp.diff(W, x)
v = d*V + (1-d)*W
v_dot = d * sp.diff(h1, x) + ((1-d) / ep) * sp.diff(W, y) + d * sp.diff(W, x)
print("\nSimplified expression for v(x,y):  ", v)


# Getting the plot of the function which shows v(x) is positive definite function
p = plot3d(v, (x, -0.5, 0.5), (y, -10, 10), xlabel='x', ylabel='y', title='Plot of v(x, y)',dpi = 500)

# plot of derivative derivative of v
p = plot3d(v_dot, (x, -0.5, 0.5), (y, -10, 10), xlabel='x', ylabel='y', title='Plot of v_dot(x, y)',dpi = 500)


# Define the Lyapunov function
def lyapunov(x):
    M = 0.147;Eq = 1.22;theta = 0.38;k = 2.3;
    return (0.5 * M * x[1]**2) - (k * Eq * np.cos(x[0] + theta)) + (k *0.9*np.cos(theta)) - 0.9*x[0] ;

# Define the range of x and y values
x = np.linspace(-0.5, 1.1, 100)
y = np.linspace(-4, 4, 100)

# Create a grid of x and y values
X, Y = np.meshgrid(x, y)

# Initialize Z
Z = np.zeros_like(X)

# Calculate the values of the Lyapunov function
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        Z[i, j] = lyapunov([X[i, j], Y[i, j]])


# Plot the region of attraction
plt.contourf(Y,X,Z, levels=[0.0, 0.15, 0.3, 0.5], linewidths=2.5)

#plotting the bound on the plot
plt.plot([-3.65, 3.65],[-0.42, -0.42] , 'k-')
plt.plot([0, 0],[1.1, 0] , 'k-')
plt.plot([-4, 0],[0, 0] , 'k-')
plt.plot([4, 0],[0, 0] , 'k-')
plt.plot([0, 0],[0, -0.42] , 'k-')
plt.fill_between(y, -0.5, -0.42, where=(y >= -4) & (y <= 4), color='white')
plt.grid()
#plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Region of Attraction')
plt.show()



## Plotting the range of perurbation parameter upto which 
#Definig the Parameters
M = 147e-3; xd1 = 0.3; xe = 0.05; xd = 2; tau = 6.6; 
Pm = 0.185; Efd = 1.;
c = 1/M*(xe+xd1)
b = (xd - xd1)/(xe + xd1)
a = (xe + xd)/(xe + xd1)
theta = 0.3; eta = 0.62; zeta = 0.05;
del_val = 0.25

k = np.sin(del_val) / ((1 - theta) * np.sqrt((1 - theta)**2 - (np.sin(del_val))**2))

# imposing restriction on k
if k < a/b:
    k = a/b
k = np.sin(del_val) / ((1 - theta) * np.sqrt((1 - theta)**2 - (np.sin(del_val))**2))



num = c * (1 - theta) * eta * zeta
den = 0.5 * b * k + b * k * lmda * zeta + (1 / (d * 4 * (1 - d) * zeta)) * (b * (1 - d) * zeta + d * np.sqrt(k**2 + zeta**2 * (lmda * k + eta * c / lmda))**2)
eps = num / den
delx = 0.001
n1 = int(1 / delx + 1)
na = np.linspace(0, 1, n1)
nxa = np.linspace(0, 1, n1)
epsx_vala = np.zeros(n1)
n = 1 + c * (1 - theta) * eta
nx = 0
epsx_valb = np.zeros(n1)

for j in np.arange(0, 1 + delx, delx):
    epsx = (c * (1 - theta) * eta * zeta) / (0.5 * b * k + b * k * lmda * zeta + (1 / (j * 4 * (1 - j) * zeta)) * (b * (1 - j) * zeta + j * np.sqrt(k**2 + zeta**2 * (lmda * k + eta * c / lmda))**2))
    epsx_valb[nx] = epsx
    nx += 1

plt.figure(1)
plt.plot(nxa, epsx_valb, '-k', color='b', linewidth=2)
plt.ylabel('\u03B5 *(d)',fontsize=16)
plt.xlabel('d',fontsize=16)
plt.title('\u03B5 *(d) vs d',fontsize=16)
plt.legend(['\u03B5 corresponding to d'],fontsize=16)
plt.grid()
plt.show()
fig = p._backend.fig
fig.savefig("D:\Sem 2\M.Tech Guide Papers/_Ep(d)_vs(d).png", dpi=1000)

