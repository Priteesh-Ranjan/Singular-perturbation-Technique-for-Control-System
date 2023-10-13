# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 10:29:05 2023

@author: Lenovo
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def fun(t, y, lmda):
    k1 = y[0]
    k2 = y[1]
    k3 = y[2]
    dk1dt = (3.6 * k2 + 180 * k2 ** 2 - 1)
    dk2dt = (1 / lmda) * (-18654 * k1 + 581 * k2 + (1.8 * k3) + (180 * k3) * k2)
    dk3dt = (1 / lmda) * (-lmda * 37308 * k2 + 1162 * k3 + 180 * k3 ** 2 - 4600)
    return [dk1dt, dk2dt, dk3dt]

y0 = [0.01, 0.01, 0.01]
tspan = np.linspace(0.01, 0.001, 1000)
lambda_values = [0.25, 0.5, 1, 2]

E1 = 12.7;
E2 =  0.004638913435227;

# Plot k1
plt.figure()
for lmda in lambda_values:
    sol = solve_ivp(fun, [tspan[0], tspan[-1]], y0, args=(lmda,), method='RK45', t_eval=tspan)
    plt.plot(sol.t, sol.y[0], label=f'lambda = {lmda}')
plt.xlabel('Time')
plt.ylabel('K1')
plt.title('Solution of K1 for different lambda values')
plt.legend()
plt.grid(True)
plt.show()

# Plot k2
plt.figure()
for lmda in lambda_values:
    sol = solve_ivp(fun, [tspan[0], tspan[-1]], y0, args=(lmda,), method='RK45', t_eval=tspan)
    plt.plot(sol.t, sol.y[1], label=f'lambda = {lmda}')
plt.plot(tspan, E1*sol.y[0]- E2, '--r', label='K2 = 17.3 * K1 - 0.0046')
plt.xlabel('Time')
plt.ylabel('K2')
plt.title('Solution of K2 for different lambda values')
plt.legend()
plt.grid(True)
plt.show()

# Plot k3
plt.figure()
for lmda in lambda_values:
    sol = solve_ivp(fun, [tspan[0], tspan[-1]], y0, args=(lmda,), method='RK45', t_eval=tspan)
    plt.plot(sol.t, sol.y[2], label=f'lambda = {lmda}')
plt.axhline(y=2.77, color='r', linestyle='--', label='k3 = 2.77')
plt.xlabel('Time')
plt.ylabel('K3')
plt.title('Solution of K3 for different lambda values')
plt.legend()
plt.grid(True)
plt.show()
