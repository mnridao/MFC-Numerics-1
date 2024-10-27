# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 19:48:13 2024

@author: guest997558
"""

class Parameters:
    
    def __init__(self):
        self.nu = 0.006 

# Base class.
class EulerianAdvection:
    
    def __init__(self):
        pass



import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def analytic(X, t, nu):
    
    return (np.exp(-(X-4*t)**2/(4*nu*(t+1))) + 
            np.exp(-(X-4*t-2*np.pi)**2/(4*nu*(t+1))))

def burg_system(u, t, k, nu):
    
    # Spatial derivative in the Fourier domain.
    u_hat = np.fft.fft(u)
    u_hat_x = 1j*k*u_hat
    u_hat_xx = -k**2*u_hat
    
    # Switching in the spatial domain.
    u_x = np.fft.ifft(u_hat_x)
    u_xx = np.fft.ifft(u_hat_xx)
    
    # ODE resolution.
    u_t = u*u_x + nu*u_xx
    return u_t.real


if __name__ == "__main__":
    
    # Define problem parameters.
    nu = 0.006
    # nu = 0.
    
    # Generate grid.
    x0, xL = 0, 10
    dx = 0.01
    nx = int((xL - x0)/dx)
    X = np.linspace(x0, xL, nx, False) # False for periodic bc.
    
    # Solver parameters.
    endtime = 10
    dt = 0.0001
    nt = int(endtime/dt)
        
    # Initial condition.
    phi0 = np.exp(-2 * (X - 0.5 * xL)**2)
    plt.plot(X, phi0)
    
    # Wavenumber discretisation. 
    k = 2*np.pi*np.fft.fftfreq(nx, d = dx)
    T = np.linspace(0, endtime, nt) #Temporal array
    U = odeint(burg_system, phi0, T, args=(k, nu,), mxstep=5000).T
    
    # Initialise matrix for trapezoidal scheme.
    A = np.zeros(shape=(nx, nx))
    np.fill_diagonal(A, 1+nu/dx)
    
    off = -nu*dt/(2*dx)
    for i in range(nx - 1):
        A[i, i+1] = off
        A[i+1, i] = off 
    
    A[0, -1] = off 
    A[-1, 0] = off
    
    Ainv = np.linalg.inv(A)
    
    #%%
        
    # Time loop.
    for i in range(15):
        
        # Initialise new values.
        phi = phi0.copy()
        # phiS = phi0.copy()
                
        # Loop over space.
        for j in range(nx):
                        
            # # Advective form upwind (not conservative).
            # phiS[j] += dt/dx *(-0.5*phi0[j]*(phi0[j]-phi0[(j-1)%nx]) 
            #                   + 0.5*nu*(phi0[(j+1)%nx]-2*phi0[j]+phi0[(j-1)%nx]))
            
            # Conservative form upwind.
            phi[j] += dt/dx *(-0.5*(phi0[j]**2-phi0[(j-1)%nx]**2) 
                              + 0.5*nu*(phi0[(j+1)%nx]-2*phi0[j]+phi0[(j-1)%nx]))
                    
            # Lax-wendroff.
        
            # # Trapezoidal scheme with intermediate predictor step.
            # phiS[j] += dt/dx *(-0.5*phi0[j]*(phi0[j]-phi0[(j-1)%nx]) 
            #                   + 0.5*nu*(phi0[(j+1)%nx]-2*phi0[j]+phi0[(j-1)%nx]))
        
        # b = dt*phi0 - dt/(2*dx) * (phiS**2 - np.roll(phiS, 1)**2)
        # phi = np.linalg.matmul(Ainv, b)
        
        # phi = phiS.copy()
        
        # Update old timestep.
        phi0 = phi.copy()
        # phi0 = phiS.copy()
                
        # Plot values for debugging.
        if i % 5 == 0:
            # plt.plot(X, phi)
            plt.plot(X, phi0)
            plt.plot(X, U[i])
            # plt.plot(X, phiS)
            # plt.plot(X, phiA)
            plt.xlim([x0, xL])
            plt.ylim([-0.1, 1.1])
            plt.show(), plt.close()
        