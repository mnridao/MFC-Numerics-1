# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 16:15:20 2024

@author: guest997558
"""

import numpy as np

class AdvectionDiffusion1D:
    """ 
    Class for handling basic Explicit Eulerian discretisations of the
    1d advection-diffusion equation.
    
    d phi/dt = -mu*phi*dphi/dx + nu*d2u/dx2
    
    """
    
    def __init__(self, grid, params, linear=False):
        """ 
        """
        
        # Consider a model class instead? that way can change params.
        self.grid = grid
        self.params = params
        
        # Default class variables (check consistent bc).
        self.linear = linear
        self.F = conservativeNonlinearFlux if not linear else linearFlux
        self.advectionScheme = upwindPeriodic
        self.diffusionScheme = centredDifferenceD2Periodic
        
        # This base class is for explicit schemes.
        self.explicit=True

    def discretise(self, phi, dt=None):
        """ 
        """
                
        return self.rhs(phi)
    
    def rhs(self, phi):
        """ 
        """
        
        spatialTerms = np.zeros_like(phi)
        for j in range(self.grid.nx):
                        
            spatialTerms[j] = (-self.params.mu*self.advectionScheme(j, self.F, self.grid) 
                               +self.params.nu*self.diffusionScheme(j, self.grid))
        
        return spatialTerms

class AdvectionDiffusionRK4(AdvectionDiffusion1D):
    """ 
    RK4 implementation for the advection-diffusion equation.
    """
    
    def __init__(self, grid, params, linear=False):
        """ 
        """
        super().__init__(grid, params, linear)
        
    def discretise(self, phi, dt):
        """ 
        """
        
        k1 = self.rhs(phi)
        k2 = self.rhs(phi+0.5*dt*k1)
        k3 = self.rhs(phi+0.5*dt*k2)
        k4 = self.rhs(phi+dt*k3)
        
        return (k1 + 2*k2 + 2*k3 + k4)/6
    
class AdvectionDiffusionImplicit:
    """
    I don't think that this is working.
    """
    def __init__(self, grid, dt, params, linear=False):
        
        # Consider a model class instead? that way can change params.
        self.grid = grid
        self.params = params
        self.linear = linear
        
        self.setupAinv(dt)
        
        # This scheme is implicit.
        self.explicit=False
        
    def setupAinv(self, dt):
        """ 
        """
        
        # If the problem is non-linear, we deal with advection explicitly.
        
        A = np.zeros(shape=(self.grid.nx, self.grid.nx))
        
        diagMain = 1 + (self.params.nu + self.linear*self.params.mu)*dt/self.grid.dx
        diagUpper = (-0.5*self.params.nu - self.linear*self.params.mu)*dt/self.grid.dx 
        diagLower = -self.params.nu*dt/self.grid.dx
        
        np.fill_diagonal(A, diagMain)
        
        for i in range(self.grid.nx - 1):
            A[i, i+1] = diagUpper
            A[i+1, i] = diagLower
        
        A[0, -1] = diagLower
        A[-1, 0] = diagUpper
        
        self.A = A # debugging.
        self.Ainv = np.linalg.inv(A)
                
    def solve(self, phi, dt):
        """
        """
        
        # Setup b vector (Ax=b).
        b = phi*dt - self.linear*0.5*(phi**2-np.roll(phi, 1)**2)*dt/self.grid.dx
        
        # Solve the equation.
        return np.linalg.matmul(self.Ainv, b)
    
### FLUX FUNCTIONS ###

def conservativeNonlinearFlux(phi):
    """ 
    """
    return 0.5*phi**2

def linearFlux(phi):
    """ 
    Linear as in linear advection.
    
    User would need to set mu as the constant advection velocity.
    """
    return phi

### ADVECTION DISCRETISATION ### 

def upwindPeriodic(j, F, grid):
    """ 
    Upwind discretisation for advection with periodic boundary conditions.
    
    Inputs
    ------
    F   : callable function
          Flux of the advection term.
    phi : np array
          prognostic variable.
    dx  : float 
          grid spacing.
    """
    
    return (F(grid.phi[j])-F(grid.phi[(j-1)%grid.nx]))/grid.dx

def centredDifferenceD1Periodic(j, F, grid):
    """ 
    Cnetred diff discretisation for advection with periodic boundary conditions.
    
    Inputs
    ------
    F   : callable function
          Flux of the advection term.
    phi : np array
          prognostic variable.
    dx  : float 
          grid spacing.
    """
    
    return 0.5*(F(grid.phi[(j+1)%grid.nx])-F(grid.phi[(j-1)%grid.nx]))/grid.dx

def upwindPeriodicOld(j, grid):
    """ 
    Upwind discretisation for advection with periodic boundary conditions.
    
    Inputs
    ------
    F   : callable function
          Flux of the advection term.
    phi : np array
          prognostic variable.
    dx  : float 
          grid spacing.
    """    
    return 0.5*(grid.phi[j]**2-grid.phi[(j-1)%grid.nx]**2)/grid.dx

def centredDifferenceD1PeriodicOld(j, grid):
    """ 
    Leapfrog discretisation for advection with periodic boundary conditions.
    
    Inputs
    ------
    F   : callable function
          Flux of the advection term.
    phi : np array
          prognostic variable.
    dx  : float 
          grid spacing.
    """
    
    return 0.25*(grid.phi[(j+1)%grid.nx]**2-grid.phi[(j-1)%grid.nx]**2)/grid.dx

### DIFFUSION DISCRETISATION ### 

def centredDifferenceD2Periodic(j, grid):
    """ 
    """
    
    return (grid.phi[(j+1)%grid.nx]-2*grid.phi[j]+grid.phi[(j-1)%grid.nx])/grid.dx**2