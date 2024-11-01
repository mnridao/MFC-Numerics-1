# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 22:18:06 2024

@author: guest997558
"""

import advectionSchemes as schemes

debug=False

class Parameters:   
    """ 
    """
    
    def __init__(self):
        
        self.nu = 0.006
        self.mu = 1
        
class Model:
    """ 
    """
    
    def __init__(self, grid, spatialScheme, params=Parameters(), linear=False, 
                 dt=None):
        """ 
        """
        
        # Useful things to know.
        self.params = params 
        self.grid = grid
        self.linear = linear
        
        # Scheme information.
        self.setSpatialScheme(spatialScheme, dt)
        
    def setSpatialScheme(self, scheme, dt=None):
        
        if scheme == "upwind":
            self.spatialScheme = schemes.AdvectionDiffusion1D(
                self.grid, self.params, linear=self.linear)
        
        elif scheme == "centredDifference":
            self.spatialScheme = schemes.AdvectionDiffusion1D(
                self.grid, self.params, linear=self.linear)
            self.spatialScheme.advectionScheme = schemes.centredDifferenceD1Periodic
        
        elif scheme == "rk4":
            self.spatialScheme = schemes.AdvectionDiffusionRK4(
                self.grid, self.params, linear=self.linear)
        
        elif scheme == "implicit":
            self.spatialScheme = schemes.AdvectionDiffusionImplicit(
                self.grid, dt, self.params, linear=self.linear)
        
    def step(self, dt):
        """ 
        """
        if debug:
            ci = self.params.mu*dt/self.grid.dx 
            di = self.params.nu*dt/self.grid.dx**2
            Pe=ci/di
            print(f"{Pe=}, \t{ci=}, \t{di=}\t, dx={self.grid.dx}, \tdt={dt}")
            
        if self.spatialScheme.explicit:
            
            # Step forward in time (only doing forward euler here).
            self.grid.phi += dt*self.spatialScheme.discretise(self.grid.phi, dt)
            
        else:
            
            # Solve a matrix equation.
            self.grid.phi = self.spatialScheme.solve(self.grid.phi, dt)
        
        return self.grid.phi
    