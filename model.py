# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 22:18:06 2024

@author: guest997558
"""

import advectionSchemes as schemes

debug=False

class Parameters:   
    """ 
    Stores the physical parameters of the viscous Burgers' equation.
    """
    
    def __init__(self):
        
        self.nu = 0.006
        self.mu = 1
        
class Model:
    """ 
    Class responsible for specifying:
        - Spatial discretisation 
        - Grid information
        - Physical parameters of the PDE
        
    Built in forward-in-time temporal discretisation and centred-difference 
    approximation for diffusion term.
    """
    
    def __init__(self, grid, spatialScheme, params=Parameters(), linear=False, 
                 dt=None):
        """ 
        Inputs
        ------
        grid : Grid1D object 
               Class with grid information.
        spatialScheme : string
                        which (advection) scheme to use.
        params : Parameters object 
                 Contains physical parameters of the PDE 
        linear : bool 
                 Flag for whether PDE is linear/nonlinear
        """
        
        # Useful things to know.
        self.params = params 
        self.grid = grid
        self.linear = linear
        
        # Scheme information.
        self.setSpatialScheme(spatialScheme, dt)
        
    def setSpatialScheme(self, scheme, dt=None):
        """ 
        Specify which advection scheme to use.
        """
        
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
        Step the model forward in time.
        """
        if self.spatialScheme.explicit:
            
            # Step forward in time (only doing forward euler here).
            self.grid.phi += dt*self.spatialScheme.discretise(self.grid.phi, dt)
            
        else:
            
            # Solve a matrix equation.
            self.grid.phi = self.spatialScheme.solve(self.grid.phi, dt)
        
        if debug:
            if self.linear:
                ci = self.params.mu*dt/self.grid.dx 
                di = self.params.nu*dt/self.grid.dx**2
            else:
                ci = self.grid.phi.max()*dt/self.grid.dx 
                di = self.params.nu*dt/self.grid.dx**2
            print(f"Pe{ci/di}, \t{ci=}, \t{di=}, \tu={self.grid.phi.max()}",)
            # print(f"c+2d={ci+2*di}")
            
        return self.grid.phi
    