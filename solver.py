# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 17:22:19 2024

@author: guest997558
"""

import numpy as np

import plotters as plt

class Solver:
    """ 
    """
    
    def __init__(self, model, dt, nt):
        """ """
        
        # Solver parameters.
        self.model   = model
        self.dt     = dt
        self.nt     = nt
        
        # Storage parameters.
        self.store = False
        self.history = None
        
        # Plotter parameters.
        self.plotResults = False
        self.plotEveryNTimesteps = 1
        self.plotter = plt.defaultPlotter
        
    def run(self):
        
        # Initialise storage array.
        if self.store:
            self.history = np.zeros(shape=(self.nt+1, *self.model.grid.phi.shape), 
                                    dtype=self.model.grid.phi.dtype)
            self.history[0, ...] = self.model.grid.phi
         
        for i in range(1, self.nt+1):
                        
            # Calculate new time step value.
            phi = self.model.step(self.dt)
                          
            if self.plotResults:
                
                if i % self.plotEveryNTimesteps == 0:
                    self.plotter(self.model.grid)
                                                                        
            # Update the old timestep value.
            self.model.grid.phi = phi.copy()
            
            # Store if necessary.
            if self.store:
                self.history[i, ...] = self.model.grid.phi
                        
    def setNewTimestep(self, dt, endtime):
        """ 
        """        
        # Set new dt and calculate new nt.
        self.dt = dt 
        self.nt = int(np.ceil(endtime/dt))