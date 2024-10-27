# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 16:55:31 2024

@author: guest997558
"""

import numpy as np

class Grid1D:
    """ 
    """
    
    def __init__(self, xbounds, nx):
        """ """
        
        self.nx = nx
        self.xbounds = xbounds
        self.periodic = True   # Remove user option for now.
                
        self.setupGrid()
    
    def setupGrid(self):
        """ 
        """
        
        # Check for periodic boundary conditions.
        numPoints = self.nx + (not self.periodic)*1
        
        # Calculate grid spacing from inputs.
        self.dx = (self.xbounds[1] - self.xbounds[0]) / self.nx
        
        # Setup array representing the grid.
        self.X = np.linspace(self.xbounds[0], self.xbounds[1], numPoints, 
                             endpoint=(not self.periodic))
        
        # Initialise the default state variable.
        self.resetFields()
    
    def resetFields(self):
        """ 
        """
        self.phi = np.zeros_like(self.X)
    
    def setNewGridBounds(self, xbounds):
        """ 
        """
        self.xbounds = xbounds 
        self.setupGrid()
    
    def setNewGridSpacing(self, dx):
        """ 
        """
        self.dx = dx 
        self.nx = int((self.xbounds[1] - self.xbounds[0]) / self.dx)