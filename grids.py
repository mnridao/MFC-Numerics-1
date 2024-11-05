# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 16:55:31 2024

@author: guest997558
"""

import numpy as np

class Grid1D:
    """ 
    Class responsible for generating the grid.
    """
    
    def __init__(self, xbounds, nx):
        """ 
        Inputs
        ------
        xbounds: list
                 Contains lower and upper bounds of the domain
        nx     : int
                 Number if grid points in the domain.
        """
        
        self.nx = nx
        self.xbounds = xbounds
        self.periodic = True   # Remove user option for now.
                
        self.setupGrid()
    
    def setupGrid(self):
        """ 
        Generates the grid (assumes periodic boundary conditions)
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
        Reset the prognostic variable values.
        """
        self.phi = np.zeros_like(self.X)
    
    def setNewGridBounds(self, xbounds):
        """ 
        Regenerates the grid based on new domain bounds.
        """
        self.xbounds = xbounds 
        self.setupGrid()
    
    def setNewGridSpacing(self, dx):
        """ 
        Reset the grid spacing (recalculates nx)
        """
        self.dx = dx 
        self.nx = int((self.xbounds[1] - self.xbounds[0]) / self.dx)