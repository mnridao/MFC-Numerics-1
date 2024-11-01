# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 22:18:51 2024

@author: guest997558
"""

import numpy as np
import matplotlib.pyplot as plt

import helpers
import plotters

if __name__ == "__main__":
    
    ## Check stability for linear case ##
    
    # Stability stuff.
    c = 1.5
    d = 1.5
    Pe = c/d
    
    u = 0.1
    nu = 0.02
    
    dx = Pe*nu/u
    dt = c*dx/u
    
    # Generate grid.
    grid = helpers.setupGrid(x0=0, xL=10, dx=dx)
    dx = grid.dx
        
    # Problem Parameters.
    params = helpers.setupParameters(nu=nu, mu=u)
            
    # Initial condition.
    ic = helpers.smoothWave
    icArgs = (grid.xbounds[1], )
    
    # Setup solver.
    solver = helpers.setupSolver(dt, endtime=1000, scheme="upwind", grid=grid, 
                                 ic=ic, icArgs=icArgs, 
                                 params=params, linear=True,
                                 plotResults=True, plotEveryN=5)
    
    c = params.mu*solver.dt/grid.dx 
    d = params.nu*solver.dt/grid.dx**2
    print(f"{Pe=}, \t{c=}, \t{d=}\t, dx={grid.dx}, \tdt={solver.dt}")
    solver.run()