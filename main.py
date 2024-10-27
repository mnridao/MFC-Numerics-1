# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 19:48:13 2024

@author: guest997558
"""

import numpy as np
import matplotlib.pyplot as plt

from grids import Grid1D
from model import Model
from solver import Solver

from referenceSolution import ReferenceSolution

if __name__ == "__main__":
        
    # Generate grid.
    xbounds = [0, 10]
    dx = 0.01
    nx = int((xbounds[1] - xbounds[0])/dx)
    grid = Grid1D(xbounds, nx)
        
    # Solver parameters.
    endtime = 30
    dt = 0.025
    nt = int(endtime/dt)
        
    # Initial condition.
    # phi0 = np.exp(-2 * (grid.X - 0.5 * grid.xbounds[1])**2)
    phi0 = np.exp(-(grid.X-3)**2/2) 
    grid.phi = phi0.copy()
    plt.plot(grid.X, grid.phi)
    
    model = Model(grid, "rk4", dt)
    solver = Solver(model, dt, nt)
    solver.plotResults=True 
    solver.plotEveryNTimesteps=5
    # solver.run()

        #%%

    sol = ReferenceSolution(model.params, grid, endtime, nt)
    sol.evaluate(phi0)
    
    #%%
    for i in range(nt):
        if i%20 == 0:
            plt.plot(sol.phi[:, i])
            plt.show(), plt.close()