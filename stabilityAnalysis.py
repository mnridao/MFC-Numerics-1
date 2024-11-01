# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 22:18:51 2024

@author: guest997558
"""

import matplotlib.pyplot as plt

import helpers
import plotters

if __name__ == "__main__":
    
    ## Check stability for linear case ##
    
    ######### USER INPUT HERE #########
    c = 0.5
    d = 0.5
    u = 1
    nu = 0.02
    # dx = 0.01
    
    scheme = "upwind"
    # scheme = "centredDifference"
    endtime=100
    plotEveryIteration=False
    plotEveryNTimesteps=1
    
    ###################################
    
    # Stability stuff.
    Pe = c/d
    dx = Pe*nu/u
    # nu = dx*u/Pe
    dt = c*dx/u
    
    # Generate grid.
    grid = helpers.setupGrid(x0=0, xL=10, dx=dx)
    dx = grid.dx
                    
    # Initial condition.
    ic = helpers.piecewiseWave
    icArgs = ()
    
    # Setup solver.
    params = helpers.setupParameters(nu=nu, mu=u)
    solver = helpers.setupSolver(dt, endtime=endtime, scheme=scheme, grid=grid, 
                                 ic=ic, icArgs=icArgs, 
                                 params=params, linear=True,
                                 plotResults=plotEveryIteration, plotEveryN=plotEveryNTimesteps)
            
    # # Add mass custom equation.
    # solver.addCustomEquation("mass", helpers.massCustom)
    
    ci = params.mu*solver.dt/grid.dx 
    di = params.nu*solver.dt/grid.dx**2
    print(f"{Pe=}, \t{ci=}, \t{di=}\t, dx={grid.dx}, \tdt={solver.dt}, \tnu={params.nu}")
    solver.run()
    
    # Plot final thing.
    plotters.defaultPlotter(solver.model.grid)
    
    # # Plot the mass.
    # mass = solver.getCustomData("mass")[2:]
    # plt.plot(mass)
    # plt.show(), plt.close()