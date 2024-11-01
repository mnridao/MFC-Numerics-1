# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 23:33:31 2024

@author: guest997558
"""

import numpy as np
import matplotlib.pyplot as plt

from grids import Grid1D
from model import Model, Parameters
from solver import Solver
from referenceSolution import ReferenceSolution
import helpers
import plotters

def runDifferentResolutionsConstD():
    """ 
    """
    err = np.zeros_like(nxs)
    for i, nx in enumerate(nxs):
        pass
    
    return err

def runDifferentResolutionsConstC(nxs, c, u0max, endtime, scheme, params, ic, icArgs=(), 
                                  debug=False):
    """ 
    """
    err = np.zeros_like(nxs, dtype=np.float64)
    for i, nx in enumerate(nxs):
        
        # Update the timestep.
        dx = 1./nx 
        dt = c*dx/u0max
        nt = int(endtime/dt)
        
        # Calculate error.
        err[i] = runResolution(dx, nx, dt, nt, endtime, ic, icArgs, debug=debug)
            
    return err

def runResolution(dx, nx, dt, nt, endtime, ic, icArgs, debug=False):
    """ 
    """
    
    print(f"c={u0max*dt/dx}, \tdt={dt}, \t dx={dx}, \t nt={nt}, \t nx={nx}, \tnt={nt}, \t, nt*dt={nt*dt}")
    
    # Update the grid.
    grid = Grid1D(xbounds, nx)
    model = Model(grid, scheme, params, dt=dt)
    model.grid.phi = ic(grid.X, *icArgs)
    
    # Run the analytic solution.
    sol = ReferenceSolution(params, grid, endtime, nt)
    sol.evaluate(model.grid.phi)
    
    # Update the solver.
    solver = Solver(model, dt, nt)
    solver.run()
    
    # Plot the last iteration.
    if debug:
        plt.plot(grid.X, model.grid.phi)
        plt.plot(grid.X, sol.phi[:, -1], 'k--')
        plt.ylim([-0.1, 1.1])
        plt.xlim(grid.xbounds)
        plt.grid()
        plt.show(), plt.close()
    
    return helpers.l2(grid.phi, sol.phi[:, -1], grid.dx)

if __name__ == "__main__":
    
    c = 0.5  # Initial Courant no.
    nxs = [10, 20, 40, 80, 160]
    
    endtime = 0.1
    err = np.zeros(shape=len(nxs))
    dts = []
    
    # Generate grid.
    xbounds = [0, 10]
    dx = 0.01
    nx = int((xbounds[1] - xbounds[0])/dx)
    grid = Grid1D(xbounds, nx)
    
    # Initial condition.
    ic = helpers.smoothWave
    icArgs=(grid.xbounds[1], )
    
    # Find initial umax.
    phi0 = ic(grid.X, *icArgs)
    u0max = np.max(phi0)
    plt.plot(grid.X, grid.phi)

    # Problem Parameters.
    params = Parameters()
    params.nu = 0.02
    # params.nu = 1e-6
    
    schemes = ["upwind", "centredDifference", "rk4"]
    schemes = ["upwind"]
    errs = np.zeros(shape=(len(schemes), len(nxs)))
    
    for i, scheme in enumerate(schemes):
        print(f"### scheme: {scheme},\t {c=},\t Pe={helpers.peclet(u0max, params.nu, grid.dx)} ###\n")

        errs[i, :] = runDifferentResolutionsConstC(nxs, c, u0max, endtime, 
                                                   scheme, params, ic, icArgs,
                                                   True)
  
    #%%
    dxs = 1./np.array(nxs)
    # plotters.plotErrors(errs, 1./np.array(nxs), schemes, grid, "", False)
    
    # Plot the error.
    err = errs[0]
    
    plt.figure()
    plt.plot(dxs, err, '-o')
    plt.yscale("log")   
    plt.xscale("log")
    plt.grid(which="both")
    plt.show(), plt.close()
    
    # Calculate experimental order of convergence.
    n, _ = np.polyfit(np.log(dxs), np.log(err), 1)
    print(n)