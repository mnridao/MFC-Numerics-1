# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 23:33:31 2024

@author: guest997558
"""

import numpy as np
import matplotlib.pyplot as plt

from referenceSolution import ReferenceSolution
import helpers
import plotters

def runDifferentCourantNumbers(cs, nxs, endtime, scheme, 
                               ic, icArgs=(), debug=False):
    """ 
    """
    
    # Iterate over different initial Courant number options.
    errs = np.zeros(shape=(len(cs), len(nxs)))
    for i, c in enumerate(cs):
        
        # Find initial umax (1 from initial condition).
        u0max=1
        
        errs[i, :] = runDifferentResolutionsConstC(nxs, c, u0max, endtime, scheme, 
                                                   helpers.setupParameters(), 
                                                   ic, icArgs, debug=debug)
    
    return errs

def runDifferentResolutionsConstC(nxs, c, u0max, endtime, scheme, params, 
                                  ic, icArgs=(), debug=False):
    """ 
    """
    err = np.zeros_like(nxs, dtype=np.float64)
    for i, nx in enumerate(nxs):
        
        # Update the timestep.
        dx = 1./nx 
        dt = c*dx/u0max
        nt = int(endtime/dt)
        
        # Calculate error.
        err[i] = runResolution(dx, nx, dt, nt, endtime, scheme, ic, icArgs, debug=debug)
            
    return err

def runResolution(dx, nx, dt, nt, endtime, scheme, ic, icArgs, debug=False):
    """ 
    """
    
    c = 1*dt/dx 
    d = 0.006*dt/dx**2
    print(f"c={c}, \td={d}, \tPe={c/d} \tdt={dt}, \t dx={dx}, \t nt={nt}, \t nx={nx}, \tnt={nt}, \t, nt*dt={nt*dt}")
        
    # Update the grid and solver.
    grid = helpers.setupGrid(x0=0, xL=10, dx=dx)
    solver = helpers.setupSolver(dt, endtime=endtime, scheme=scheme, grid=grid, 
                                 ic=ic, icArgs=icArgs, linear=False)
    
    # Run the analytic solution.
    sol = ReferenceSolution(solver.model.params, grid, endtime, nt=500)
    sol.evaluate(solver.model.grid.phi)
    
    # Run the solver.
    solver.run()
    
    # Plot the last iteration.
    if debug:
        model = solver.model 
        plt.plot(model.grid.X, model.grid.phi)
        plt.plot(model.grid.X, sol.phi[:, -1], 'k--')
        plt.ylim([-0.1, 1.1])
        plt.xlim(model.grid.xbounds)
        plt.grid()
        plt.show(), plt.close()
    
    return helpers.l2(grid.phi, sol.phi[:, -1], grid.dx)

def calculateOrders(dxs, errs):
    
    n = np.zeros_like(dxs, dtype=np.float64)
    for i, err in enumerate(errs):
        ni, _ = np.polyfit(np.log(dxs), np.log(err), 1)
        n[i] = ni 
    return n

def displayOrders(n, cs, scheme):
    
    print(scheme)
    for (ni, c) in zip(n, cs):
        print(f"{c=}, \tn={ni}")

def runConvergenceExperiment():
    """ 
    """
    
    # Initial condition for this experiment.
    ic = helpers.smoothWave
    icArgs = (10,)
    
    # Courant no. ranges.
    cs = [1, 0.5, 0.25]
    nxs = [10, 20, 40, 80, 160]
    endtime = 0.1 
    
    # Errors for upwind scheme compared to Fourier model.
    errsU = runDifferentCourantNumbers(cs, nxs, endtime, "upwind", 
                                       ic, icArgs, debug=False)
    
    # Plot the upwind errors.
    x = 1e-2, 3e-2
    y0 = 3e-4
    y1 = y0*(x[1]/x[0])
    orderLine = [x, [y0, y1], "first-order"]
    
    # Plot the error.
    dxs = 1./np.array(nxs)
    plotters.plotErrors(errsU, dxs, cs, orderLine, "upwindConvergence.png", saveFig=True)
    
    # Calculate the order.
    n = calculateOrders(dxs, errsU)
    displayOrders(n, cs, "upwind")
    
    # Centred difference not working.
    
if __name__ == "__main__":
    
    runConvergenceExperiment()