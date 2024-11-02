# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 01:00:06 2024

@author: guest997558
"""

import helpers
import plotters

def generateSolutionPlots():
    """ 
    """
    
    # Setup the grid.
    xbounds=[0, 10]
    dx = 0.01
    grid = helpers.setupGrid(*xbounds, dx)
    
    # Initial condition.
    ic=helpers.bump1 
    icArgs=()
    
    # Setup the solver.
    endtime=12
    solver = helpers.setupSolver(dt=0.001, endtime=endtime, scheme="upwind", grid=grid, 
                                 ic=ic, icArgs=icArgs, linear=False)
    solver.store=True
    solver.run()
    
    # Plot the solution.
    plotters.plotSolutionFancy(solver.model.grid, endtime, ic, icArgs, "solution.pdf", saveFig=True)
    plotters.plotSolutionVeryFancy(solver, "solutionEvolution.pdf", saveFig=True)
 
if __name__ == "__main__":
    
    generateSolutionPlots()