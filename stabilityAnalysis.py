# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 22:18:51 2024

@author: guest997558
"""

import helpers
import plotters

def stabilityDemo(c, d, scheme, endtime, ybounds=[-0.1, 1.1], plot=True):
    """ 
    """
   
    # Keep these fixed.
    u = 1
    nu = 0.02
   
    # Calculate stability stuff.
    Pe = c/d
    dx = Pe*nu/u
    dt = c*dx/u
    
    # Generate grid.
    grid = helpers.setupGrid(x0=0, xL=10, dx=dx)
    dx = grid.dx
        
    ic = helpers.smoothWave
    icArgs = (grid.xbounds[1],)
    
    # Setup solver.
    params = helpers.setupParameters(nu=nu, mu=u)
    solver = helpers.setupSolver(dt, endtime=endtime, scheme=scheme, grid=grid, 
                                 ic=ic, icArgs=icArgs, 
                                 params=params, linear=True)
    
    # Add mass to solver.
    solver.addCustomEquation("mass", helpers.massCustom)
    
    ci = params.mu*solver.dt/grid.dx 
    di = params.nu*solver.dt/grid.dx**2
    print(f"{Pe=}, \t{ci=}, \t{di=}\t, dx={grid.dx}, \tdt={solver.dt}, \tnu={params.nu}")
    solver.run()
    
    # Plot final thing.
    if plot:
        plotters.plotSolutionAtTime(grid, ic, icArgs, ybounds, endtime, figname="", saveFig=False)
    
    return solver

def runStabilityDemo(cases, scheme, ic, icArgs, endtime, ybounds, figname, saveFig):
    """ 
    """
    
    # I promise you that using lists here is not less safe than arrays.
    grids, labels = [], [] 
    
    # Iterate through the different stability cases.
    for (c, d) in cases:
        
        # Run the solver for the current stability case.
        sol = stabilityDemo(c, d, scheme, endtime, plot=False)
        
        # For the plots.
        grids.append(sol.model.grid)
        labels.append(f"c={c}, d={d}")
    
    # Plot the stability results.
    plotters.plotSolutionsAtTime(grids, labels, ic, icArgs, 
                                 ybounds, endtime, figname, saveFig)
    
    return grids, labels
    
def runMassAnalysis(cases, scheme, ic, icArgs, endtime, figname, saveFig, plot=True):
    """ 
    """
    
    solvers, labels = [], []
    
    # Iterate through the stability cases.
    for (c, d) in cases:
        
        # Run the solver for the current stability case.
        sol = stabilityDemo(c, d, scheme, endtime, plot=False)
        
        # For the plots.
        solvers.append(sol)
        labels.append(f"c={c}, d={d}")
    
    # Plot the mass evolution.
    if plot:
        plotters.plotMasses(solvers, labels, figname, saveFig)
    
    return solvers, labels
    
def runStabilityExperiment():
    """ 
    """
    
    # Initial condition for this experiment.
    ic = helpers.smoothWave
    icArgs = (10,)
    
    # These are hard coded here.
    casesUpwind = [(0.8, 0.4), (0.5, 0.1)]
    casesCentred = [(0.7, 0.1), (0.5, 0.3), (0.4, 0.6)]
    
    # Stability demo plots for upwind advection.
    runStabilityDemo(casesUpwind, "upwind", ic, icArgs, endtime=1.5, 
                     ybounds=[-0.1, 1.1], 
                     figname="upwindStabilityDemo.pdf", saveFig=True)
    
    # Stability demo plots for centred advection - endtimes split for clarity.
    runStabilityDemo(casesCentred[:-1], "centredDifference", ic, icArgs, endtime=40, 
                     ybounds=[-0.5, 2], 
                     figname="centredStabilityDemo1.pdf", saveFig=True)
    
    # This case for centred difference blows up much quicker than the others.
    runStabilityDemo([casesCentred[-1]], "centredDifference", ic, icArgs, endtime=0.6, 
                     ybounds=[-0.5, 2], 
                     figname="centredStabilityDemo2.pdf", saveFig=True)
    
    # Mass evolution to check stability for upwind scheme.
    runMassAnalysis(casesUpwind, "upwind", ic, icArgs, endtime=10, 
                    figname="upwindMassEvolution.pdf", saveFig=True)
    
    # Mass evolution for first two cases of centred difference (different endtimes).
    sols1, labels1 = runMassAnalysis(casesCentred[:-1], "centredDifference", 
                                     ic, icArgs, endtime=150, 
                                     figname="", saveFig=False, plot=False)
    sols2, labels2 = runMassAnalysis([casesCentred[-1]], "centredDifference", 
                                     ic, icArgs, endtime=10, 
                                     figname="", saveFig=False, plot=False)
    
    sols = sols2 + sols1[::-1]
    labels = labels2 + labels1[::-1]
    plotters.plotMasses(sols, labels, "centredeMassEvolution.pdf", True)
    
if __name__ == "__main__":
        
    runStabilityExperiment()
           