# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 22:15:06 2024

@author: guest997558
"""
import numpy as np
import matplotlib.pyplot as plt

def defaultPlotter(grid):
    """
    For debugging.
    """
    plt.plot(grid.X, grid.phi)
    plt.xlim(grid.xbounds)
    plt.ylim([-0.1, 1.1])
    plt.xlabel("X")
    plt.ylabel("phi")
    plt.show(), plt.close()
    
def plotSolutionAtTime(grid, ic, icArgs, ybounds, t, figname, saveFig):
    """ 
    """
    
    # Initial condition.
    phi0 = ic(grid.X, *icArgs)
    
    plt.plot(grid.X, grid.phi, label="numerical")
    plt.plot(grid.X, phi0, 'k--', label="IC")
    plt.xlim(grid.xbounds)
    plt.ylim(ybounds)
    plt.xlabel("X")
    plt.ylabel(r"$\phi$")
    plt.legend()
    plt.grid(which="both")
    
    xpos = 0.8*(grid.xbounds[1]-grid.xbounds[0])
    ypos = 0.8*ybounds[1]
    plt.text(xpos, ypos, f"t={t:.2f} s")
    
    if saveFig:
        plt.tight_layout()
        figname = figname if figname else "plot.pdf"
        plt.savefig(figname)
    
def plotSolutionsAtTime(grids, labels, ic, icArgs, ybounds, t, figname, saveFig):
    """ 
    """
    
    plt.figure()
    for i, grid in enumerate(grids):
        
        # Initial condition.
        if i == 0:
            phi0 = ic(grid.X, *icArgs)
            plt.plot(grid.X, phi0, "k--", label="IC")
        
        # Numerical solution.
        plt.plot(grid.X, grid.phi, label=labels[i])
        
    plt.xlim(grid.xbounds)
    plt.ylim(ybounds)
    plt.xlabel("X")
    plt.ylabel(r"$\phi$")
    plt.legend()
    plt.grid(which="both")
    
    xpos = 0.8*(grid.xbounds[1]-grid.xbounds[0])
    ypos = 0.8*ybounds[1]
    plt.text(xpos, ypos, f"t={t:.2f} s")
    
    if saveFig:
        plt.tight_layout()
        figname = figname if figname else "plot.pdf"
        plt.savefig(figname)
    
def plotMasses(solvers, labels, figname, saveFig):
    """ 
    """
    
    plt.figure()
    for i, sol in enumerate(solvers):
        
        # Extract mass information from solver.
        mass = sol.getCustomData("mass")
        t = np.arange(0, sol.dt*(sol.nt+1), sol.dt)
        
        # Remove last element of t if necessary.
        if mass.shape != t.shape:
            t = t[:-1]
        
        # Plot mass evolution with time for the current solver.
        plt.plot(t, np.abs(mass), label=labels[i])
        
    # Analytical mass for linear burger equation for problem I have.
    plt.plot(t, np.ones_like(mass)*5, "k--", label="analytic")    
    
    plt.legend()  
    plt.grid(which="both")
    plt.yscale("log")
    plt.xlabel("Time [s]")
    plt.ylabel("Absolute Mass")
    plt.xlim([t[0], t[-1]])
    
    if saveFig:
        plt.tight_layout()
        plt.savefig(figname)
    
def plotErrors0ld(errs, dxs, schemes, grid, figname, saveFig):
    """ 
    """
    
    plt.figure()
    for err, s in zip(errs, schemes):
        plt.plot(dxs, err, '-o', label=s)
    plt.yscale("log")
    plt.xscale("log")
    plt.grid(which="both")
    
    plt.legend()
    plt.xlabel("$\Delta x")
    plt.ylabel(r"$l_2$ error")
    
    if saveFig:
        plt.tight_layout()
        plt.savefig(figname)
    
def plotErrors(errs, dxs, labels, orderLine, figname="", saveFig=False):
    """ 
    """
    
    plt.figure()
    for i, err in enumerate(errs):
        plt.plot(dxs, err, '-o', label=f"c={labels[i]}")
    
    # Plot order line.
    x, y, order = orderLine
    plt.plot(x, y, 'k--', label=order)
        
    plt.yscale("log")   
    plt.xscale("log")
    plt.grid(which="both")
    plt.legend()
    plt.xlabel(r"$\Delta x$")
    plt.ylabel(r"$l_2$ error")
    plt.xlim([dxs[-1], dxs[0]])
    
    if saveFig:
        plt.tight_layout()
        plt.savefig(figname)
    
def plotAmplifcationFactors(cs, ds, maxAs, cases, scheme, 
                            figname="", saveFig=False):
    """ 
    """
    
    # Unstable is unstable.
    maxAs[maxAs>1] = np.nan
    
    # Create grid.
    XS, YS = np.meshgrid(ds, cs)
        
    plt.figure()
    plt.contourf(XS, YS, maxAs)
    
    # Plot case locations.
    for (c, d) in cases:
        plt.plot(d, c, '*k', markersize=8)
    
    plt.xlabel("d")
    plt.ylabel("c")
    plt.grid()
    plt.xlim([ds[0], ds[-1]])
    plt.ylim([cs[0], cs[-1]])
    plt.title(scheme)
    
    if saveFig:
        plt.tight_layout()
        plt.savefig(figname)
        
def plotSolutionFancy(grid, t, ic, icArgs, figname="", saveFig=False):
    
    # Initial condition
    phi0 = ic(grid.X, *icArgs)
    
    plt.figure()
    plt.plot(grid.X, phi0, "k--", label="IC")
    plt.plot(grid.X, grid.phi, label="numerical")
    
    plt.text(8, 0.9, f"t={t:.2f} s")
    plt.xlim(grid.xbounds)
    plt.grid(which="both")
    plt.ylim([-0.1, 1.1])
    plt.legend()
    plt.xlabel("X")
    plt.ylabel("phi")
    
    if saveFig:
        plt.tight_layout()
        plt.savefig(figname)

def plotSolutionVeryFancy(solver, figname="", saveFig=False):
    
    # Plot the solution with time.
    t = np.arange(0, solver.dt*(solver.nt+1), solver.dt)
    XS, YS = np.meshgrid(solver.model.grid.X, t)
    
    plt.figure()
    plt.contourf(XS, YS, solver.history)
    cbar = plt.colorbar()
    cbar.set_label(r"$\phi$")
    
    plt.xlabel("X [m]")
    plt.ylabel("Time [s]")
    
    if saveFig:
        plt.tight_layout()
        plt.savefig(figname)