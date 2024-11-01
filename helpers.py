# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 19:11:12 2024

@author: guest997558
"""

import numpy as np

from grids import Grid1D
from solver import Solver
from model import Model, Parameters

### MISC ###

def setupParameters(nu=None, mu=None):
    params = Parameters()
    if nu:
        params.nu = nu 
    if mu:
        params.mu = mu
    return params

def setupGrid(x0, xL, dx):
    """ 
    """
    # Calculate no. of grid points.
    nx = int((xL - x0)/dx)
    
    # Recalculate dx for consistency.
    dx = (xL-x0)/nx
    
    return Grid1D([x0, xL], nx)
    
def setupSolver(dt, endtime, scheme, grid, ic, icArgs, params=Parameters(), 
                linear=False, plotResults=False, plotEveryN=1):
    """ 
    """
    # Calculate the no. of time steps.
    nt = int(np.floor(endtime/dt))
    
    # Recalculate dt for consistency.
    dt = endtime/nt
    
    # Setup model and then solver.
    model = Model(grid, scheme, params, linear, dt)
    model.grid.phi = ic(grid.X, *icArgs)    
    solver = Solver(model, dt, nt)
    solver.plotResults=plotResults
    solver.plotEveryNTimesteps=plotEveryN
    
    return solver

### STABILITY STUFF ###

def mass():
    """ 
    Should not increase ideally.
    """

### CONVERGENCE STUFF ###

def l2(phi, phiA, dx):
    """ 
    """
    return (np.sqrt(np.sum(dx*np.power(phi - phiA, 2)))/
            np.sqrt(np.sum(dx*np.power(phiA, 2))))

def order(err, dxs):
    """ 
    """
    n, _ = np.polyfit(np.log(dxs), np.log(err), 1)
    return n

### INITIAL CONDITIONS ###

def bump1(X):
    """ 
    """
    return np.exp(-(X-3)**2/2) 

def bump2(X, xL):
    """ 
    """
    return np.exp(-2 * (X - 0.5 * xL)**2)

def smoothWave(X, xL):
    """ 
    """
    return 0.5*(1-np.cos(2*np.pi*(X/xL)))