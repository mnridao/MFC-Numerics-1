# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 22:15:06 2024

@author: guest997558
"""

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
    
def plotErrors(errs, dxs, schemes, grid, figname, saveFig):
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