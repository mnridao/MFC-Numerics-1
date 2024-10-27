# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 22:15:06 2024

@author: guest997558
"""

import matplotlib.pyplot as plt

def defaultPlotter(grid):
    """ 
    """
    plt.plot(grid.X, grid.phi)
    plt.xlim(grid.xbounds)
    plt.ylim([-0.1, 1.1])
    plt.xlabel("X")
    plt.ylabel("phi")
    plt.show(), plt.close()