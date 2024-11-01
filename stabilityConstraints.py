# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 13:34:52 2024

@author: guest997558
"""

import numpy as np
import matplotlib.pyplot as plt

def upwindAmp(kdx, c, d):
    A = 1 + (2*d+c)*(np.cos(kdx) - 1) - 1j*c*np.sin(kdx)
    return A*A.conjugate()

def centredAmp(kdx, c, d):
    A = 1 + 2*d*(np.cos(kdx)-1) - 1j*c*np.sin(kdx)
    return A*A.conjugate()

if __name__ == "__main__":
    
    # Create ranges of parameters.    
    cs = np.linspace(0, 1.5, 400)
    ds = np.linspace(0, 1.5, 400)
    kdxs = np.linspace(0, 2*np.pi, 50)
    
    amp = centredAmp
    
    maxAs = np.zeros(shape=(cs.shape[0], ds.shape[0]))
    for i, c in enumerate(cs):
        
        for j, d in enumerate(ds):
            
            maxA = 0
            for kdx in kdxs:
                
                A = amp(kdx, c, d).real
                if A > maxA:
                    maxA = A 
            
            maxAs[i, j] = maxA 
    
        #%%
    maxAs[maxAs>1] = np.nan
    
    XS, YS = np.meshgrid(ds, cs)
    plt.contourf(XS, YS, maxAs)
    # plt.plot()
    # plt.contour(XS, YS, maxAs, levels=[1], colors='black', linewidths=0.5)
    # plt.colorbar()
    # plt.plot(ds, 1-2*ds, 'k--', label=r"c+2d = 1")
    # plt.plot([0.5, 0.5], [0, 1], 'k--')
    plt.xlabel("d")
    plt.ylabel("c")
    # plt.legend()
    plt.grid()
    plt.xlim([ds[0], ds[-1]])
    plt.ylim([cs[0], cs[-1]])
    plt.title("Centred Difference")
    plt.show()
                