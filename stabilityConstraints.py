# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 13:34:52 2024

@author: guest997558
"""

import numpy as np

import plotters

def upwindAmp(kdx, c, d):
    A = 1 + (2*d+c)*(np.cos(kdx) - 1) - 1j*c*np.sin(kdx)
    return A*A.conjugate()

def centredAmp(kdx, c, d):
    A = 1 + 2*d*(np.cos(kdx)-1) - 1j*c*np.sin(kdx)
    return A*A.conjugate()

def findMaxAs(cs, ds, kdxs, amp):
    """ 
    """
    
    maxAs = np.zeros(shape=(cs.shape[0], ds.shape[0]))
    
    # Consider a range of Courant numbers.
    for i, c in enumerate(cs):
        
        # Consider a range of diffusion coefficients.
        for j, d in enumerate(ds):
            
            # Find the maximum amplification factor for a range of kdx.
            maxA = 0
            for kdx in kdxs:
                A = amp(kdx, c, d).real
                if A > maxA:
                    maxA = A 
            maxAs[i, j] = maxA 
    return maxAs
        
def runNumericalStabilityConstrants():
    """ 
    """
    
    # These are hard coded here.
    casesUpwind = [(0.8, 0.4), (0.5, 0.1)]
    casesCentred = [(0.7, 0.1), (0.5, 0.3), (0.4, 0.6)]
    
    # Create ranges of parameters.    
    cs = np.linspace(0, 1.5, 200)
    ds = np.linspace(0, 1.5, 200)
    kdxs = np.linspace(0, 2*np.pi, 50)
    
    # Stability constraints for centred difference.
    maxACentred = findMaxAs(cs, ds, kdxs, centredAmp)
    plotters.plotAmplifcationFactors(cs, ds, maxACentred, casesCentred, "Centred Difference",
                            figname="centredDifferenceStabilityConstraint.pdf", saveFig=True)

    # Stability constraints for upwind.
    maxAUpwind = findMaxAs(cs, ds, kdxs, upwindAmp)
    plotters.plotAmplifcationFactors(cs, ds, maxAUpwind, casesUpwind, "Upwind",
                            figname="upwindStabilityConstraint.pdf", saveFig=True)

if __name__ == "__main__":
    
    runNumericalStabilityConstrants()
        