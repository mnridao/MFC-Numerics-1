# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 19:48:13 2024

@author: guest997558
"""

from stabilityConstraints import runNumericalStabilityConstrants
from stabilityAnalysis import runStabilityExperiment

if __name__ == "__main__":
        
    # Find numerical stability constraints for each advection scheme.
    runNumericalStabilityConstrants()
    
    # Perform a stability analysis for each advection scheme.
    runStabilityExperiment()
    
    # Perform an order of convergence analysis for each scheme.
    # ...