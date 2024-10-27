# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 22:05:37 2024

@author: guest997558
"""

import numpy as np
from scipy.integrate import odeint

class ReferenceSolution:
    """
    Solves the 1D Burgers' equation using the FFT for spatial derivatives.
    Equation: ∂u/∂t + μ * u * ∂u/∂x = ν * ∂²u/∂x²
    """

    def __init__(self, params, grid, endtime, nt):
        """
        """
        self.params = params
        self.k = 2 * np.pi * np.fft.fftfreq(grid.nx, d=grid.dx)  # Wavenumber discretization
        self.T = np.linspace(0, endtime, nt)  # Time array
        self.phi = None  # Placeholder for the solution to be computed

    def evaluate(self, phi0):
        """
        """
        args = (self.k, self.params.mu, self.params.nu)
        # Solve the system using odeint and transpose to shape (nt, nx)
        self.phi = odeint(self.burgerSystem, phi0, self.T, args=args, mxstep=5000).T

    def burgerSystem(self, u, t, k, mu, nu):
        """
        """
        # Spatial derivatives in Fourier domain
        u_hat = np.fft.fft(u)
        u_hat_x = 1j * k * u_hat
        u_hat_xx = -k**2 * u_hat

        # Transform back to spatial domain
        u_x = np.fft.ifft(u_hat_x).real
        u_xx = np.fft.ifft(u_hat_xx).real

        # ODE resolution in the spatial domain
        u_t = -mu * u * u_x + nu * u_xx
        return u_t