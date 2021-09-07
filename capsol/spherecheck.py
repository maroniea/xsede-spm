import numpy as np
import math


from numpy import (pi, sinh, cosh, tanh, arccosh,
                   exp, log10, arctanh)


epsilon_0 = 8.854e-12 # Permittivity


# All uses SI units by default...
def sphereCap(Rtip, d, eps_r, h_sample):
    """To give a dimensionally correct answer, input all values
        in base SI units."""
    alpha = np.arccosh(1 + d/Rtip + h_sample/(eps_r*Rtip))
    return (4 * np.pi * epsilon_0 * Rtip * sum_sinh(alpha))

def sphereMetalCap(Rtip, d):
    """To give a dimensionally correct answer, input all values
        in base SI units.

        Parameters:

        Rtip : Tip radius (m)
        d : Tip-sample separation (m)
        
        Returns: Capacitance (F)"""
    alpha = np.arccosh(1 + d/Rtip)
    return (4 * np.pi * epsilon_0 * Rtip * sum_sinh(alpha))


def coth(x):
    """The hyperpolic cotanget of x."""
    return (1 / tanh(x))

def int_sum_sinh(alpha, eps):
    """Determines the number of terms needed to be within eps
    of the sum, using the integral test."""
    if alpha < 5:
        n_min = np.ceil(2 / alpha * arctanh(exp(-eps * alpha / sinh(alpha))))
    else:
        n_min = 1
    return n_min


def sum_sinh(alpha, eps=1e-8):
    """This calculates the infinite sum, sinh(alpha) / sinh(alpha * n),
    for n = 1 to infinity.
    We manually require at least 4 terms so that the derivative is
    numerically stable. We use math.fsum to give a numerically stable sum."""
    # For alpha greater than 37, the sum is 1 to within double precision,
    # so short-circuit the calculation to avoid overflow errors.
    if alpha > 37:
        return 1
    else:
        summand = lambda n: sinh(alpha) / sinh(alpha * n)

        N_max = max(4, int_sum_sinh(alpha, eps))
        terms = summand(np.arange(1, N_max + 1))
        return math.fsum(terms)

