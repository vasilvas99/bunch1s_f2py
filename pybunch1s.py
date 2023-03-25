"""
    pybunch1s is a high-level pythonic interface module for the bunch1s models. Requires the compiled fortran bunch1s module
    to be availabe.

    Download and extract the pre-compiled fortran code for your platform from:

    Linux: https://nightly.link/vasilvas99/bunch1s_f2py/workflows/compile-bunch1s/main/build_unix.zip?h=5e973cebe0192e1902ec70ffb6946e18b78a24b2
    Windows: https://nightly.link/vasilvas99/bunch1s_f2py/workflows/compile-bunch1s/main/build_windows.zip?h=5e973cebe0192e1902ec70ffb6946e18b78a24b2

    and place it next to this file.
"""
try:
    import bunch1s as b1s
except:
    raise ModuleNotFoundError(
        """
    pybunch1s is a high-level interface module. It requires the fortran source for bunch1s to be compiled with f2py first.
    To obtain the pre-compiled fortran module download it for your platform:
    
    Linux: https://nightly.link/vasilvas99/bunch1s_f2py/workflows/compile-bunch1s/main/build_unix.zip?h=5e973cebe0192e1902ec70ffb6946e18b78a24b2
    Windows: https://nightly.link/vasilvas99/bunch1s_f2py/workflows/compile-bunch1s/main/build_windows.zip?h=5e973cebe0192e1902ec70ffb6946e18b78a24b2
    
    Or compile it manually with Cmake by following the instructions at: https://github.com/vasilvas99/bunch1s_f2py#compiling-manually

    Once you have obtained the .pyd (Windows) or .so (Linux) compiled bunch1s module, put in the same directory as this module (pybunch1s).
    """
    )

import numpy as np


def gpmm2(y: np.ndarray, bdef: float, p1: float, p2: float) -> np.ndarray:
    """
    13.09.02
    corresponds to MMII, non-dimensionalized!

    Last correction of the model:

    It is changed to be in analogy with LW2 with two constants - K and U

    and two exponents - p and n, previous        versions of the model in the new

    frame should be classified in p = 0 and n_old = n_new - 1

    (c) V. Tonchev, 04.01.2012

    Args:
        y (np.ndarray): step positions
        bdef (float): bunch definition
        p1 (float): interaction exponent
        p2 (float): interaction exponent

    Returns:
        np.ndarray: step velocities
    """
    return b1s.gpmm2(y, bdef, p1, p2)


def g1smm(y: np.ndarray, bdef: float, p1: float, p2: float) -> np.ndarray:
    """

    One-sided MM0 model non-dimensionalized!

    (c) V. Tonchev, 15.01.2012

    Args:
        y (np.ndarray): step positions
        bdef (float): bunch definition
        p1 (float): interaction exponent
        p2 (float): interaction exponent

    Returns:
        np.ndarray: step velocities
    """
    return b1s.g1smm(y, bdef, p1, p2)


def g1slw(y: np.ndarray, b: float, U: float, n: float) -> np.ndarray:
    """
    Model: Popkov, Krug (LW), PRB 73, 235430 (2006)

    Args:
        y (np.ndarray): step positions
        b (np.ndarray): TODO
        U (np.ndarray): TODO
        n (np.ndarray): interaction order (exponent) (en-1 in fortran)

    Returns:
        np.ndarray: step velocities
    """
    par = np.array([b, U, 0, 0, n])
    return b1s.g1slw(y, par)


def gkrug(y: np.ndarray, b: float, U: float, n: float) -> np.ndarray:
    """
    Model: Popkov, Krug (LW), PRB 73, 235430 (2006)

    Args:
        y (np.ndarray): step positions
        b (np.ndarray): TODO
        U (np.ndarray): TODO
        n (np.ndarray): interaction order (exponent) (en-1 in fortran!!)

    Returns:
        np.ndarray: step velocities
    """
    par = np.array([b, U, 0, 0, n])
    return b1s.gkrug(y, par)


def g_pk2(y: np.ndarray, K: float, U: float, ro: float, n: float) -> np.ndarray:
    """
    Evolves from the Model of Popkov, Krug, PRB 73, 235430 (2006), LW2

    the stabilization part is the same but the destabilization part is the stabilization one with an opposite sign

    and evenually different power.

    Args:
        y (np.ndarray): step positions
        K (float): TODO
        U (float): TODO
        ro (float): TODO
        n (float): interaction order (exponent) (en-1 in fortran!!)

    Returns:
        np.ndarray: step velocities
    """
    par = np.array([K, U, 0, ro, n])
    return b1s.g_pk2(y, par)


def g_mm0(y: np.ndarray, b: float, U: float, n: float) -> np.ndarray:
    """
    r - fixed to 1 !!!!
    corresponds to MMI as introduced by VT in 2002 but slightly changed
    to meet the analogy with the model of Popkov and Krug thus

    it still has 2 parameters - b and U, like in the PK model         plus

    the two powers r and n, note the difference - here

    in the stabilization part the terrace widths are raised to (n+1)

    (c) VT, May 2008, Lexington KY
    Args:
        y (np.ndarray): step positions
        b (float): TODO
        U (float): TODO
        n (float): interaction order (exponent) (n-1 in fortran!!)

    Returns:
        np.ndarray: step velocities
    """
    par = np.array([b, U, 0, 0, n])
    return b1s.g_mm0(y, par)


def g_mm1(y: np.ndarray, alpha: float, beta: float, rho: float, n: float) -> np.ndarray:
    """
    corresponds to MMI

    Args:
        y (np.ndarray): step positions
        alpha (float): attraction parameter
        beta (float): repulsion parameter
        rho (float): TODO
        n (float): interaction order (exponent) (n-1 in fortran!!)

    Returns:
        np.ndarray: step velocities
    """
    par = [rho, beta, alpha, 0, n]
    return b1s.g_mm1(y, par)


def gise2(
    y: np.ndarray, gamma: np.ndarray, l0l: float, dpl: float, dml: float, n: float
) -> np.ndarray:
    """
    TODO: Summary

    Args:
        y (np.ndarray): step positions
        gamma (np.ndarray): TODO
        l0l (float): initia vicinal distance
        dpl (float): TODO
        dml (float): TODO
        n (float): interaction order (exponent) (n-1 in fortran!!)

    Returns:
        np.ndarray: step velocities
    """
    par = np.ndarray(gamma, l0l, dpl, dml, n)
    return b1s.gise2(y, par)
