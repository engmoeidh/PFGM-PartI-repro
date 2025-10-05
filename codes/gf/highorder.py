import numpy as np

def deriv4_central(phi, h):
    """Fourth-order centered first derivative; interior only."""
    d = np.zeros_like(phi)
    d[2:-2] = ( -phi[:-4] + 8*phi[1:-3] - 8*phi[3:-1] + phi[4:] ) / (12*h)
    return d

def deriv4_one_sided_left(phi, h):
    # 4th-order forward difference for first point
    return (-25*phi[0] + 48*phi[1] - 36*phi[2] + 16*phi[3] - 3*phi[4])/(12*h)

def deriv4_one_sided_right(phi, h):
    # 4th-order backward difference for last point
    n=len(phi)
    return (25*phi[n-1] - 48*phi[n-2] + 36*phi[n-3] - 16*phi[n-4] + 3*phi[n-5])/(12*h)

def dphi_dx_4th(phi, h):
    d = deriv4_central(phi, h)
    d[0]  = deriv4_one_sided_left(phi, h)
    d[1]  = ( -3*phi[0] - 10*phi[1] + 18*phi[2] - 6*phi[3] + phi[4])/(12*h)  # 3rd order near boundary
    d[-1] = deriv4_one_sided_right(phi, h)
    d[-2] = ( 3*phi[-1] + 10*phi[-2] - 18*phi[-3] + 6*phi[-4] - phi[-5])/(12*h)
    return d

def simpson(y, x):
    """Composite Simpson rule (even number of intervals required)."""
    if len(x)%2==0:  # make odd number of points by dropping last if needed
        x=x[:-1]; y=y[:-1]
    h = x[1]-x[0]
    S = y[0] + y[-1] + 4*np.sum(y[1:-1:2]) + 2*np.sum(y[2:-2:2])
    return S * h/3.0
