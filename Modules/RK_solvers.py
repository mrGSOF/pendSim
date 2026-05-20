"""
 * Created on: 19 May 2026
 * Author: Guy Soffer
 * Copyright (C) 2026 Guy Soffer
"""

from Modules.vecLib import scaleV, addV
from Modules.matLib import MxV, inv, scaleV, addV
from Modules.Euler_solvers import FWD_Euler

def RK2(fdiff, X, u, t, dt):
    """ """
    X0 = X
    f1_dx = fdiff(X, u, t)

    f2_X0 = FWD_Euler(X0, f1_dx, 0.5*dt)
    f2_dx = fdiff(f2_X0, u, t)

    X1 = FW_EulerD(X0, f2_dx, dt)
    return X1      

def RK4(fdiff, X, u, t, dt):
    """Y[n+1] = (1/6)*(k1 +2k2 +2*k3 +k4)"""
    t0, h N, = t, dt, len(X)
    
    #k1 = h*f(y, x)
    k1 = scaleV(fdiff(X, u, t), h)
    
    #k2 = h*f(y +0.5*k1, x +0.5*h)  
    X2 = addV(X, [0.5*k1]*N)
    k2 = scaleV(fdiff(X2, u, t +0.5*h), h)

    #k3 = h*f(y +0.5*k2, x +0.5*h)  
    X3 = addV(X, [0.5*k2]*N)
    k3 = scaleV(fdiff(X3, u, t +0.5*h), h)

    #k4 = h*f(y +k3, x +h)  
    X4 = addV(X, [k3]*N)
    k4 = scaleV(fdiff(X4, u, t +h), h)

    #y_new = y +(1/6)*(k1 +2*k2 +2*k3 +k4)
    k = addv(addv(addV(k1, k2), k3), k4)
    return addV(X, scaleV(k, (1/6)))

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matLib import zeros
    import math
    
    assy = Assembly(
        pend_mass_kg=0.1, pend_rod_m=1.0,
        pend_friction=0.01, pend_viscosity=0.01,
        pend_omega_rps=0.0, pend_theta_rad=0.0,
        cart_mass_kg=1.0,
        cart_friction=0.1, cart_viscosity=0.2,
        cart_vel_mps=0.0, cart_pos_m=0.0,
        encoderLines = 4096, g_mps2=9.81
        ).print()

    dt = 0.01
    STEPS  = int(10/dt)
    force  = [0.01]*STEPS
    theta  = [0]*STEPS
    omega  = [0]*STEPS
    pos    = [0]*STEPS
    count  = [0]*STEPS
    time   = [0]*STEPS

    for frame in range(0, STEPS):
#        force[frame] = ctrl.step()            #< Step the controller
        assy.step( force=force[frame], dt=dt ) #< Step the cart+pendulum assembly

        ### Store the results
        time[frame]  = assy.time
        theta[frame] = assy.pend.theta_rad
        omega[frame] = assy.pend.omega_rps
        pos[frame]   = assy.cart.pos_m
        #count[frame]  = assy.getEncoderCount()/1024
    plt.plot(time, theta)
    plt.plot(time, pos)
    #plt.plot(time, count)
    plt.show()
