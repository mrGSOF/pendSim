"""
 * Created on: 19 May 2026
 * Author: Guy Soffer
 * Copyright (C) 2026 Guy Soffer
"""

from Modules.vecLib import scaleV, addV

def FWD_Euler(X, X_dot, dt):
    """ """
    return addV(X, scaleV(X_dot, dt))

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from math import pi, cos, sin
    
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
