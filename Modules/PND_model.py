"""
 * Created on: 9 April 2026
 * Author: Guy Soffer
 * Copyright (C) 2026 Guy Soffer
"""

from math import pi, cos, sin
from Modules.vecLib import scaleV, addV
from Modules.matLib import MxV, inv, scaleV, addV

_2pi = 2*pi
def sign(val) -> int:
    if val < 0.0:
        return -1
    return 1

def _FW_EulerInt(X, X_dot, dt):
    """ """
    return addV(X, scaleV(X_dot, dt))

def _RK2(fdiff, X, dt):
    """ """
    Fx, f1_cartState, f1_pendState = X
    X0 = f1_cartState +f1_pendState
    f1_dx = fdiff(Fx,
                  f1_cartState,
                  f1_pendState
                 )

    f2_X0 = _FW_EulerInt(X0, f1_dx, 0.5*dt)
    f2_dx = fdiff(Fx,
                  f2_X0[0:2],
                  f2_X0[2:4]
                 )

    X1 = _FW_EulerInt(X0, f2_dx, dt)
    return X1      

##def _RK4(fdiff, X, dt):
##    """ """
##    return

class Encoder():
    def __init__(self, counts):

        self.counts = counts
        self.radToCounts = counts/_2pi

    def getCount(self, rad) -> int:
        return (rad%_2pi)*self.radToCounts
        
        
class Pendulum():
    def __init__(
        self,
        mass_kg, radius_m,
        friction=0.01, viscosity=0.01,
        omega_rps=0.0, theta_rad=0.0):
        
        self.mass_kg  = mass_kg
        self.radius_m = radius_m
        self.setState(
            #omega_dot_rps2, #< Angular acceleration
            omega_rps,      #< Angular velocity
            theta_rad       #< Angle
        )
        
    def setState(self, omega_rps, theta_rad) -> None:
        #self.omega_dot_rps2 = omega_dot_rps2
        self.omega_rps      = omega_rps
        self.theta_rad      = theta_rad

    def getState(self) -> list:
        return [#self.omega_dot_rps2,
                self.omega_rps,
                self.theta_rad,                
                ]

    def __str__(self) -> str:
        s = "## ROTOR STATE ##\n"
        #s += "omega_dot_rps2: %.4f\n"%self.omega_dot_rps2
        s += "omega_rps     : %.4f\n"%self.omega_rps
        s += "theta-rad     : %.4f\n"%self.theta_rad
        return s

    def print(self) -> None:
        print(self.__str__())

class Cart():
    def __init__(self, mass_kg, friction, viscosity,
                 vel_mps=0.0, pos_m=0.0):
        self.mass_kg = mass_kg
        self.friction = friction
        self.viscosity = viscosity
        self.setState(
            #acc_mps2  #< Angular acceleration
            vel_mps,  #< Angular velocity
            pos_m,    #< Angle
        )

    def setState(self, vel_mps, pos_m ) -> None:
        #self.acc_mps2 = acc_mps2
        self.vel_mps  = vel_mps
        self.pos_m    = pos_m

    def getState(self) -> list:
        return [#self.acc_mps2,
                self.vel_mps,
                self.pos_m
                ]

    def __str__(self) -> str:
        s = "## CART STATE ##\n"
        #s += "acc_mps2 : %.4f\n"%self.acc_mps2
        s += "vel_mps  : %.4f\n"%self.vel_mps
        s += "pos_m    : %.4f\n"%self.pos_m
        return s
    
    def print(self):
        print(self.__str__())
        return self
        
class Assembly():
    def __init__(
        self,
        pend_mass_kg=0.1, pend_rod_m=1.0,
        pend_friction=0.01, pend_viscosity=0.01,
        pend_omega_rps=0.0, pend_theta_rad=0.0,
        cart_mass_kg=1.0,
        cart_friction=0.1, cart_viscosity=0.2,
        cart_vel_mps=0.0, cart_pos_m=0.0,
        encoderLines = 4096, g_mps2=9.81):
        self.time = 0.0
        self.g_mps2 = g_mps2
        self.pend      = Pendulum(pend_mass_kg, pend_rod_m,
                                  pend_friction, pend_viscosity,
                                  pend_omega_rps, pend_theta_rad)
        
        self.cart      = Cart(cart_mass_kg, cart_friction, cart_viscosity,
                              cart_vel_mps, cart_pos_m)
        
        self.encoder   = Encoder(counts=encoderLines)

    def step(self, forceX, forceY, dt):
        self.time += dt
        Fx = forceX
        #self.FW_Euler(Fx, dt)
        self.RK2(Fx, dt)

    def _diff(self, u, cartState, pendState):
        Fx = u
        Vx = cartState[0] #<vel_mps
        v = self.cart.viscosity
        mp = self.pend.mass_kg
        p = self.pend.getState()
        mc = self.cart.mass_kg
        r = self.pend.radius_m
        g = self.g_mps2
        Sx = sin(p[1])
        Cx = cos(p[1])
        w = p[0]
        w2 = w*w

        C = [[mc+mp, 0.0, -mp*r*Cx, 0.0],
             [0.0,   1.0,    0.0,   0.0],
             [-Cx,   0.0,     r,    0.0],
             [0.0,   0.0,    0.0,   1.0]
            ]

        D = [Fx -mp*r*Sx*w2 -v*Vx,
             Vx,
             g*Sx,
             w]

        invC = inv(C)

        accX, velX, omega_dot, omega = MxV(invC, D)

        return (accX, velX, omega_dot, omega)
        
    def FW_Euler(self, forceX, dt):
        cartState = self.cart.getState()
        pendState = self.pend.getState()
        dX = self._diff(forceX,
                        cartState,
                        pendState
                       )
        X0 = cartState + pendState
        X1 = _FW_EulerInt(X0, dX, dt)

        velX, posX  = X1[0:2]
        omega,theta = X1[2:4]
        
        self.cart.setState(velX, posX)        
        self.pend.setState(omega, theta)        

    def RK2(self, forceX, dt):
        X = (forceX, self.cart.getState(), self.pend.getState())
        X1 = _RK2(self._diff, X, dt)
        cartState = X1[0:2]
        pendState = X1[2:4]
        self.cart.setState(*cartState)        
        self.pend.setState(*pendState)        

    def getPend_rad(self):
        return self.rotor.theta_rad
    
    def getEncoderCount(self) -> int:
        return self.encoder.getCount(self.rotor.theta_rad)

    def __str__(self):
        s = str(self.pend) +"\n"
        s += str(self.cart) +"\n"
        s += "Gravity (m/s2): %.4f"%(self.g_mps2)
        s += "\n"
        return s

    def print(self):
        print(self.__str__())
        return self

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from matLib import zeros
    
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
