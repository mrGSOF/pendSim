"""
 * Created on: 9 April 2026
 * Author: Guy Soffer
 * Copyright (C) 2026 Guy Soffer
"""

from math import pi, cos, sin
from Modules.vecLib import scaleV, addV#, subV, crossV3, polarV2
from Modules.matLib import MxV, inv3#, rotateV2, getM

_2pi = 2*pi
def sign(val) -> int:
    if val < 0.0:
        return -1
    return 1

def _FW_EulerInt(X, X_dot, dt):
    return addV(X, scaleV(X_dot, dt))

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
        omega_dot_rps2=0.0, omega_rps=0.0, theta_rad=0.0):
        
        self.mass_kg  = mass_kg
        self.radius_m = radius_m
        self.setState(
            theta_rad,      #< Angle
            omega_rps,      #< Angular velocity
            omega_dot_rps2  #< Angular acceleration
        )
        
    def setState(self, theta_rad, omega_rps, omega_dot_rps2) -> None:
        self.omega_dot_rps2 = omega_dot_rps2
        self.omega_rps      = omega_rps
        self.theta_rad      = theta_rad

    def getState(self) -> list:
        return [self.theta_rad,
                self.omega_rps,
                self.omega_dot_rps2
                ]

    def __str__(self) -> str:
        s = "## ROTOR STATE ##\n"
        s += "omega_dot_rps2: %.4f\n"%self.omega_dot_rps2
        s += "omega_rps     : %.4f\n"%self.omega_rps
        s += "theta-rad     : %.4f\n"%self.theta_rad
        return s

    def print(self) -> None:
        print(self.__str__())

class Cart():
    def __init__(self, mass_kg, friction, viscosity,
                 accX_mps2=0.0, velX_mps=0.0, posX_m=0.0,
                 accY_mps2=0.0, velY_mps=0.0, posY_m=0.0,
                 ):
        self.mass_kg = mass_kg
        self.friction = friction
        self.viscosity = viscosity
        self.setState(
            posX_m,    #< Position-X
            velX_mps,  #< Velocity-X
            accX_mps2, #< Acceleration-X
            posY_m,    #< Position-Y
            velY_mps,  #< Velocity-Y
            accY_mps2, #< Acceleration-Y
        )

    def setState(self, posX_m, velX_mps, accX_mps2, posY_m, velY_mps, accY_mps2) -> None:
        self.accX_mps2 = accX_mps2
        self.velX_mps  = velX_mps
        self.posX_m    = posX_m
        self.accY_mps2 = accY_mps2
        self.velY_mps  = velY_mps
        self.posY_m    = posY_m

    def getState(self) -> list:
        return [self.posX_m,
                self.velX_mps,
                self.accX_mps2,
                self.posY_m,
                self.velY_mps,
                self.accY_mps2
                ]

    def __str__(self) -> str:
        s = "## CART STATE ##\n"
        s += "accX_mps2 : %.4f\n"%self.accX_mps2
        s += "velX_mps  : %.4f\n"%self.velX_mps
        s += "posX_m    : %.4f\n"%self.posX_m
        s += "accY_mps2 : %.4f\n"%self.accY_mps2
        s += "velY_mps  : %.4f\n"%self.velY_mps
        s += "posY_m    : %.4f\n"%self.posY_m
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
                                  0.0, pend_omega_rps, pend_theta_rad)
        
        self.cart      = Cart(cart_mass_kg, cart_friction, cart_viscosity,
                              0.0, cart_vel_mps, cart_pos_m)
        
        self.encoder   = Encoder(counts=encoderLines)

    def step(self, forceX, forceY, dt):
        self.time += dt
        Fx = forceX -self.cart.viscosity*self.cart.velX_mps
        Fy = forceY -self.cart.viscosity*self.cart.velY_mps
        #self.FW_Euler(Fx, Fy, dt)
        self.RK2(Fx, Fy, dt)
        
    def _diff(self, u, cartState, pendState):
        Fx, Fy = u
        mp = self.pend.mass_kg
        mc = self.cart.mass_kg
        g = self.g_mps2
        r = self.pend.radius_m
        Vx = cartState[1]
        Vy = cartState[4]

        Sx = sin(pendState[0])
        Cx = cos(pendState[0])
        if (pendState[1] > 100):
            pendState[1] = 100
        elif (pendState[1] < -100):
            pendState[1] = -100
        w = pendState[0]
        w2 = w**2

        C = [[mc+mp, 0.0, -mp*r*Cx],
             [0.0, mc+mp, mp*r*Sx],
             [-Cx, Sx, r]
            ]

        D = [Fx -mp*r*Sx*w2,
             Fy -mp*r*Cx*w2,
             Sx*Vx*(1-w) -Cx*Vy*(1-w) +g*Sx
             ]

##        C = [[mc+mp, 0.0, -mp*r*Cx],
##             [0.0, mc+mp, 0.0],
##             [-Cx, 0.0, r]
##            ]
##
##        D = [Fx -mp*r*Sx*w2,
##             Fy,
##             g*Sx
##             ]

        accX, accY, omega_dot = MxV(inv3(C), D)
        return (accX, accY, omega_dot)
        
    def FW_Euler(self, forceX, forceY, dt):
        cartState = self.cart.getState()
        pendState = self.pend.getState()
        A1 = self._diff((forceX, forceY),
                        cartState,
                        pendState
                       )
        V0 = (cartState[1], cartState[4], pendState[1])
        V1 = _FW_EulerInt(V0, A1, dt)

        X0 = (cartState[0], cartState[3], pendState[0])
        X1 = _FW_EulerInt(X0, V1, dt)

        accX, accY, omega_dot = A1
        velX, velY, omega = V1
        posX, posY, theta = X1
        
        self.cart.setState(posX, velX, accX, posY, velY, accY)        
        self.pend.setState(theta, omega, omega_dot)        

    def RK2(self, forceX, forceY, dt):
        u = (forceX, forceY)
        f1_cartState = self.cart.getState()
        f1_pendState = self.pend.getState()
        f1_A1 = self._diff(u,
                           f1_cartState,
                           f1_pendState
                          )
        
        f1_V0 = (f1_cartState[1], f1_cartState[4], f1_pendState[1])

        f2_V1 = _FW_EulerInt(f1_V0, f1_A1, 0.5*dt)
        f1_X0 = (f1_cartState[0], f1_cartState[3], f1_pendState[0])
        f2_X1 = _FW_EulerInt(f1_X0, f2_V1, 0.5*dt)

        accX, accY, omega_dot = f1_A1
        velX, velY, omega     = f2_V1
        posX, posY, theta     = f2_X1
        
        f2_cartState = [posX, velX, accX, posY, velY, accY]
        f2_pendState = [theta, omega, omega_dot]
        f2_A1 = self._diff(u,
                           f2_cartState,
                           f2_pendState
                          )
        V1 = _FW_EulerInt(f1_V0, f2_A1, dt)
        X1 = _FW_EulerInt(f1_X0, V1, dt)
        
        accX, accY, omega_dot = f2_A1
        velX, velY, omega     = V1
        posX, posY, theta     = X1
        self.cart.setState(posX, velX, accX, posY, velY, accY)        
        self.pend.setState(theta, omega, omega_dot)        

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
    posX   = [0]*STEPS
    posY   = [0]*STEPS
    count  = [0]*STEPS
    time   = [0]*STEPS

    for frame in range(0, STEPS):
#        force[frame] = ctrl.step()            #< Step the controller
        assy.step( forceX=force[frame], forceY=0.0, dt=dt ) #< Step the cart+pendulum assembly

        ### Store the results
        time[frame]  = assy.time
        theta[frame] = assy.pend.theta_rad
        omega[frame] = assy.pend.omega_rps
        posX[frame]   = assy.cart.posX_m
        posY[frame]   = assy.cart.posY_m
        #count[frame]  = assy.getEncoderCount()/1024
    plt.plot(time, theta)
    plt.plot(time, posX)
    plt.plot(time, posY)
    #plt.plot(time, count)
    plt.show()
