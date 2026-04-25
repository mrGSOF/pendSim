"""
 * Created on: 9 April 2026
 * Author: Guy Soffer
 * Copyright (C) 2026 Guy Soffer
"""

from math import pi, cos, sin
#from vecLib import subV, scaleV, crossV3, polarV2
#from matLib import MxV, rotateV2, getM
#
#def crossV2(V2a, V2b) -> float:
#    return (crossV3(list(V2a) +[0], list(V2b) +[0]))[2]

_2pi = 2*pi
def sign(val) -> int:
    if val < 0.0:
        return -1
    return 1

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
                 acc_mps2=0.0, vel_mps=0.0, pos_m=0.0):
        self.mass_kg = mass_kg
        self.friction = friction
        self.viscosity = viscosity
        self.setState(
            pos_m,    #< Angle
            vel_mps,  #< Angular velocity
            acc_mps2  #< Angular acceleration
        )

    def setState(self, pos_m, vel_mps, acc_mps2) -> None:
        self.acc_mps2 = acc_mps2
        self.vel_mps  = vel_mps
        self.pos_m    = pos_m

    def getState(self) -> list:
        return [self.pos_m,
                self.vel_mps,
                self.acc_mps2
                ]

    def __str__(self) -> str:
        s = "## CART STATE ##\n"
        s += "acc_mps2 : %.4f\n"%self.acc_mps2
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
                                  0.0, pend_omega_rps, pend_theta_rad)
        
        self.cart      = Cart(cart_mass_kg, cart_friction, cart_viscosity,
                              0.0, cart_vel_mps, cart_pos_m)
        
        self.encoder   = Encoder(counts=encoderLines)

    def step(self, force, dt):
        self.time += dt
        f = force - self.cart.viscosity*self.cart.vel_mps
        mp = self.pend.mass_kg
        p = self.pend.getState()
        mc = self.cart.mass_kg
        c = self.cart.getState()
        r = self.pend.radius_m
        g = self.g_mps2
        Sx = sin(p[0])
        Cx = cos(p[0])
        w2 = p[1]**2
        D  = mc +mp*(Sx**2)
        
        #acc = (f -mp*Sx*(r*(p[1]**2) -g*Cx)) / (mc +mp*(Sx**2))
        acc = (f -mp*Sx*(r*w2 -g*Cx)) / D
        vel = c[1] +acc*dt
        pos = c[0] +vel*dt
        self.cart.setState(pos, vel, acc)        
        
        #omega_dot = (Cx*(f -(p[1]**2)*mp*r*Sx) +g*Sx*(mp+mc)) / (r*(mc +mp*(Sx**2)))
        omega_dot = (Cx*(f -w2*mp*r*Sx) +g*Sx*(mp+mc)) / (r*D)
        omega = p[1] +omega_dot*dt
        theta = p[0] +omega*dt
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
