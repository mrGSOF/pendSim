#!/usr/bin/python
"""
 * Demo_PNDL_openloop.py
 * Created on: 13 April 2026
 * Author:     Guy Soffer (gsoffer@yahoo.com)
 * Copyright (C) 2026 Guy Soffer
"""
import os, time
from math import pi
from Modules.PND_model_XY import Assembly
#from Modules.PND_model import Assembly
from Modules.Mouse import Mouse
from GSOF_Cockpit.Generic import Map as MAP

from GSOF_Cockpit.GraphicsLib import imageLoad, getScreen, init, fillScreen, update
from GSOF_Cockpit.Text import Text
from GSOF_Cockpit import Pygame_Colors as COLOR
from GSOF_Cockpit.Clock_base import Clock

def wait(t):
    while t>time.time():
        pass
    return

# Initialise screen.
BG_color = COLOR.DARK
screen_size=(960,640)
pos = (0, 0)
init()
screen = getScreen(screen_size)
fillScreen( screen, COLOR.WHITE )
   
# Initialise Dials.
path = './skin'
background = Text(  screen=screen, pos=pos, size=screen_size, color=BG_color, name='' )

assy = Assembly(
    pend_mass_kg=0.1, pend_rod_m=0.5,
    pend_friction=0.02, pend_viscosity=0.01,
    pend_omega_rps=0.0, pend_theta_rad=-0.1,
    cart_mass_kg=0.2,
    cart_friction=0.1, cart_viscosity=3.0,
    cart_vel_mps=0.0, cart_pos_m=0.0,
    encoderLines = 4096, g_mps2=9.81
    ).print()

view = MAP.Map(
               screen, pos=pos,
               size=screen_size,
               kp=0.9,
               bodyImage   = imageLoad("%s/Background_960x640.png"%path),
               #mapImage    = imageLoad("%s/Grid_BackgroundWhite300x300.png"%path),
               markerImage = imageLoad("%s/Assy_640x640.png"%path),
              )

ctrl = Mouse(screen_size)
clock = Clock()

fps = 20
dt = 1/fps
decimation = 50
_dt = 1/(fps*decimation)
radToDeg = 180/pi
M_TO_PXL = 300*assy.pend.radius_m

cartPosX_pxl = M_TO_PXL*assy.cart.posX_m
cartPosY_pxl = M_TO_PXL*assy.cart.posY_m
print("DT: %1.4f seconds"%(dt))
print("FPS: %d frames per second"%(fps))
print("DEC: %1.2f calcs per draw"%(decimation))

T0 = time.time()
t0 = T0
while True:
    ###Loop to update gauges
    ct = ctrl.getMouse()
    for i in range(0,decimation):
        forceX = 0.5*(ct["posX"] -cartPosX_pxl) 
        forceY = 0.5*(ct["posY"] -cartPosY_pxl) 
        assy.step(forceX, forceY, _dt)
        cartPosX_pxl = M_TO_PXL*assy.cart.posX_m
        cartPosY_pxl = M_TO_PXL*assy.cart.posY_m
    
    view.update(
        x=cartPosX_pxl,
        y=cartPosY_pxl,
        deg=-radToDeg*assy.pend.theta_rad
        )
    
    view.draw()
    update()
    T0 += dt
    wait(T0) #< More accurate than clock.tick(Fs=fps) or clock.tick_busy_loop(Fs=fps)
    print("%1.4f"%(time.time() -T0))
