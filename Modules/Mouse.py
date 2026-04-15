"""
 * Data.py
 * Created on: 6 Jan 2025
 * Author:     Guy Soffer
 * Copyright (C) 2025 Guy Soffer
"""

from math import pi, atan2
import sys
import pygame
from GSOF_Cockpit.GraphicsLib import getMouse

radToDeg = 180/pi
class Mouse():
    """Data source to drive gauges screen"""
    def __init__(self, screen_size):
        self.zeroPoint = (int(screen_size[0]/2), int(screen_size[0]/2))
        self.cnt = 0
        self.pos = (0,0)
        self.pos_Z1 = self.pos
        self.dir    = 0.0
        self.dir_Z1 = self.dir

    def getMouse(self) -> dict:
        """Generate and return new set of data"""
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                print('Exiting....')
                sys.exit()   # end program.

        # We have data.
        self.cnt += 1

        # Sample the mouse position
        self.pos_Z1 = self.pos
        self.pos = (getMouse())["pos"]

        velX = self.pos[0] -self.pos_Z1[0]
        velY = self.pos[1] -self.pos_Z1[1]
        self.dir_Z1 = self.dir
        self.dir = atan2(-velY, velX)

        return {'cnt': self.cnt,
                'posX':(self.pos[0] -self.zeroPoint[0]),
                'posY':(self.pos[1] -self.zeroPoint[1]),
                'velX': velX,
                'velY': velY,
                'deg': -radToDeg*self.dir +90,
               }
