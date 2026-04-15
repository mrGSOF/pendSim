# Pendulum on Cart Simulation
Simple physics simulation of pendulum on a cart. Demostrating movement of pendulum and cart.

## Theory of operation
Some excellent tutorial on youtube for how to derive the equations of motion:\
https://www.youtube.com/watch?v=5qJY-ZaKSic

https://www.youtube.com/watch?v=qjhAAQexzLg

https://www.youtube.com/watch?v=CAHIqiSthaM

## Why?
1. To be an interactive tool and support tutorials.
2. To interact with the controller during development (Software In The Loop). This eliminates the needs for any hardware during early stages of controller develop.

## What does it simulate / demonstrate
TBD

### Free body cart with pendulum
![alt text](./figures/cartPendulum.gif "Free body cart with pendulum)")

## Future plans
TBD

Requires Python to run the PND_model (inclues dedicated unit-test). For visual support install pyGame and GSOF_Cockpit as well.

http://python.org/

http://www.pygame.org

https://github.com/mrGSOF/GSOF_Cockpit

## Running instructions
- Install requirements `pip install -r requirements.txt`
- Clone and install GSOF_Cockpit (`pip setup.py`)
- Clone bldcSim
- run `python Demo_PNDL_simple.py`

Interactive operation is supported using the mouse.
