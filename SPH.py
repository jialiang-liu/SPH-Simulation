#!/usr/bin/env python3

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import math
import numpy as np

# Parameters for the glass
radius = 150
height = 500

# Parameters for particles
n = 1000
mass = 0.25
density = 1
h = 10  # Kernal radius
visc = 1  # Viscosity constant
dt = 0.005  # Integration timestep
damp = -0.5

poly6 = mass * 315 / (65 * math.pi * (h ** 9))
spiky = visc * mass * 45 / (math.pi * (h ** 6))
lapla = visc * mass * 45 / (math.pi * (h ** 6))

# Store values
location = np.zeros((n, 2))
velosity = np.zeros((n, 2))
force = np.zeros((n, 2))
rho = np.zeros((n))
pressure = np.zeros((n))