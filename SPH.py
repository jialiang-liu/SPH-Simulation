#!/usr/bin/env python3

'''from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *'''
import math
import numpy as np
from spatialhash import SpatialHash
import matplotlib.pyplot as plt
from matplotlib import animation

# Parameters for the glass
diameter = 300
height = 500

#ylim = [50, 470]
'''def Xlim(y):
	if y > height:
		return [10, diameter - 10]
	temp = y / height * 50
	return [10 + temp, diameter - 10 - temp]'''
#xlim = [50, 250]
xlim = [10, diameter - 10]
ylim = [0, height]

# Parameters for particles
n = 1000
mass = 1
k = 20
rho_0 = 1
h = 5  # Kernal radius
visc = 1  # Viscosity constant
dt = 0.01  # Integration timestep
damp = -0.5
fg = np.array([0, -0.1])

# Kernels
# For density...
def poly6():
	global mass, h
	return mass * 315 / (64 * math.pi * (h ** 9))
# For pressure...
def spiky():
	global mass, h
	return -mass *  45 / (math.pi * (h ** 6))
# For viscosity...
def laplacian():
	global mass, h
	return mass * 45 / (math.pi * (h ** 6))

# Store values
location = np.zeros((n, 2))
velosity = np.zeros((n, 2))
force = np.zeros((n, 2))
rho = np.zeros((n))
pressure = np.zeros((n))
hashmap = SpatialHash(diameter, height, h)

def bound(i):
	global location, velosity, xlim, ylim, damp
	if location[i, 0] < xlim[0]:
		velosity[i, 0] *= damp
		location[i, 0] = xlim[0]
	elif location[i, 0] > xlim[1]:
		velosity[i, 0] *= damp
		location[i, 0] = xlim[1]
	if location[i, 1] < ylim[0]:
		velosity[i, 1] *= damp
		location[i, 1] = ylim[0]
	elif location[i, 1] > ylim[1]:
		velosity[i, 1] *= damp
		location[i, 1] = ylim[1]

def Pressure():
	global pressure, hashmap, rho, rho_0, h, k, n
	for i in range(n):
		rho[i] = 0
		for j in hashmap.neighbours(location[i]):
			r = location[j] - location[i]
			r2 = r.dot(r)
			if r2 <= h ** 2:
				rho[i] += poly6() * (h ** 2 - r2) ** 3
		pressure[i] = k * (rho[i] - rho_0)

def Force():
	global n, hashmap, location, h, visc, mass, fg, force, pressure, velosity
	kspiky = spiky()
	klaplacian = laplacian()
	for i in range(n):
		fp = np.array([0., 0.])
		fv = np.array([0., 0.])
		for j in hashmap.neighbours(location[i]):
			if i == j:
				continue
			r = location[j] - location[i]
			rl = math.sqrt(r.dot(r))
			if 0 < rl < h:
				fp += -r / rl * (pressure[i] + pressure[j]) / 2 / rho[j] * kspiky * ((h - rl) ** 2)
				fv += visc * mass * (velosity[j] - velosity[i]) / rho[j] * klaplacian * (h - rl)
		force[i] = fp + fv + fg

def euler():
	global dt, velosity, force, rho, location, hashmap, n
	for i in range(n):
		velosity[i] += dt * force[i] / rho[i]
		location[i] += dt * velosity[i]
		bound(i)
		hashmap.move(i, location[i])

def animate(self):
	global scat, location
	Pressure()
	Force()
	euler()
	scat.set_offsets(location)

i = 0
#for y in range(ylim[0] + 50, ylim[1] - 50, h):
#	for x in range(xlim[0] + 100, xlim[1] - 100, h):
for y in range(ylim[0], ylim[0] + 200, h):
    for x in range(xlim[0], xlim[1], h):
        if i < n:
            location[i] = x + np.random.normal(0, 0.25), y + np.random.normal(0, 0.25)
            #if y == 150#(ylim[1] - 100 - ylim[0]) / 2:
            velosity[i] = [100, 0]
            hashmap.move(i, location[i])
        i += 1
    pass

figure = plt.figure()
ax = plt.axes(xlim = (0, diameter), ylim = (0, height))
scat = ax.scatter([], [])
anim = animation.FuncAnimation(figure, animate, frames = range(5000), interval = 1000 / 60, blit = False, repeat = False, fargs = ())
Writer = animation.writers['ffmpeg']
writer = Writer(fps = 15, metadata = dict(artist = 'Me'), bitrate = 1800)
anim.save('finaltemp2.mp4')
plt.show()