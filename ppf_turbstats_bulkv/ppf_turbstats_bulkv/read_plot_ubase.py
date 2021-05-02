import matplotlib as mpl
from matplotlib import pylab
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from numpy import pi,cos,arange, ones
#############################
#############################
'''
read and Ubase, Umean, etc
Author: Pratik Aghor
'''
#############################
# filename = 'data/Ubase.asc'

# load file, skip 1st row meant for chflow
Ubase = np.loadtxt('Ubase.asc', skiprows=1)
Umean = np.loadtxt('Umean.asc', skiprows=1)
Umeany = np.loadtxt('Umeany.asc', skiprows=1)

# Wbase = np.loadtxt('data/Wbase.asc', skiprows=1)
Ny = len(Ubase)
a = -1; b = 1;
# y = np.linspace(b, a, Ny) # Cheb goes from 1 to -1
r = (b - a) / 2
c = (b + a) / 2
y = c + r*cos(pi*arange(0,Ny)/(Ny-1)) # y = Chebyshev-Gauss-Lobatto grid

U0 = np.ones(Ny) - y*y
U0y = -2.0*y
#############################
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()  # Get current axes

ax.plot(U0, y, '-b', label=r'$(1-y^{2})$')  # Plot analytical
ax.plot(Ubase, y, '--r', label="Ubase")  # Plot numerical

ax.set_xlabel(r'$U$', fontsize=20)  # Set x label
ax.set_ylabel(r'$y$', fontsize=20)  # Set y label
ax.legend(loc=1)
fig.savefig('Ubase_ppf.png')

#############################
fig = plt.figure(2)  # Create a figure instance
ax = fig.gca()  # Get current axes

ax.plot(U0, y, '-b', label=r'$(1-y^{2})$')  # Plot analytical
ax.plot(Umean, y, '--r', label="Umean")  # Plot numerical

ax.set_xlabel(r'$Umean$', fontsize=20)  # Set x label
ax.set_ylabel(r'$y$', fontsize=20)  # Set y label
ax.legend(loc=1)
fig.savefig('Umean_ppf.png')

#############################
fig = plt.figure(3)  # Create a figure instance
ax = fig.gca()  # Get current axes

ax.plot(U0y, y, '-b', label=r'$(-2y)$')  # Plot analytical
ax.plot(Umeany, y, '--r', label=r'$Umean_{y}$')  # Plot numerical

ax.set_xlabel(r'$Umean_{y}$', fontsize=20)  # Set x label
ax.set_ylabel(r'$y$', fontsize=20)  # Set y label
ax.legend(loc=1)
fig.savefig('Umeany_ppf.png')

#############################
