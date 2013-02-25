#! /usr/bin/python

import scipy as sc
import scipy.sparse as sparse
import scipy.sparse.linalg
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

from mayavi import mlab
from tvtk.tools import visual

import os,sys

# Turn off gui

mlab.options.offscreen = True

# Grid size
N = 1001
# Set step size
h = 0.1

# Calculate dt
dt = h**2
k = dt/h
NumOfTimeSteps = int(500)

# Create mesh for plotting
X, Y = np.mgrid[:N,:N]

# Create tmp folder for pngs
os.system("mkdir /tmp/1234")

# Create A matrix

diag = np.zeros(N*N) - 4
ldiag = np.ones(N*N)

# Shcheme matrix for central differences
A = sparse.dia_matrix(([diag,ldiag,ldiag,ldiag,ldiag],[0,1,-1,-N,N]),shape=(N*N,N*N))
I = sparse.identity(N*N)

#print A.todense()

# Init U vectors
U_2 = np.zeros(N*N)
U_1 = np.zeros(N*N)

# Initial wave
U_1[(N*N)/2] = 0.2

fig = mlab.figure(size=(1000,1000))
camera = fig.scene.camera
#fig.scene.disable_render = True
#fig.scene.off_screen_rendering = True
mlab.view(focalpoint=(0,0,0),elevation=70,distance=70)
#camera.yaw(120)

#visual.set_viewer(fig)

U_new = U_1
s = mlab.surf(X,Y,U_new.reshape(N,N),extent=[0,100,0,100,0,50],vmax=0.2,vmin=-0.2)
## Run forloops to calculate
for n in range(1,NumOfTimeSteps+1):

    # Calculate the new vector
    U_new = 2*U_1 + k*A*U_1 - U_2
    # Set rand_values to zero
    # First N values
    U_new[:N] = 0
    # Last N values
    U_new[N*N-N:] = 0
    # Every Nth value
    U_new[::N] = 0
    # Every Nth +1 value
    U_new[N-1::N] = 0

    s.mlab_source.scalars = U_new.reshape(N,N)
    #fig.scene.render()
    # Render picture
    fname = '/tmp/1234/tmp_%05d.png' %n
    mlab.savefig(fname,size=(1000,1000))
    print "Pic %s of %s created" % (n,NumOfTimeSteps)

    # Update old vectors
    U_2 = U_1
    U_1 = U_new

mlab.close()
# Generate video of pngs and clear tmp_folders
os.system("ffmpeg -y -r 20 -sameq -i /tmp/1234/tmp_%05d.png movie.mp4")
os.system("rm /tmp/1234/tmp*.png")
os.system("rmdir /tmp/1234")

print "Finished!"
