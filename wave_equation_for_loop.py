#! /usr/bin/python

import scipy as sc
import scipy.sparse as sparse
import scipy.sparse.linalg
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import os,sys

# Grid size
N = 100

# Set step size
dx = 0.1
dy = 0.1

# Calculate dt
dt = dx**2
NumOfTimeSteps = int(500)

# Initial conditions
u = np.zeros((N+1,N+1))

U_old2 = u

initial_u = np.copy(u)
initial_u[N/2,N/2] = 0.2
U_old1 = initial_u

# Create mesh for plotting
X, Y = np.mgrid[:N+1,:N+1]

# Create tmp folder for pngs
os.system("mkdir /tmp/1234")

# Run forloops to calculate
for n in range(1,NumOfTimeSteps+1):
    Unew = np.zeros((N+1,N+1))
    for m in range(1,N):
        for k in range(1,N):
            Unew[m,k] = 2*U_old1[m,k]-U_old2[m,k]+(dt/dx)*(U_old1[m+1,k]-2*U_old1[m,k]+U_old1[m-1,k])+(dt/dy)*(U_old1[m,k+1]-2*U_old1[m,k]+U_old1[m,k-1])

    # Render picture
    fname = '/tmp/1234/tmp_%05d.png' %n
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    surf = ax.plot_surface(X,Y,Unew, rstride=5, cstride=5, cmap=cm.ocean)
    ax.view_init(elev=30, azim=n/2)
    ax.set_zlim3d(-0.3,0.3)
    ax.set_axis_off()
    plt.savefig(fname, dpi=300)
    plt.clf()
    plt.close()
    print "Pic %s of %s created" % (n,NumOfTimeSteps)

    U_old2 = U_old1
    U_old1 = Unew

# Generate video of pngs and clear tmp_folders
os.system("ffmpeg -y -r 20 -b 1800 -i /tmp/1234/tmp_%05d.png movie.mp4")
os.system("rm /tmp/1234/tmp*.png")
os.system("rmdir /tmp/1234")

print "Finished!"
