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
import math

# Turn off gui
mlab.options.offscreen = True

# Grid size
N = 101
# Set step size
h = 0.01

# Calculate dt
dt = h**2
k = dt/h

# Number of steps
NumOfTimeSteps = int(500)

# Create mesh for plotting
X, Y = np.mgrid[:N,:N]

# Create tmp folder for pngs
os.system("mkdir /tmp/1234")

# Create A matrix

diag = np.zeros(N*N) - 4
ldiag = np.ones(N*N)

# Shcheme matrix for central differences
#A = sparse.dia_matrix(([diag,ldiag,ldiag,ldiag,ldiag],[0,1,-1,-N+1,N]),shape=(N*N,N*N))
I = sparse.identity(N*N)

#print A.todense()
method = [(0,0,-4),(1,0,1),(-1,0,1),(0,1,1),(0,-1,1)]

def two_to_one(rows,cols,i,j):
    return i*cols+j

def one_to_two(rows,cols,n):
    i = int(math.floor(float(n)/float(cols)))
    j = n-cols*i
    return (i,j)

mask = np.ones((N,N))
mask[0][:]=0
mask[N-1][:]=0
for i in xrange(30,40):
    for j in xrange(30,40):
        mask[i,j]=0
for i in xrange(N):
    mask[i][0]=0
    mask[i][N-1]=0

A = sparse.lil_matrix((N*N,N*N))
for i in xrange(N*N):
    (x,y) = one_to_two(N,N,i)
    if mask[x][y]:
        for (dx,dy,c) in method:
            j = two_to_one(N,N,x+dx,y+dy)
            if i<N*N and j<N*N:
                A[i,j] = c

A = sparse.dia_matrix(A)
comp = (2*I +k*A)
mask = mask.reshape(N*N,1)

# Initial waves
U_1 = np.zeros((N,N))

def gauss(x,y):
    return 1.0/(2.0*math.pi)*math.exp(-(x^2 + y^2)/2)

x_0 = N/2.0
y_0 = N/2.0

for i in xrange(-5,5):
    for j in xrange(-5,5):
        U_1[x_0+j,y_0+i]=gauss(j/5,i/5)

U_1 = U_1.reshape(N*N,1)
U_2 = U_1

fig = mlab.figure(size=(1000,1000))
camera = fig.scene.camera
#fig.scene.disable_render = True
#fig.scene.off_screen_rendering = True
mlab.view(focalpoint=(0,0,0),elevation=70,distance=70)
#camera.yaw(120)

#visual.set_viewer(fig)

U_new = U_1
s = mlab.surf(X,Y,U_new.reshape(N,N),extent=[0,100,0,100,0,50],vmax=0.2,vmin=-0.2)

# Run forloops to calculate
i = 100
for n in range(1,NumOfTimeSteps+1):

    # Calculate the new vector
    U_new = (2*I + (1.0/100)*k*A)*U_1 - U_2
    #U_new = comp*U_1 - U_2
    U_new = np.multiply(U_new,mask)

    if (n%i==0):
        s.mlab_source.scalars = U_new.reshape(N,N)
        #fig.scene.render()
        # Render picture
        fname = '/tmp/1234/tmp_%05d.png' % int(n/i)
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
os.system("vlc movie.mp4")

print "Finished!"
