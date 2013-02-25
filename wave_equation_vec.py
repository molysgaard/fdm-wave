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
import math

# Grid size
N = 400

# Set step size
h = 0.01
# wave speed, must be less than one to make sense
c = 0.005
dt = 0.001

k = dt/h
NumOfTimeSteps = int(10000)
plotSteps = 20

print c*dt*h*h

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

method = [(0,0,-4),(1,0,1),(-1,0,1),(0,1,1),(0,-1,1)]
#method = [(0,0,-4),(1,0,1),(-1,0,1),(0,1,1),(0,-1,1),(1,1,0),(1,-1,0),(-1,1,0),(-1,-1,0)]

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
comp = (2*I +k*c**2*A)
mask = mask.reshape(N*N,1)

# Init U vectors
U_1 = np.zeros((N,N))

# Initial wave
def gauss(x,y):
    return 1/(2*math.pi)*math.exp(-(x**2 + y**2)/2)

# the origin of the initial wave
x_0 = N/2
y_0 = N/2

for i in xrange(-20,20):
    for j in xrange(-20,20):
        U_1[x_0+j,y_0+i]=gauss(j/5.0,i/5.0)

U_1 = U_1.reshape(N*N,1)
U_2 = U_1

colors = [[(0, 0.9, 1, 1.) for x in range(N)] for y in range(N)]

## Run forloops to calculate
i = 0
for n in range(1,NumOfTimeSteps+1):

    # Calculate the new vector
    #U_new = (2*I + k*A)*U_1 - U_2
    U_new = comp*U_1 - U_2
    U_new = np.multiply(U_new,mask)

    if n%plotSteps==0:
        ## Render picture

        fname = '/tmp/1234/tmp_%05d.png' %i
        i+=1

        #plt.imshow(U_new.reshape(N,N))
        #plt.savefig(fname, dpi=100)
        #plt.clf()

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection='3d')
        surf = ax.plot_surface(X,Y,np.array(U_new.reshape(N,N)), linewidth=0, shade=True, antialiased=False, facecolors=colors)
        #surf = ax.plot_surface(X,Y,np.array(U_new.reshape(N,N)), rstride=5, cstride=5, cmap=cm.ocean)
        #ax.view_init(elev=30, azim=n/2)
        ax.set_zlim3d(-0.1,0.1)
        plt.savefig(fname, dpi=150)
        plt.clf()
        plt.close()
        print "Pic %s of %s created" % (n,NumOfTimeSteps)

    U_2 = U_1
    U_1 = U_new

## Generate video of pngs and clear tmp_folders
#os.system("ffmpeg -y -r 20 -sameq -i /tmp/1234/tmp_%05d.png movie.mp4")
os.system("ffmpeg -y -r 20 -sameq -i /tmp/1234/tmp_%05d.png movie.mp4")

os.system("rm /tmp/1234/tmp*.png")
os.system("rmdir /tmp/1234")
os.system("vlc movie.mp4")

print "Finished!"
