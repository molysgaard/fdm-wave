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


# Set step size
dx = 0.1
dy = 0.1
# Grid size
N = 100 # 1x1 metres room
# wave speed, must be less than one to make sense
c = 340 # sound speed 340 m/s
dt = 1.0/100000000 # we check the pressure each 1/1000 second

NumOfTimeSteps = 400
plotSteps = 2

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

def mult(xs,a):
    ret = []
    for (x,y,z) in xs:
        ret.append((x,y,a*z))
    return ret

# second order compass method
cross = mult([(1,0,1),(-1,0,1)],dt/(dx*dx)) + mult([(0,1,1),(0,-1,1)],dt/(dy*dy))

diag = mult([(1,1,1),(-1,-1,1),(-1,1,1),(1,-1,1)],dt/(dx*dx+dy*dy))

#method = cross+[(0,0,-4*dt/(dx*dy))]
#method = diag+[(0,0,-4*dt/(dx*dx+dy*dy))]

method = mult(cross+diag+[(0,0,-4*dt/(dx*dy))]+[(0,0,-4*dt/(dx*dx+dy*dy))], 0.5)

#p = 1
#one = mult([(1,0,1.0),(-1,0,1.0),(0,1,1.0),(0,-1,1.0)],16.0)
#two = mult([(2,0,1.0),(-2,0,1.0),(0,2,1.0),(0,-2,1.0)],-1.0)
#gress = mult(one + two,p*p/12.0)
#method = gress + [(0,0,-5.0*p*p)]

print method

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
        for (ddx,ddy,coeff) in method:
            j = two_to_one(N,N,x+ddx,y+ddy)
            if i<N*N and j<N*N:
                A[i,j] += coeff

A = sparse.dia_matrix(A)
comp = 2*I + c*c*A
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

#os.system("rm /tmp/1234/tmp*.png")
#os.system("rmdir /tmp/1234")
os.system("vlc movie.mp4")

print "Finished!"
