#! /usr/bin/python

import scipy as sc
import scipy.sparse as sparse
import scipy.sparse.linalg
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib import cm

import os,sys
import math

# Grid size (1D)
N = 1000.0

dx = 0.1
dt = dx**2

c = 0.5

print "Speed of scheme (dx/dt): %s" % (dx/dt)
print "c: %s" % c

s = c**2*(dt/dx)**2

print "Computed s: %s" % s

# Number of Time steps
NumberOfTimeSteps = 6400

# Create A matrix for solving
diag = np.ones(N)*(-2)
ldiag = np.ones(N)
A = sparse.dia_matrix(([diag,ldiag,ldiag],[0,1,-1]),shape=(N,N))
I = sparse.identity(N)

# Function to start wave
def start_wave(vector):
    length = len(vector)
    for x in range(length):
        if x > 990:
            print x
        #vector[x] = np.exp(-(x-length/2.0)**2/sigma)
        vector[x] = 0.3*np.sin(20*np.pi*x/length)
    return vector

# Start wave
U_1 = start_wave(np.zeros(N))
U_2 = U_1

# Create dir
os.system("mkdir /tmp/1234")

# Calculations
PlotSteps = 30
for i in range(NumberOfTimeSteps):

    U_new = (2*I + s*A)*U_1 - U_2

    if i%PlotSteps==0:
        fname = '/tmp/1234/tmp_%05d.png' %(i/PlotSteps)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_ylim(-1,1)
        plt.plot(U_new)
        plt.savefig(fname,dpi=300)
        print "Pic %s of %s created" % (i/PlotSteps,NumberOfTimeSteps/PlotSteps)
        plt.clf()
        plt.close()

    # Update vecotrs
    U_2 = U_1
    U_1 = U_new
    #U_1[N/2.0] = 0.3*np.sin(np.pi*i/NumberOfTimeSteps)
    U_1[0] = 0
    U_1[-0] = 0


# Generate video of pngs and clear tmp_folders
os.system("ffmpeg -y -r 20 -sameq -i /tmp/1234/tmp_%05d.png movie.mp4")
#os.system("rm /tmp/1234/tmp*.png")
#os.system("rmdir /tmp/1234")
os.system("vlc movie.mp4")

print "Finished!"
