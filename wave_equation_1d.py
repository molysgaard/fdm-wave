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

def main():
    # Grid size (1D)
    N = 1000.0

    dx = 0.1
    dt = 0.01

    c = 0.2

    print "Speed of scheme (dx/dt): %s" % (dx/dt)
    print "c: %s" % c

    s = c**2*(dt/dx)**2

    print "Computed s: %s" % s

    # Number of Time steps
    NumberOfTimeSteps = 10000

    # Create A matrix for solving
    diag = np.ones(N)*(-2)
    ldiag = np.ones(N)
    A = sparse.dia_matrix(([diag,ldiag,ldiag],[0,1,-1]),shape=(N,N))
    I = sparse.identity(N)

    # Function to start wave
    def start_wave(vector):
        length = len(vector)
        for x in range(length):
            vector[x] = np.sin(2*np.pi*x/length)
        return vector

    # Start wave
    U_2 = start_wave(np.zeros(N))

    # Taylor approximation for correct start
    U_1 = 0.5*(2*I+s*A)*U_2
    U_1[0] = 0
    U_1[-0] = 0
    U_exact_start = U_2

    # Create dir
    os.system("mkdir /tmp/1234")

    # Calculations
    PlotSteps = 2000

    # Error
    error = np.zeros(NumberOfTimeSteps)
    for i in range(2,NumberOfTimeSteps):

        U_new = (2*I + s*A)*U_1 - U_2
        U_exact = U_exact_start*np.cos(2*c*np.pi*i*(dt/dx)/N)

        # Add error to errorvector
        error[i] = max(abs(U_new - U_exact))

        if i%PlotSteps==0:
            print i
        #    fname = '/tmp/1234/tmp_%05d.png' %(i/PlotSteps)
        #    fig = plt.figure()
        #    ax = fig.add_subplot(111)
        #    ax.set_ylim(-1.5,1.5)
        #    plt.plot(U_new)
        #    plt.plot(U_exact)
        #    plt.savefig(fname,dpi=80)
        #    print "Pic %s of %s created" % (i/PlotSteps,NumberOfTimeSteps/PlotSteps)
        #    plt.clf()
        #    plt.close()

        # Update vecotrs
        U_2 = U_1
        U_1 = U_new
        U_1[0] = 0
        U_1[-0] = 0

    # Plot error vs time
    plt.plot(error)
    plt.title("Error plot \n N:%s,c:%s,dt:%s,dx:%s" % (N,c,dt,dx))
    plt.xlabel("Time")
    plt.ylabel("Max error")
    plt.savefig("N%s_c%s_dx%s_t%s.png" % (N,c,dx,NumberOfTimeSteps))

    ### Generate video of pngs and clear tmp_folders
    #os.system("ffmpeg -y -r 10 -sameq -i /tmp/1234/tmp_%05d.png movie.mp4")
    #os.system("mv movie.mp4 N%s_c%s_dx%s_t%s.mp4" % (N,c,dx,NumberOfTimeSteps))
    ##os.system("rm /tmp/1234/tmp*.png")
    ##os.system("rmdir /tmp/1234")
    #os.system("vlc N%s_c%s_dx%s_t%s.mp4" % (N,c,dx,NumberOfTimeSteps))

    print "Finished!"

main()
