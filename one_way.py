#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import glumpy
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from scipy.misc import imsave

# schema can be note_first, note_second, paper_first, paper_second
schema = "paper_first"

num_measures = 16000
printStep = 100

snapshot_times = [420,2100,3480,4630,9400,10300,12100]

h = 0.2
Nr = 300
Nc = 100

#Nr = 100
#Nc = 50

initial_row = 25
initial_col = Nc/2

k = 0.01
p = k/h

cs = np.zeros((Nr,Nc))
for r in xrange(Nr):
    for c in xrange(Nc):
        #if float(r)+10.0*np.sin(float(c)/100.0*2.0*np.pi+np.pi/2.0) > 40:
        #    if float(r)+10.0*np.sin(float(c)/100.0*np.pi+np.pi/2.0) > 80:
        #        cs[r,c] = 3.3333333
        #    else:
        #        cs[r,c] = 6.6666666
        #else:
        #    cs[r,c] = 10.0
        if float(r)+10.0*np.sin(float(c)/100.0*2.0*np.pi+np.pi/2.0) > 100:
            if float(r)+20.0*np.sin(float(c)/100.0*np.pi+np.pi/2.0) > 220:
                cs[r,c] = 0.33333333
            else:
                cs[r,c] = 0.66666666
        else:
            cs[r,c] = 1.0

savemask = np.ones((Nr,Nc))
for r in xrange(Nr):
    for c in xrange(Nc):
        if cs[r,c] < 0.90 and c%4==0:
            savemask[r,c] = 0.0
        elif cs[r,c] < 0.60 and (c%4==0 or c%4==1):
            savemask[r,c] = 0.0
savemask = savemask.reshape(Nr*Nc,1)

cs = np.squeeze(np.asarray(cs.reshape(Nr*Nc,1)))

###################################################################################
# the matrix for our difference method
def create_scheme():
    d = cs**2*k*k/(h*h)

    M = sparse.dia_matrix(([-4*d,d,d,d,d],[0,-1,1,-Nc,Nc]),shape=(Nr*Nc,Nr*Nc))
    I = sparse.identity(Nr*Nc)
    A = 2*I + M
    A = sparse.lil_matrix(A)
    
    z = np.zeros(Nr*Nc)
    for i in xrange(Nr):
        A[i*Nc,:] = z # column 1
        A[i*Nc-1,:] = z # column Nc
        A[0,:] = z # row 1

    return sparse.csc_matrix(A)

A = create_scheme()

def add_absorbing_simple(m):
    for i in xrange(Nr):
        # simplest discretization from note
        m[i*Nc,i*Nc] = 1-a*p
        m[i*Nc,i*Nc+1] = a*p

    return m

def add_absorbing_paper_first(m,n):
    for i in xrange(Nc):
        # discretization from the paper
        # row 1
        c = cs[i]
        m[i,i] = (c/h - 1/k)
        m[i,i+Nc] = -(c/h+1/k)

        n[i,i] = -(c/h+1/k)
        n[i,i+Nc] = (c/h - 1/k)

        # row Nr
        c = cs[Nr*Nc-1]
        j = (Nr-1)*Nc+i
        m[j,j] = -(c/h - 1/k)
        m[j,j-Nc] = (c/h+1/k)

        n[j,j] = (c/h+1/k)
        n[j,j-Nc] = -(c/h - 1/k)

    for i in xrange(Nr):
        # discretization from the paper
        # column 1
        c = cs[i*Nc]
        m[i*Nc,i*Nc] = (c/h - 1/k)
        m[i*Nc,i*Nc+1] = -(c/h + 1/k)

        n[i*Nc,i*Nc] = -(c/h + 1/k)
        n[i*Nc,i*Nc+1] = (c/h - 1/k)
        # column Nc
        c = cs[(i+1)*Nc-1]
        m[(i+1)*Nc-1,(i+1)*Nc-1] = -(c/h - 1/k)
        m[(i+1)*Nc-1,(i+1)*Nc-2] = (c/h+1/k)

        n[(i+1)*Nc-1,(i+1)*Nc-1] = (c/h+1/k)
        n[(i+1)*Nc-1,(i+1)*Nc-2] = -(c/h - 1/k)

    return m, n

def add_absorbing_paper_second(m,n):
    for i in xrange(Nr):
        # discretization from the paper
        m[i*Nc,i*Nc] = 1/h - 1/k
        m[i*Nc,i*Nc+1] = -(1/h+1/k)

        n[i*Nc,i*Nc] = -(1/h+1/k)
        n[i*Nc,i*Nc+1] = 1/h - 1/k

    return m, n


if schema=="note_first":
    A = add_absorbing_simple(A)
elif schema=="paper_first":
    B = sparse.identity(Nr*Nc)
    A, B = add_absorbing_paper_first(A, B)
    B = sparse.csr_matrix(B)
elif schema=="paper_second":
    B = sparse.identity(Nr*Nc)
    A, B = add_absorbing_paper_second(A, B)
    B = sparse.csr_matrix(B)

def create_mask():
    m = np.ones((Nr,Nc))
    m[Nr-1,:] = 0 # row Nr
    return m.reshape(Nr*Nc,1)

def create_old_mask():
    m = np.ones((Nr,Nc))
    m[:,0] = 0
    m[:,Nc-1] = 0
    m[0,:] = 0
    m[Nr-1,:] = 0
    return m.reshape(Nr*Nc,1)

mask = create_mask()
old_mask = create_old_mask()



###################################################################################
# code to create the initial state
def create_gauss_wave(initial, x_0, y_0):
    def gauss(x,y):
        return 1.0/(2.0*np.pi)*np.exp(-(x**2 + y**2)/2.0)

    for i in xrange(-20,20):
        for j in xrange(-20,20):
            initial[x_0+j,y_0+i]=gauss(float(j)/4.0,float(i)/4.0)
    return initial

uold = create_gauss_wave(np.zeros((Nr,Nc)), initial_row, initial_col).reshape(Nr*Nc,1)
uold = np.matrix(uold)
unow = spsolve(B,0.5*A*uold).reshape(Nr*Nc,1)

####################################################################################
# beginning of rendering code

Z = unow.reshape(Nr,Nc).astype(np.float32)

fig = glumpy.figure( (2*Nc,2*Nr) )
image = glumpy.image.Image(Z, interpolation='nearest',
                           colormap=glumpy.colormap.Grey)

@fig.event
def on_draw():
    pass
    #fig.clear()
    #image.draw(0,0,0,fig.width,fig.height)

counter = 0
measurements = np.zeros(num_measures)
@fig.event
def on_idle(dt):
    global uold, unow, image, counter

    rhs = A*unow - np.multiply(uold,old_mask)

    if schema=="note_first":
        unew = rhs
    elif schema=="paper_first":
        unew = spsolve(B, rhs).reshape(Nr*Nc,1)

    #unew = np.multiply(unew,mask)

    uold = unow
    unow = unew

    if counter == len(measurements):
        plt.plot(measurements)
        plt.show()
        exit()

    measurements[counter] = unow[2*Nc+Nc/2]
    print counter
    if counter in snapshot_times or counter%printStep==0:
        tmp = unow-unow.min()
        imsave('shots/shot%04d.png' % counter, np.multiply(tmp,savemask).reshape(Nr,Nc))


    #if counter%printStep==0:
    #    Z = (3000*unew+100).reshape(Nr,Nc).astype(np.uint8)
    #    image = glumpy.image.Image(Z, interpolation='nearest',
    #                colormap=glumpy.colormap.LightBlue)

    #    image.update()
    #    fig.redraw()
    #    #print np.sum(np.abs(Z))
    counter += 1

glumpy.show()
