#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import glumpy
import scipy.sparse as sparse

printStep = 1

h = 0.2
N = 200
c = 1.0
k = 0.01

a = c*c
p = k/h


###################################################################################
# the matrix for our difference method
def create_scheme():
    d = np.ones(N*N)
    M = sparse.dia_matrix(([-4*d,d,d,d,d],[0,-1,1,-N,N]),shape=(N*N,N*N))
    I = sparse.identity(N*N)
    A = 2*I + c*c*k*k/(h*h)*M
    A = sparse.lil_matrix(A)
    
    z = np.zeros(N*N)
    for i in xrange(N):
        A[i*N,:] = z
        #A[i,:] = z
        #A[i*N-1,:] = z
        #A[N-i*N,:] = z

    return sparse.csc_matrix(A)

A = create_scheme()


def add_absorbing_simple(m):
    for i in xrange(N):
        # simplest discretization from note
        m[i*N,i*N] = 1-a*p
        m[i*N,i*N+1] = a*p

    return m

def add_absorbing_paper(m,n):
    for i in xrange(N):
        # discretization from the paper
        m[i*N,i*N] = 1/h - 1/k
        m[i*N,i*N+1] = -(1/h+1/k)

        n[i*N,i*N] = -(1/h+1/k)
        n[i*N,i*N+1] = 1/h - 1/k

    return m, n


#A = add_absorbing_simple(A)

B = sparse.identity(N*N)
A, B = add_absorbing_paper(A, B)

def create_mask():
    m = np.ones((N,N))
    m[0,:] = 0
    m[N-1,:] = 0
    #m[:,0] = 0 # comment out for absorbing boundaries
    m[:,N-1] = 0
    return m.reshape(N*N,1)

def create_old_mask():
    m = np.ones((N,N))
    for i in xrange(1,N-1):
        m[i,0] = 0
    return m.reshape(N*N,1)


mask = create_mask()
old_mask = create_old_mask()



###################################################################################
# code to create the initial state
def create_gauss_wave(initial, x_0, y_0):
    def gauss(x,y):
        return 1/(2*np.pi)*np.exp(-(x**2 + y**2)/2)

    for i in xrange(-20,20):
        for j in xrange(-20,20):
            initial[x_0+j,y_0+i]=gauss(j/4,i/4)
    return initial

unow = create_gauss_wave(np.zeros((N,N)), N/2, N/2).reshape(N*N,1)
uold = unow

####################################################################################
# beginning of rendering code

Z = unow.reshape(N,N).astype(np.float32)

fig = glumpy.figure( (2*N,2*N) )
image = glumpy.image.Image(Z, interpolation='nearest',
                           colormap=glumpy.colormap.Grey)

@fig.event
def on_draw():
    fig.clear()
    image.draw(0,0,0,fig.width,fig.height)

i = 0
@fig.event
def on_idle(dt):
    global uold, unow, image, i

    #unew = A*unow - np.multiply(uold,old_mask)
    unew = sparse.linalg.spsolve(B, A*unow - np.multiply(uold,old_mask))
    unew = np.multiply(unew,mask)

    uold = unow
    unow = unew

    if i==printStep:
        i=0
        Z = unew.reshape(N,N).astype(np.float32)
        Z = (3000*Z+100).astype(np.uint8)
        image = glumpy.image.Image(Z, interpolation='nearest',
                    colormap=glumpy.colormap.LightBlue)

        image.update()
        fig.redraw()
        print np.sum(np.abs(Z))
    else:
        i+=1

glumpy.show()
