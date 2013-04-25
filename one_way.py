#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import glumpy
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve

# schema can be note_first, note_second, paper_first, paper_second
schema = "paper_first"

printStep = 1

h = 0.2
Nr = 100
Nc = 200
c = 1.0
k = 0.01

a = c*c
p = k/h


###################################################################################
# the matrix for our difference method
def create_scheme():
    d = np.ones(Nr*Nc)
    M = sparse.dia_matrix(([-4*d,d,d,d,d],[0,-1,1,-Nc,Nc]),shape=(Nr*Nc,Nr*Nc))
    I = sparse.identity(Nr*Nc)
    A = 2*I + c*c*k*k/(h*h)*M
    A = sparse.lil_matrix(A)
    
    z = np.zeros(Nr*Nc)
    for i in xrange(Nr):
        A[i*Nc,:] = z # column 1
        A[i*Nc-1,:] = z # column Nc
        A[0,:] = z # row 1
        #A[i,:] = z
        #A[i*N-1,:] = z
        #A[N-i*N,:] = z

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
        m[i,i] = 1/h - 1/k
        m[i,i+Nc] = -(1/h+1/k)

        n[i,i] = -(1/h+1/k)
        n[i,i+Nc] = 1/h - 1/k

    for i in xrange(Nr):
        # discretization from the paper
        # column 1
        m[i*Nc,i*Nc] = 1/h - 1/k
        m[i*Nc,i*Nc+1] = -(1/h+1/k)

        n[i*Nc,i*Nc] = -(1/h+1/k)
        n[i*Nc,i*Nc+1] = 1/h - 1/k
        # column Nc
        m[(i+1)*Nc-1,(i+1)*Nc-1] = -(1/h - 1/k)
        m[(i+1)*Nc-1,(i+1)*Nc-2] = (1/h+1/k)

        n[(i+1)*Nc-1,(i+1)*Nc-1] = (1/h+1/k)
        n[(i+1)*Nc-1,(i+1)*Nc-2] = -(1/h - 1/k)

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
    for i in xrange(1,Nr-1):
        m[i,0] = 0 # column 1
        m[i,Nc-1] = 0 # column Nc
    for i in xrange(1,Nc-1):
        m[0,i] = 0 # row 1
    return m.reshape(Nr*Nc,1)


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

unow = create_gauss_wave(np.zeros((Nr,Nc)), Nr/2, Nc/2).reshape(Nr*Nc,1)
unow = np.matrix(unow)
print unow.shape
uold = unow
print A.shape, B.shape

####################################################################################
# beginning of rendering code

Z = unow.reshape(Nr,Nc).astype(np.float32)

fig = glumpy.figure( (2*Nc,2*Nr) )
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


    rhs = A*unow - np.multiply(uold,old_mask)

    if schema=="note_first":
        unew = rhs
    elif schema=="paper_first":
        unew = spsolve(B, rhs).reshape(Nr*Nc,1)

    unew = np.multiply(unew,mask)

    uold = unow
    unow = unew

    if i==printStep:
        i=0
        Z = (3000*unew+100).reshape(Nr,Nc).astype(np.uint8)
        image = glumpy.image.Image(Z, interpolation='nearest',
                    colormap=glumpy.colormap.LightBlue)

        image.update()
        fig.redraw()
        print np.sum(np.abs(Z))
    else:
        i+=1

glumpy.show()
