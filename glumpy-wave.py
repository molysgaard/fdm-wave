#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy, glumpy
import scipy.sparse as sparse
from scipy.misc import imread

# Set step size
dx = 0.1
dy = 0.1
# Grid size
Nx = 800
Ny = 200
c = 340
dt = 1.00/100000

I = sparse.identity(Ny*Nx)

# Create A matrix
def mult(xs,a):
    ret = []
    for (x,y,z) in xs:
        ret.append((x,y,a*z))
    return ret

# second order cross method
cross = mult([(1,0,1),(-1,0,1)],dt/(dx*dx)) + mult([(0,1,1),(0,-1,1)],dt/(dy*dy))

# second order diagonal method
diag = mult([(1,1,1),(-1,-1,1),(-1,1,1),(1,-1,1)],dt/(dx*dx+dy*dy))

#method = cross+[(0,0,2-4*c*c*dt/(dx*dy))]
#method = diag+[(0,0,2-4*c*c*dt/(dx*dx+dy*dy))]
#method = mult(cross+diag+[(0,0,2-4*c*c*dt/(dx*dy))]+[(0,0,2-4*c*c*dt/(dx*dx+dy*dy))], 0.5)

p = dt/(dx*dy)*c
one = mult([(1,0,1.0),(-1,0,1.0),(0,1,1.0),(0,-1,1.0)],16.0)
two = mult([(2,0,1.0),(-2,0,1.0),(0,2,1.0),(0,-2,1.0)],-1.0)
gress = mult(one + two,p*p/12.0)
method = gress + [(0,0,2-5.0*p*p)]

print method

def two_to_one(rows,cols,i,j):
    return i*cols+j

def one_to_two(rows,cols,n):
    i = int(numpy.floor(float(n)/float(cols)))
    j = n-cols*i
    return (i,j)

#mask = numpy.ones((N,N))
#mask[0][:]=0
#mask[N-1][:]=0
#for i in xrange(30,40):
#    for j in xrange(30,40):
#        mask[i,j]=0
#for i in xrange(N):
#    mask[i][0]=0
#    mask[i][N-1]=0

mask = imread('cave.png')
mask = 1-mask

A = sparse.lil_matrix((Ny*Nx,Ny*Nx))
for i in xrange(Ny*Nx):
    (x,y) = one_to_two(Ny,Nx,i)
    if mask[x,y]:
        for (ddx,ddy,coeff) in method:
            j = two_to_one(Ny,Nx,x+ddx,y+ddy)
            if i<Ny*Nx and j<Ny*Nx:
                A[i,j] += coeff

A = sparse.dia_matrix(A)
#comp = 2*I + c*c*A
comp = A
mask = mask.reshape(Ny*Nx,1)

# Init U vectors
unow = numpy.zeros((Ny,Nx))

# Initial wave
def gauss(x,y):
    return 1/(2*numpy.pi)*numpy.exp(-(x**2 + y**2)/2)

# the origin of the initial wave
x_0 = 95
y_0 = 74

for i in xrange(-20,20):
    for j in xrange(-20,20):
        unow[x_0+j,y_0+i]=gauss(j/4.0,i/4.0)

unow = unow.reshape(Ny*Nx,1)
uold = unow
############################

plotSteps = 2

Z = unow.reshape(Ny,Nx).astype(numpy.float32)

fig = glumpy.figure( (2*Nx,2*Ny) )
image = glumpy.image.Image(Z, interpolation='nearest',
                           colormap=glumpy.colormap.Grey)

#image_mask = glumpy.image.Image(mask.astype(numpy.uint8), interpolation='nearest',
#        colormap=glumpy.colormap.Grey)

@fig.event
def on_draw():
    fig.clear()
    image.draw(0,0,0,fig.width,fig.height)
    #image_mask.draw(0,0,0,fig.width,fig.height)

i = 0
printStep = 10

@fig.event
def on_idle(dt):
    global uold, unow, image, i

    unew = comp*unow - uold
    unew = numpy.multiply(unew,mask)

    uold = unow
    unow = unew

    if i==printStep:
        i=0
        Z = unew.reshape(Ny,Nx).astype(numpy.float32)
        Z = (3000*Z+100).astype(numpy.uint8)
        image = glumpy.image.Image(Z, interpolation='nearest',
                    colormap=glumpy.colormap.LightBlue)

        image.update()
        fig.redraw()
        print numpy.sum(numpy.abs(Z))
    else:
        i+=1

glumpy.show()
