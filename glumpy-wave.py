#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy, glumpy
import scipy.sparse as sparse
from scipy.misc import imread
from difference import Scheme

# Set step size
dx = 0.1
dy = 0.1
# Grid size
Nx = 800
Ny = 200
c = 340
dt = 1.00/100000

scheme = Scheme(dx,dy,dt,c,Nx,Ny,'cave.png')
scheme.method = scheme.cross4
scheme.setup()

comp = scheme.A
mask = scheme.mask

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

Z = unow.reshape(Ny,Nx).astype(numpy.float32)

fig = glumpy.figure( (2*Nx,2*Ny) )
image = glumpy.image.Image(Z, interpolation='nearest',
                           colormap=glumpy.colormap.Grey)

@fig.event
def on_draw():
    fig.clear()
    image.draw(0,0,0,fig.width,fig.height)

i = 0
printStep = 4

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
