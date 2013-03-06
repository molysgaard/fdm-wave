#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
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

class Scheme(object):
    def __init__(self, dx, dy, dt, c, Nx, Ny, mask_path):
        self.dx = dx
        self.dy = dy
        self.dt = dt
        self.c = c
        self.Nx = Nx
        self.Ny = Ny
        self.generate_methods()
        self.mask = self.mask_from_image(mask_path)

    def setup(self):
        self.A = self.A_from_method(self.method)
        self.mask = self.mask.reshape(Nx*Ny,1)

    def mask_from_image(self, path):
        mask = imread('cave.png')
        return 1-mask

    def A_from_method(self, method):
        Nx = self.Nx
        Ny = self.Ny

        def two_to_one(rows,cols,i,j):
            return i*cols+j

        def one_to_two(rows,cols,n):
            i = int(numpy.floor(float(n)/float(cols)))
            j = n-cols*i
            return (i,j)

        A = sparse.lil_matrix((Ny*Nx,Ny*Nx))
        for i in xrange(Ny*Nx):
            (x,y) = one_to_two(Ny,Nx,i)
            if self.mask[x,y]:
                for (ddx,ddy,coeff) in method:
                    j = two_to_one(Ny,Nx,x+ddx,y+ddy)
                    if i<Ny*Nx and j<Ny*Nx:
                        A[i,j] += coeff

        return sparse.dia_matrix(A)

    def generate_methods(self):
        dx = self.dx
        dy = self.dy
        dt = self.dt
        c = self.c

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

        self.cross2 = cross+[(0,0,2-4*c*c*dt/(dx*dy))]
        self.diag2 = diag+[(0,0,2-4*c*c*dt/(dx*dx+dy*dy))]
        self.crossdiag2 = mult(cross+diag+[(0,0,2-4*c*c*dt/(dx*dy))]+[(0,0,2-4*c*c*dt/(dx*dx+dy*dy))], 0.5)

        p = dt/(dx*dy)*c
        one = mult([(1,0,1.0),(-1,0,1.0),(0,1,1.0),(0,-1,1.0)],16.0)
        two = mult([(2,0,1.0),(-2,0,1.0),(0,2,1.0),(0,-2,1.0)],-1.0)
        gress = mult(one + two,p*p/12.0)
        self.cross4 = gress + [(0,0,2-5.0*p*p)]
