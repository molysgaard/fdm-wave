#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def main():
    num_measures = 16000
    printStep = 10

    h = 0.2/2.0
    Nr = 300*2
    Nc = 100*2

    level_2 = 100*2
    level_3 = 220*2

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
            if float(r)+20.0*np.sin(float(c)/200.0*2.0*np.pi+np.pi/2.0) > level_2:
                if float(r)+40.0*np.sin(float(c)/200.0*np.pi+np.pi/2.0) > level_3:
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
        return A
        #A = A.todense()
        
        #A[0:Nr,:]=0
        #A[(Nr-1)*Nc:Nr*Nc,:]=0
        #A[0:Nc:Nr*Nc,:]=0
        #A[Nc:Nc:Nr*Nc,:]=0

        #return sparse.csc_matrix(A)

    A = create_scheme()

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

    B = sparse.identity(Nr*Nc)
    A, B = add_absorbing_paper_first(A, B)
    B = sparse.csr_matrix(B)

    def create_old_mask():
        m = np.ones((Nr,Nc))
        m[:,0] = 0
        m[:,Nc-1] = 0
        m[0,:] = 0
        m[Nr-1,:] = 0
        return m.reshape(Nr*Nc,1)

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

    gs = gridspec.GridSpec(2, 1, height_ratios=[1,4])
    counter = 0
    measurements = np.zeros(num_measures)
    while counter < num_measures:
        rhs = A*unow - np.multiply(uold,old_mask)

        unew = spsolve(B, rhs).reshape(Nr*Nc,1)

        uold = unow
        unow = unew

        measurements[counter] = unow[5*Nc+Nc/2]
        print counter
        if counter%printStep==0:
            plt.subplot(gs[0])
            plt.plot(measurements[0:counter])
            plt.xlim((0,num_measures))
            plt.ylim((-0.013,0.022))

            plt.subplot(gs[1])
            tmp = unow-unow.min()+0.02
            tmp = np.multiply(tmp,savemask)
            tmp = tmp.reshape(Nr,Nc)
            tmp = np.transpose(tmp)
            image = plt.imshow(tmp)
            image.interpolation = 'nearest'
            plt.savefig('shots/shot%05d.png' % (counter/printStep), bbox_inches=0, dpi=120)

        counter += 1

main()
