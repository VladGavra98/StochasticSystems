# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 17:54:22 2021

Symmetrical aircraft responses for a rigid aircraft in symmetrical atmospheric turbulence conditions.

Cessna  Ce500 Citation I - configuration start (2)

@author: vladg
"""
import numpy as np
import control.matlab as cm
import control as c
np.set_printoptions(precision=5)

class Cessna:

    def __init__(self,V=59.9,muc=113,sigma_u=2,sigma_w=3,Lg=150):
        """
        Cessna contrusctor

        Input:   V, muc, sigma_u, sigma_w, Lg

        """

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                   # Prelimianry checks.
          # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if muc < 90 or muc > 130:
            print("\n\n   Error in the input parameters ! \n\n")

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                          # AIRCRAFT FLIGHT CONDITION 'Start(2)'.
          # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        self.V     = V
        g0    = 9.80665
        m     = 6035
        self.muc   = muc
        twmuc = 2*muc
        KY2   = 0.893
        c     = 2.022
        S     = 24.2
        lh    = 5.5

         # TURBULENCE PARAMETERS

        # sigma_u = 2
        # sigma_w = 3
        # Lg      = 150
        self.sigma_u = sigma_u
        self.sigma_w = sigma_w

        self.sigma_ug   = sigma_u/V;
        self.sigma_ag   = sigma_w/V;

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                              # AIRCRAFT SYMMETRIC AERODYNAMIC DERIVATIVES
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CX0  = 0.0000
        CZ0  =-1.250
        Cm0  =  0.0000

        CXu  =-0.2510
        CZu  =-2.5000
        Cmu  = 0.0000

        CXa  = 0.5120
        CZa  = -5.160
        Cma  = -0.430

        CXq  = 0.0000
        CZq  =-3.8600
        Cmq  = -7.040

        CXd  = 0.0000
        CZd  =-0.6238
        Cmd  = -1.553

        CXfa = 0.0000
        CZfa =-1.4700
        Cmfa = -3.750

        CZfug = 0.000
        Cmfug = -Cm0*lh/c;

        CZfag= CZfa-CZq;
        Cmfag=  Cmfa-Cmq;


         #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         #                 CALCULATION OF AIRCRAFT SYMMETRIC STABILITY DERIVATIVES
         #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        xu   = (V/c)*(CXu/twmuc);
        xa   = (V/c)*(CXa/twmuc);
        xt   = (V/c)*(CZ0/twmuc);
        xq   = 0;
        xd   = (V/c)*(CXd/twmuc);
        xug  = xu;
        xfug = 0;
        xag  = xa;
        xfag = 0;

        zu   = (V/c)*( CZu/(twmuc-CZfa));
        za   = (V/c)*( CZa/(twmuc-CZfa));
        zt   = (V/c)*(-CX0/(twmuc-CZfa));
        zq   = (V/c)*((CZq+twmuc)/(twmuc-CZfa));
        zd   = (V/c)*( CZd/(twmuc-CZfa));
        zug  = zu;
        zfug = (V/c)*( CZfug/(twmuc-CZfa));
        zag  = za;
        zfag = (V/c)*( CZfag/(twmuc-CZfa));

        mu   = (V/c)*(( Cmu+CZu*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
        ma   = (V/c)*(( Cma+CZa*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
        mt   = (V/c)*((-CX0*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
        mq   = (V/c)*(Cmq+Cmfa*(twmuc+CZq)/(twmuc-CZfa))/(twmuc*KY2);
        md   = (V/c)*((Cmd+CZd*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
        mug  = mu;
        mfug = (V/c)*(Cmfug+CZfug*Cmfa/(twmuc-CZfa))/(twmuc*KY2);
        mag  = ma;
        mfag = (V/c)*(Cmfag+CZfag*Cmfa/(twmuc-CZfa))/(twmuc*KY2);


        # Matrix Dimensions

        n  = 4        # 4 orignal a/c states
        m  = 1      # 1 input (elevator deflection)
        p  = (5,2)  # 5 outputs

        # #     Creating the different c-matrices (c1, c2 &c3) for symmetrical flight
        # #c1 matrix
        # c1      = np.zeros((n,n))
        # c1[0,0] = -2*muc*(c/V)
        # c1[1,1] = (CZfa - 2*muc)*(c/V)
        # c1[2,2] = -(c/V)
        # c1[3,1] = Cmfa*(c/V)
        # c1[3,3] = -2*muc*KY2*((c/V)**2)

        # #c2 matrix
        # c2      = np.zeros((n,n))
        # c2[0,0] = -CXu
        # c2[0,1] = -CXa
        # c2[0,2] = -CZ0
        # c2[0,3] = -CXq*(c/V)
        # c2[1,0] = -CZu
        # c2[1,1] = -CZa
        # c2[1,2] = -CX0
        # c2[1,3] = -(CZq + 2*muc)*(c/V)
        # c2[2,3] = -(c/V)
        # c2[3,0] = -Cmu
        # c2[3,1] = -Cma
        # c2[3,3] = -Cmq*(c/V)

        # #c3 matrix
        # c3 = np.zeros((n,m))
        # c3[0,0] = -CXd
        # c3[1,0] = -CZd
        # c3[3,0] = -Cmd

         #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         #                         Basic Aircraft  Symmetric Model
         #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        A_s = np.array([[xu, xa, xt, 0],
                        [zu, za, zt, zq],
                        [0, 0, 0, V/c],
                        [mu, ma, mt, mq]])

        B_s = np.array([[xd],
                        [zd],
                        [0],
                        [md]])



         #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         #                      Symmetric Response to Turbulence
         #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        A_12 = np.zeros((n,3))
        A_12 = np.array([[xug, xag, 0],
                         [zug-zfug*(c/Lg), zag, zfag*(c/V)],
                         [0,                0,            0],
                         [mug-mfug*(c/Lg), mag, mfag *(c/V)]])

        A_22 = np.zeros((3,3))

        A_22 = np.array([[-V/Lg,      0,          0],
                         [0,          0,          1],
                         [0, -V**2/(Lg**2), -2*(V/Lg)]])


         #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         #                Complete Aircraft Dynamics Model (symmetric + turbulence)
         #                         A = [A_11  A_12]
         #                             [0     A_22]
         #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        self.A  = np.zeros((n+3,n+3))
        self.B  = np.zeros((n+3,3))

        self.A[:n,:n] = A_s
        self.A[:n,n:] = A_12
        self.A[n:,n:] = A_22


        self.B[:n,:1] = B_s
        self.B[:,1:]  = np.array([[0,                                   0],
                            [zfug*(c/V)*self.sigma_ug*np.sqrt(2*V/Lg), zfag*(c/V)*self.sigma_ag*np.sqrt(3*V/Lg)],
                            [0,                                        0],
                            [mfug*(c/V)*self.sigma_ug*np.sqrt(2*V/Lg),  mfag*(c/V)*self.sigma_ag*np.sqrt(3*V/Lg)],
                            [self.sigma_ug*np.sqrt(2*V/Lg),                0],
                            [0,                                  self.sigma_ag*np.sqrt(3*V/Lg)],
                            [0,                           (1-2*np.sqrt(3))*self.sigma_ag*np.power(V/Lg,1.5)]])

        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #                         Aircraft State Space
        #                         az = V*q - V * alpha
        #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        p = 5
        self.C = np.array([[1, 0, 0, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0, 0, 0],
                           [0, 0, 1, 0, 0, 0, 0],
                           [0, 0, 0, 1, 0, 0, 0],
                           [-V*zu/g0, -V*za/g0, zt/g0, (V*V/c - V*zq)/g0, -V*(zug-zfug*(c/Lg))/g0, -V*zag/g0, -V*zfag*(c/V)/g0]])

        self.D      = np.zeros((p,3))
        self.D[1,:] = -V*np.array([zd, zfug*(c/V)*self.sigma_ug*np.sqrt(2*V/Lg), zfag*(c/V)*self.sigma_ag*np.sqrt(3*V/Lg)])/g0


        #Check the \g0 !!!

    def state_space(self):
        """
        Method to return the Cessna continuous-time State-Space model
        for the selected configuration paramters.

        If a different configuration is deisred -> initiliase the aircraft with new parameters

        Input : -

        Output: control.matlab.StatePsace Object (LTI system)
                x  = [u ,alpha ,theta, qc/V,  u_g,  alpha_g,  alpha*_g]
                u  = [delta_e , w1 , w3]
                y  = [u ,alpha ,theta, qc/V,  n_z ]


        """
        # prelimianry checks

        if(np.shape(self.A) != (7,7)):
            print("Wrong matrix A\n")

        if(np.shape(self.B) != (7,3)):
            print("Wrong matrix B\n")

        if(np.shape(self.C) != (5,7)):
            print("Wrong matrix C\n")

        if(np.shape(self.D) != (5,3)):
            print("Wrong matrix D\n")

        system = cm.minreal(cm.ss(self.A,self.B,self.C,self.D))

        return system





































