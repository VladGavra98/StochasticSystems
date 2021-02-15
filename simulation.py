# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 20:16:15 2021


Simulation of the Cessna Landing model for stochastic inputs w1(t) and w3(t)

@author: vladg
"""
import numpy as np
import scipy as sp
import scipy.signal
import control.matlab as cm
import control
import  Cessna_model
import numpy.random
import matplotlib.pyplot as plt

np.set_printoptions(precision=3)
# np.random.seed(4)

  #  Data and constants:
V     = 59.9
m     = 6035
muc   = 113
g0    = 9.80665
c     = 2.022

# TIME AXIS INPUT VECTOR DEFINITION
dt    = 0.01               # sec
T     = 200                # sec
t     = np.arange(0,T,dt)  # sec - check for lickage
N     = len(t)             # number of samples

#Selected turbulence input:
    # 1 for w1  = horizontal
    # 2 for w3  = vertical

windex = 1  # CHANGE THIS

# Number of Monte Carlo Iterations:
Niter  = 1
EPS    = 10e-12


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                       # HELPER FUNCTIONS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def plotting(yout,t,u,title):
    """
    Plots the time traces for the output variables [u, alpha, theta,q,n_z] according to a turbulence model

    The model is selected beforeahnd, when the yout response is simulated.

    """
    # Get each a/c state + the load factor

    plt.figure(str(title))

    y_u     = yout[:,0]           # -
    y_alpha = yout[:,1]       # rad
    y_theta = yout[:,2]       # rad
    y_q     = yout[:,3]      # rad/s
    y_nz    = yout[:,4]      # g

    plt.tight_layout()
    plt.subplot(5,1,1);
    plt.plot(t,y_u * V)
    plt.xlabel('t [s]');
    plt.ylabel(r'$u$ [m/s]');
    plt.title('Airspeed Deviation');
    plt.grid("True")

    # There might a mistake which resutls in weired AoA variations
    plt.subplot(5,1,2);
    plt.plot(t,y_alpha*180/np.pi)
    plt.xlabel('t [s]');
    plt.ylabel(r'$\alpha$ [deg]');
    plt.title('Angle of Attack');
    plt.grid("True")


    plt.subplot(5,1,3);
    plt.plot(t,y_theta*180/np.pi)
    plt.xlabel('t [s]');
    plt.ylabel(r'$\theta$ [deg]');
    plt.title('Pitch Angle');
    plt.grid("True")


    plt.subplot(5,1,4);
    plt.plot(t,y_q* 180/np.pi* (c/V))
    plt.xlabel('t [s]');
    plt.ylabel(r'$q$ [deg/s]');
    plt.title('Pitch Rate');
    plt.grid("True")


    plt.subplot(5,1,5);
    plt.plot(t,y_nz)
    plt.xlabel('t [s]');
    plt.ylabel(r'$n_z$ [deg/s]');
    plt.title('Load Factor');
    plt.grid("True")

    print(title + str('plotted \n\n '))




    #  MODEL IMPORT
cessna    = Cessna_model.Cessna()
model_ss  = cessna.state_space()




#     # INPUT VECTOR DEFINITION
# nn = np.zeros((1,N));                    # input elevator
# w  = np.random.randn(N)/np.sqrt(dt);     # scaled input hor. turbulence the sqrt(dt) because of lsim
u  = np.zeros((3,N))                     # input vector definition (vertical

    # Pass on the turbulence input (horizonatal/ vertical):
# u[windex,:] = w

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                # STABILITY & PRELIMIANRY CHECKS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# pole,zero = control.pzmap(model_ss,grid=True)
H_matrix   = cm.ss2tf(model_ss)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                # SIMULATION OF MOTION VARIABLES
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# ---------------------------------- MONTE CARLO LOOP ---------------------------------------------

ysum  = np.zeros((N,5))

for n in range(Niter):

        # REDEFINE INPUT VECTOR
    w  = np.random.randn(N)/np.sqrt(dt);     # scaled input by sqrt(dt) because of lsim
    u[windex,:] = w
    yaux = cm.lsim(model_ss,np.transpose(u),t)[0]       # tranpose u because that's how lsim wants it.....
    ysum = ysum + yaux

if Niter > 1:
    yout  = ysum/Niter
else:
    yout = ysum

if yout.shape[1] !=5:
    print("Number of model states retreived is: " + str(yout.shape[1]))


    # Get each a/c state + the load factor   -- NOT NEEDED HERE
# y          = {}
# y['u']     = yout[:,0]           # -
# y['alpha'] = yout[:,1]           # rad
# y['theta'] = yout[:,2]           # rad
# y['q']     = yout[:,3]           # rad/s
# y['nz']    = yout[:,4]           # g



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                # Plotting Time Results
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


if np.all(u[2,:] == 0. ) and np.all(u[0,:] == 0. ) :
    print("\n\n   Plot response to horizontal turbulence... \n\n")
    plotting(yout, t, u, title="Horizontal_turbulence")



if np.all(u[1,:] == 0. ) and np.all(u[0,:] == 0. ):
    print("\n\n   Plot response to vertical turbulence.. \n\n")
    plotting(yout, t, u, title="Vertical_turbulence")





# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                # SPECTRAL ANALYSIS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

W = 1        # noise intensity
W = W/dt     #CD scaling - this will be needed later for Lyapunov


# Define frequnecy axis:
freq  = (1/T) * np.arange(1,N//2)
omega = 2*np.pi*freq


remove_states = False          #KEEP THIS !!!
output_names  = ['u','alpha','theta','q','nz']

# Dictionaries to store tfsand power spectral functions
H           = {}   # transfer functions dict
y           = {}   # time response dict
Syy_ana     = {}   # spectral FUNCTION dict
Syy_exp     = {}   # experimental PERIODDOGRAM
Syy_expfilt = {}   # experimental FILTERED PERIODOGRAM

mag_ana     = {}
mag_exp     = {}
mag_filt    = {}

for i in range(len(output_names)):

    if remove_states:
        #might generate numerical instabilities
        aux  = cm.minreal(H_matrix[i,windex])     # [*,1] since we deal with w1(t) input
        num  = aux.num[0][0]
        rnum = num[np.abs(num) > EPS]

        H[output_names[i]] = (cm.tf(rnum,aux.den))
        print("States reduced: " + str(len(num) - len(rnum)))

        # print(" Old num: " +str(num), "\n Reduced num:" +str(rnum))

    else:

        # Isolate output variable from model SS:
        # this should be more robust
         # H[output_names[i]] = cm.minreal(cm.ss2tf(model_ss.A,model_ss.B[:,windex],model_ss.C[i,:],model_ss.D[i,windex]))

         # or by using H_matrix already defined in the previous section
         H[output_names[i]] = cm.minreal(H_matrix[i,windex])

    y[output_names[i]] = yout[:,i]


# -------------------------------- ANALYTICAL METHOD ----------------------------------------------


Suu = 1

for i in range(len(output_names)):
    # w,mag,phase = control.bode(H[output_names[i]],omega)
    # Syy_ana[output_names[i]] = mag **2
    Syy_ana[output_names[i]] =  H[output_names[i]] ** 2

        #signal.bode is different from cm.bode !!!!


# ------------------------------- EXPERIMENTAL METHOD ---------------------------------------------
# --------------------------------  (Monte Carlo)   -----------------------------------------------

Wexp = sp.fft.fft(w,N)

    # Generate PERIODOGRAMS ():

for i in range(len(output_names)):
        Y = sp.fft.fft(y[output_names[i]],N)
        Syy_exp[output_names[i]] = 1/N * Y * np.conjugate(Y)                # discrete
        Syy_exp[output_names[i]] = np.real(Syy_exp[output_names[i]]) * dt   # continuous


    # Generate FILTERED PERIODOGRAMS:

for i in range(len(output_names)):
    Syy_expfilt[output_names[i]]  = np.copy((Syy_exp[output_names[i]]))      #initialize with first version of the vector
    Syy_expfilt[output_names[i]][1:len(omega)-1] = 0.25 * Syy_exp[output_names[i]][:len(omega)-2]\
                                              +0.5 * Syy_exp[output_names[i]][1:len(omega)-1]\
                                         +0.25 * Syy_exp[output_names[i]][1:len(omega)-1]


    # Plotting

for i in range(len(output_names)):

    # Calcualte experimental Syy magnitude
    mag_exp[output_names[i]]  = np.absolute(Syy_exp[output_names[i]][1:N//2])    # experimental (no filter)
    mag_filt[output_names[i]] = np.absolute(Syy_expfilt[output_names[i]][1:N//2])   # WITH filter


    #Calcualte analytical Syy magnitude
    w,h = sp.signal.freqs(Syy_ana[output_names[i]].num[0][0],Syy_ana[output_names[i]].den[0][0],omega)
    mag_ana[output_names[i]] = np.absolute(np.real(h))

    plotname = str("Syy_" + output_names[i])
    plt.figure(plotname)
    plt.loglog(omega,mag_exp[output_names[i]],'-',  c = 'blue',  label = 'Experimental')
    plt.loglog(omega,mag_filt[output_names[i]],'-', c = 'green', label = 'Exp. Filtered')
    plt.loglog(omega,mag_ana[output_names[i]],'--', c = 'red'  , label = 'Analytical')

    plt.grid(True)
    plt.legend(loc='best')
    plt.xlabel(r"$\omega$ [rad/s]")
    plt.ylabel(r"$|S_{yy}|$ [$rad^2$/Hz]")
    # plt.ylim(EPS/10,10)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                    # VARIANCE CALCULATIONS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sigma_ana  = {}
sigma_freq = {}
sigma_time = {}


# ---------------------------------- 0 Lyapunov Equation ------------------------------------------
    # !!!!!!!!!! BBUGGED !!!!!
W   = 1
Bin = model_ss.B[:,windex]
Q   = np.dot(Bin,W*np.transpose(Bin))

# solve the eq for steady state Cxx
Cxx_ss = cm.lyap(model_ss.A,Q)

for i in range(len(output_names)-1):   #sigma-nz will not be calcualted
    sigma_ana[output_names[i]] = (Cxx_ss[i,i])


# ---------------------------------- 1 Power Spectral Density -------------------------------------

for i in range(len(output_names)):

    sigma_aux = (1/np.pi) * sp.integrate.simps(mag_ana[output_names[i]],omega)
    sigma_freq[output_names[i]] = sigma_aux


# ---------------------------------- 2 Transfer Function ------------------------------------------


# --------------------------------- 3 Time domanin Variance ---------------------------------------

for i in range(len(output_names)):
    sigma_time[output_names[i]] = np.var(y[output_names[i]])




print(sigma_ana,sigma_time,sigma_freq,sep='\n')  #time and freq seem to be close enough (same Order of Magnitude)

plt.show()


print("\n\n\n    Done! \n\n")