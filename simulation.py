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
import Cessna_model
import numpy.random
import matplotlib.pyplot as plt
# matplotlib.rcParams['text.usetex'] = True

plt.rc('font', size=14) #controls default text size
plt.rc('axes', titlesize=18) #fontsize of the title
plt.rc('axes', labelsize=14) #fontsize of the x and y labels
plt.rc('xtick', labelsize=12) #fontsize of the x tick labels
plt.rc('ytick', labelsize=12) #fontsize of the y tick labels
plt.rc('legend', fontsize=13) #fontsize of the legend



np.set_printoptions(precision=3)
# np.random.seed(3)

  #  Data and constants:
V     = 59.9
m     = 6035
muc   = 113
g0    = 9.80665
c     = 2.022

# TIME AXIS INPUT VECTOR DEFINITION
dt    = 0.01               # sec
T     = 10000         # sec
t     = np.arange(0,T,dt)  # sec - check for lickage
N     = len(t)             # number of samples

# Selected turbulence input:
    # 1 for w1  = horizontal
    # 2 for w3  = vertical

turb_lst = ['none','horizontal','vertical']
windex       = 1            # CHANGE THIS
plottingflag = True         # Switch for plotting in time (False) or frequency (True)
combined_plot= True

# Number of Monte Carlo Iterations:
Niter  = 1
EPS    = 10e-12

print("Turbulence modelled: " + str(turb_lst[windex]))
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                       # HELPER FUNCTIONS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def plotting(yout_lst,t,u,title,figure=None):
    """
    Plots the time traces for the output variables [u, alpha, theta,q,n_z] according to a turbulence model

    The model is selected beforeahnd, when the yout response is simulated.

    """
    # Get each a/c state + the load factor

    fig = plt.figure(str(title))

    y_u     = yout[:,0]           # -
    y_alpha = yout[:,1]       # rad
    y_theta = yout[:,2]       # rad
    y_q     = yout[:,3]      # rad/s
    y_nz    = yout[:,4]      # g

    plt.tight_layout()
    plt.subplot(5,1,1);
    plt.plot(t,y_u)
    # plt.xlabel('t [s]');
    plt.ylabel(r'$\hat{u}$ [-]');
    # plt.title('Airspeed Deviation');
    plt.grid("True")
    plt.xlim(0,T)

    # There might a mistake which resutls in weired AoA variations
    plt.subplot(5,1,2);
    plt.plot(t,y_alpha* 180/np.pi)
    # plt.xlabel('t [s]');
    plt.ylabel(r'$\alpha$ [deg]');
    # plt.title('Angle of Attack');
    plt.grid("True")
    plt.xlim(0,T)


    plt.subplot(5,1,3);
    plt.plot(t,y_theta* 180/np.pi)
    # plt.xlabel('t [s]');
    plt.ylabel(r'$\theta$ [deg]');
    # plt.title('Pitch Angle');
    plt.grid("True")
    plt.xlim(0,T)


    plt.subplot(5,1,4);
    plt.plot(t,y_q * 180/np.pi)
    # plt.xlabel('t [s]');
    plt.ylabel(r'$q$ [deg]');
    # plt.title('Pitch Rate');
    plt.grid("True")
    plt.xlim(0,T)


    plt.subplot(5,1,5);
    plt.plot(t,y_nz)
    plt.xlabel('t [s]');
    plt.ylabel(r'$n_z$ [-]');
    # plt.title('Load Factor');
    plt.grid("True")
    plt.xlim(0,T)

    print(title + str(' plotted \n\n '))

    return fig

def plotting_combined(yout_h,yout_v,t,u,title,figure=None):
    """
    Plots the time traces for the output variables [u, alpha, theta,q,n_z] according to a turbulence model

    The model is selected beforeahnd, when the yout response is simulated.

    """
    # Get each a/c state + the load factor

    fig = plt.figure(str(title))

    y_uh     = yout_h[:,0]           # -
    y_alphah = yout_h[:,1]       # rad
    y_thetah = yout_h[:,2]       # rad
    y_qh     = yout_h[:,3]      # rad/s
    y_nzh    = yout_h[:,4]      # g

    y_uv     = yout_v[:,0]           # -
    y_alphav = yout_v[:,1]       # rad
    y_thetav = yout_v[:,2]       # rad
    y_qv     = yout_v[:,3]      # rad/s
    y_nzv    = yout_v[:,4]      # g

    plt.tight_layout()

    plt.subplot(5,1,1)
    plt.plot(t,y_uh,label="Horizontal",c='b')
    plt.plot(t,y_uv,label="Vertical",c='orange')
    # plt.xlabel('t [s]');
    plt.ylabel(r'$\hat{u}$ [-]');
    # plt.title('Airspeed Deviation');
    plt.grid("True")
    plt.legend()
    plt.xlim(0,T)

    # There might a mistake which resutls in weired AoA variations
    plt.subplot(5,1,2)
    plt.plot(t,y_alphah* 180/np.pi,c='b')
    plt.plot(t,y_alphav* 180/np.pi,c='orange')
    # plt.xlabel('t [s]');
    plt.ylabel(r'$\alpha$ [deg]');
    # plt.title('Angle of Attack');
    plt.grid("True")
    plt.xlim(0,T)


    plt.subplot(5,1,3)
    plt.plot(t,y_thetah* 180/np.pi,c='b')
    plt.plot(t,y_thetav* 180/np.pi,c='orange')
    # plt.xlabel('t [s]');
    plt.ylabel(r'$\theta$ [deg]');
    # plt.title('Pitch Angle');
    plt.grid("True")
    plt.xlim(0,T)


    plt.subplot(5,1,4)
    plt.plot(t,y_qv * 180/np.pi,c='orange')
    plt.plot(t,y_qh * 180/np.pi,c='b')

    # plt.xlabel('t [s]');
    plt.ylabel(r'$q$ [deg]');
    # plt.title('Pitch Rate');
    plt.grid("True")
    plt.xlim(0,T)


    plt.subplot(5,1,5)
    plt.plot(t,y_nzv,c='orange')
    plt.plot(t,y_nzh,c='b')
    plt.xlabel('t [s]');
    plt.ylabel(r'$n_z$ [-]');
    # plt.title('Load Factor');
    plt.grid("True")
    plt.xlim(0,T)

    print(title + str(' plotted \n\n '))

    return fig
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                       #  MODEL IMPORT
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


cessna    = Cessna_model.Cessna()
model_ss  = cessna.state_space()


H_matrix   = cm.ss2tf(model_ss)
model_LTI  = sp.signal.lti(model_ss.A,model_ss.B,model_ss.C,model_ss.D)
model_dLTI = sp.signal.dlti(model_LTI,dt)

#     # INPUT VECTOR DEFINITION
# nn = np.zeros((1,N));                    # input elevator
w   = np.random.randn(N)/np.sqrt(dt);     # scaled input hor. turbulence the sqrt(dt) because of lsim
u   = np.zeros((3,N))                       # input vector definition (vertical
u_h = np.zeros((3,N))
u_v = np.zeros((3,N))

    # Pass on the turbulence input (horizonatal/ vertical):
u[windex,:] = w

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                # STABILITY & PRELIMIANRY CHECKS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# pole,zero = control.pzmap(model_ss,grid=True)

# H_matrix [output, input]


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                # SIMULATION OF MOTION VARIABLES
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# ---------------------------------- MONTE CARLO LOOP ---------------------------------------------
if Niter > 1:
    ysum  = np.zeros((N,5))

    for n in range(Niter):

            # REDEFINE INPUT VECTOR

        yaux = cm.lsim(model_ss,np.transpose(u),t)[0]       # tranpose u because that's how lsim wants it.....
        ysum = ysum + yaux

    yout  = ysum/Niter
else:
    w  = np.random.randn(N)/np.sqrt(dt)     # scaled input by sqrt(dt) because of lsim
    u[windex,:] = w
    yout = cm.lsim(model_ss,np.transpose(u),t)[0]


# Just once - combined (for plotting)
u_h[1,:] = w
yout_h = cm.lsim(model_ss,np.transpose(u_h),t)[0]

u_v[2,:] = w
yout_v = cm.lsim(model_ss,np.transpose(u_v),t)[0]

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

if not plottingflag:
    if np.any(u[1,:] != 0. ) and np.all(u[0,:] == 0. ) :
        print("\n\n   Plot response to horizontal turbulence... \n\n")
        plotting(yout, t, u, title="Horizontal_turbulence")

    if np.any(u[2,:] != 0. ) and np.all(u[0,:] == 0. ):
        print("\n\n   Plot response to vertical turbulence.. \n\n")
        plotting(yout, t, u, title="Vertical_turbulence")

    if combined_plot:
        print("\n\n   Plot combined response to turbulence.. \n\n")
        fig = plotting_combined(yout_h,yout_v, t, u, title="Combined")



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                # SPECTRAL ANALYSIS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

W = 1        # noise intensity
W = W/dt     # CD scaling - this will be needed later for Lyapunov


#      Define frequnecy axis:
freq  = (1/T) * np.arange(1,N)
omega = 2*np.pi*freq


remove_states = False        #KEEP THIS !!!
output_names  = ['u','alpha','theta','q','nz']

#    Dictionaries to store tfsand power spectral functions
H           = {}   # transfer functions dict
y           = {}   # time response dict
Syy_ana     = {}   # spectral FUNCTION dict
Syy_exp     = {}   # experimental PERIODDOGRAM
Syy_expfilt = {}   # experimental FILTERED PERIODOGRAM

mag_ana     = {}
mag_exp     = {}
mag_filt    = {}

h  = {}    # cross transfer fucntion ??

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
         H[output_names[i]] = (H_matrix[i,windex])

    y[output_names[i]] = yout[:,i]


# -------------------------------- ANALYTICAL METHOD ----------------------------------------------


Suu = 1

for i in range(len(output_names)):
    # w,mag,phase = control.bode(H[output_names[i]],omega)
    # Syy_ana[output_names[i]] = mag **2
    Syy_ana[output_names[i]] =  H[output_names[i]]

        #signal.bode is different from cm.bode !!!!


        #  Calcualte analytical Syy magnitude
for i in range(len(output_names)):
    w,mag = sp.signal.freqs(Syy_ana[output_names[i]].num[0][0],Syy_ana[output_names[i]].den[0][0],omega)
    # mag = cm.bode(Syy_ana[output_names[i]],omega,Plot=False)[0]
    mag_ana[output_names[i]] = np.power((np.abs(mag)),2)


# ------------------------------- EXPERIMENTAL METHOD ---------------------------------------------


Wexp = sp.fft.fft(w,N)

    # Generate PERIODOGRAMS ():

for i in range(len(output_names)):
        Y = sp.fft.fft(y[output_names[i]],N)
        Syy_exp[output_names[i]] = 1/N * Y * np.conjugate(Y)                # discrete
        Syy_exp[output_names[i]] = np.real(Syy_exp[output_names[i]]) * dt   # continuous

        #np.real is taken to force the variable to be lsited as object dtype = float for futre plotting
        # the imaginary part should be and it is (verified) 0 for all omega


    # Generate FILTERED PERIODOGRAMS:

for i in range(len(output_names)):
    Syy_expfilt[output_names[i]]  = np.copy((Syy_exp[output_names[i]]))      # initialize with first version of the vector
    Syy_expfilt[output_names[i]][1:len(omega)-1] = 0.25 * Syy_exp[output_names[i]][:len(omega)-2]\
                                              +0.5 * Syy_exp[output_names[i]][1:len(omega)-1]\
                                         +0.25 * Syy_exp[output_names[i]][1:len(omega)-1]




    # Plotting :

label_lst =[r'$S_{\hat{u} \hat{u}}$', r'$S_{\alpha \alpha}$', r'$S_{\theta \theta}$', r'$S_{~\frac{q \cdot c}{V}\frac{q \cdot c}{V}}$', r'$S_{n_{z} n_{z}}$']
unit_lst = [r'$~[\frac{1}{Hz}$]', r'$~\left [ \frac{rad^2}{Hz} \right ]$', r'$~\left [ \frac{rad^2}{Hz} \right ]$', r'$~\left [ \frac{rad^2}{Hz} \right ]$', r'$~\left [ \frac{1}{Hz} \right ]$']

if plottingflag:
    for i in range(len(output_names)):

        # Calcualte experimental Syy magnitude
        mag_exp[output_names[i]]  = np.absolute(Syy_exp[output_names[i]][1:N//2])    # experimental (no filter)
        mag_filt[output_names[i]] = np.absolute(Syy_expfilt[output_names[i]][1:N//2])   # WITH filter


        plotname = str("Syy_" + output_names[i] + '_h')
        plt.figure(plotname)
        plt.loglog(omega[1:N//2],mag_exp[output_names[i]],'-',  c = 'blue',  label = 'Experimental')
        plt.loglog(omega[1:N//2],mag_filt[output_names[i]],'-', c = 'green', label = 'Experimental - Filtered')
        plt.loglog(omega[1:],mag_ana[output_names[i]][1:],'--', c = 'red'  , label = 'Analytical')

        plt.grid(True)
        plt.legend(loc='lower center')
        plt.xlabel(r'$\omega~ \left [ \frac{rad}{s} \right ]$')
        plt.ylabel(label_lst[i] + unit_lst[i] )
        plt.tight_layout()
        # plt.ylim(EPS/10,10)



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                    # VARIANCE CALCULATIONS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sigma_lyap = {}
sigma_freq = {}
sigma_time = {}
sigma_tran = {}      # numerical integration of cross transfer fucntions

# ---------------------------------- 0 Lyapunov Equation ------------------------------------------
    # !!!!!!!!!! BBUGGED !!!!!
W   = 1
Wc  = W/dt
Bin = model_ss.B[:,windex]
Q   = np.dot(Bin,W*np.transpose(Bin))

# solve the eq for steady state Cxx
Cxx_ss = cm.lyap(model_ss.A,Q)

for i in range(len(output_names)-1):   #sigma-nz will not be calcualted
    sigma_lyap[output_names[i]] = (Cxx_ss[i,i])


# ---------------------------------- 1 Power Spectral Density -------------------------------------

for i in range(len(output_names)):
    sigma_aux = (1/np.pi) * sp.integrate.simps(mag_ana[output_names[i]],omega)
    sigma_freq[output_names[i]] = sigma_aux


# ---------------------------------- 2 Transfer Functions BUGGED DUE TO LSIM-----------------------------------------

#  CALCULATION OF PRODUCT MATRIX OF IMPULSE RESPONSES

    # define a longer time:
ext       = 1
textended = np.arange(0,T,dt)
x0        = model_ss.B[:,windex].reshape(7)
u_impulse = np.zeros((ext*N*ext,3))

    # Model impulse response:

h_matrix  = sp.signal.lsim(model_LTI,u_impulse,textended,x0)[1]
# h_matrix  = sp.signal.impulse(model_LTI,x0,textended)[1]


for i in range(len(output_names)):
        h[output_names[i]+output_names[i]] = np.power(h_matrix[:,i],2)
        # h[output_names[i]+output_names[i]]  = sp.signal.impulse(H[output_names[i]],x0,textended)[1]
        # h[output_names[i]+output_names[i]] = np.power(cm.lsim(H[output_names[i]],u_impulse[:,0],textended,x0[:5])[0],2)


      # Initilization of sigmas:
for k in range(len(output_names)):
      sigma_tran[output_names[k]+output_names[k]] = 0

      # Manual integration:
for epoch in range(N):
    for k in range(len(output_names)):
        sigma_tran[output_names[k]+output_names[k]]+= dt*h[output_names[k]+output_names[k]][epoch]

    # SciPy integration:
# for k in range(len(output_names)):
#         sigma_tran[output_names[k]+output_names[k]] = sp.integrate.simps(h[output_names[k]+output_names[k]],textended)

# Cuu_lst = np.zeros((N))
# for n in range(1,N):
#     Cuu_lst[n] = h[output_names[2]+output_names[2]][n]*dt + Cuu_lst[n-1]


# --------------------------------- 3 Time Domanin Variance ---------------------------------------

for i in range(len(output_names)):
    sigma_time[output_names[i]] = np.var(y[output_names[i]])


avg_sig_lst = []


# --------------------------------- PRINTING VARIANCES ---------------------------------------

print("Variable \t Time     PSD  ")

for i in range(len(output_names)):

    avg_sig= (sigma_time[output_names[i]]+sigma_freq[output_names[i]])/2

    print(output_names[i],end=' ')
    # print( " %.3f  %.3f  %.3f" %(sigma_time[output_names[i]]/avg_sig,sigma_freq[output_names[i]]/avg_sig,sigma_tran[output_names[i]+output_names[i]]/avg_sig))
    print( "\t %.5e \t %.3e  " %(sigma_time[output_names[i]],sigma_freq[output_names[i]]))

    avg_sig_lst.append(avg_sig)

plt.show()


print("\n\n\n    Done! \n\n")