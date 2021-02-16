# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 14:06:55 2021

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
import Cessna_model


def foo(x):
    return 1-x*x

def myfilter(y):

    y_filt2 = np.zeros(len(y))
    y_filt2[0] = y[0]
    y_filt2[-1] = y[-1]
    for i in range(1,len(y)-1):

        y_filt2[i] = 0.25 * y[i-2] + 0.5 * y[i] + 0.25 * y[i+1]
    return y_filt2


# N = 10
# x = np.arange(0,N,0.1)
# y = foo(x) + np.random.randn(100)

# y_filt1 = np.zeros(len(y))
# y_filt1[0] = y[0]
# y_filt1[-1] = y[-1]
# y_filt1[1:len(y)-1] = 0.25 * y[:len(y)-2] + 0.5 * y[1:len(y)-1] + 0.25 * y[1:len(y)-1]


# # y_filt2 = y
# # for i in range(1,len(y)-1):
# #     y_filt2[i] = 0.25 * y[i-2] + 0.5 * y[i] + 0.25 * y[i+1]

# y_filt2 = myfilter(y)
# plt.plot(x,y)
# plt.plot(x,y_filt1)
# # plt.plot(x,y_filt2)




# Testing the model:


cessna = Cessna_model.Cessna()
testmodel = cessna.state_space()

print(cm.pole(testmodel))
