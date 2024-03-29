
from scipy.integrate._ivp.common import validate_max_step
import numpy as np
import pandas as pd
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
import math 
from comp_model import *

#the equations model the differential terms of this model, calculate how much these variables are going to change in each iterations
#find the ks that give the diff eq to 0

#Time array.
t_factor = 1 # Time factor for graphs.
time = 100/t_factor # Time of simulation depending on t_factor.
sampling_rate = 10*t_factor #number of samples per time factor units.
time_array = np.linspace(0, time, math.floor(time * sampling_rate + 1))


#Initial conditions
z0 = [3.1968, 140.3708, 1.4717,  100, 300, 298.31, 0.7114, 0.628942699789106, 1.0084, 0.7221, 1.3593, 1.0084, 150, 3,	140,	0.7205,1,1, 5.36792170e-02,  3.28039657e-03, 2.47918457e-02,  6.30086107e-04, 2.27306535e-04,  2.42721849e-04, 1] 
#Constant parameters for steady state

#buscar valores que tengan sentido y evidentemente sea positive
# Call odeint with the correct function and arguments
x = odeint(comp_model, z0, time_array)

#Get solution of the differential equation.
#sol = solve_ivp(com_model_basic, t_span = (time_array[0], time_array[-1]), t_eval = time_array, y0 = z0, method = 'RK45') #shows how the steady state values have changed over time. 


## Plot all variables.
#plt.figure(1)
#for i in range(0, 10):
#  plt.subplot(2, 5,i+1)
#  plt.plot(time_array, sol.y[i, :])
#  plt.ylabel(str(i))
#plt.show()


#plt.figure(0)
#plt.plot(time_array, sol.y[4, :])
#plt.show()


plt.figure(1)
plt.plot(time_array, x[:, 4])
#plt.ylim([0.99, 1.01])
plt.show()



#plt.figure(1)
#plt.subplot(6,1,1)
#plt.subplot(6,1,2)
#plt.plot(time_array, x[:, 20])
#plt.subplot(6,1,3)
#plt.plot(time_array, x[:, 21])
#plt.subplot(6,1,4)
#plt.plot(time_array, x[:, 22])
#plt.subplot(6,1,5)
#plt.plot(time_array, x[:, 23])
#plt.subplot(6,1,6)
#plt.plot(time_array, x[:, 24])
#plt.show()

#for i in range(0, 3):
#  plt.figure(i)
#  plt.plot(time_array, sol.y[6+i, :])
#  plt.show()

