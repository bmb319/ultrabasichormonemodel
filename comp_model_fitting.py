
from inhibsynHAtoHA import *
from VMATH import *
from VHNMT import *
from VHTDC import *
from VHAT import *
from VHTL import *
from HTin import *
from inhibRHAtoHA import *
from fireha import *
from activ_E2_to_ha_syn_neuron import *
from activ_E2_to_ha_R_neuron import *
from inhib_cp_to_ERalpha import *
from BPin import *
from BE2in import *

#Initial conditions
#z0 = [3.1968, 140.3708, 1.4717,  100, 300, 298.31, 0.7114, 0.628942699789106, 1.0084, 0.7221, 1.3593, 1.0084, 150, 3,	140,	0.7205,1,1, 5.36792170e-02,  3.28039657e-03, 2.47918457e-02,  6.30086107e-04, 2.27306535e-04,  2.42721849e-04, 1] 


import numpy as np

#this equation will find the values of x so that dz is 0
#set the parameter that do not need to be fitted to values
#the goal is to find the parameters that we do not know yet lik b5, b.., d3.. these must always be positive
#then z[] values you substitute them for their values in equilibrium, which are found on "main" as the initial conditions


def comp_model_fitting(x): 

  dz = np.zeros(1) 
  z = np.array([0.7114])
  #q = x[0]
  #fitting z6
  g0HH =  4.940277777777778 # 10
  z6 =  0.7114
  b8 = 100 #Histamine bound to autoreceptors produce G∗ 
  b9 = 961.094 #T∗ facilitates the reversion of G∗ to G
  b10 = 8.792311243082075 #G∗ produces T∗ #(fitted)
  b11 = 85 #66.2992 #decay coefficient of T∗ 
  z8 = 1.0084
  z7 = 0.628942699789106


  dz[0]  = b8*z8**2*(g0HH - z[0]) - b9*z7*z[0] #(fitted)






  #si hay una t valriable, ponemos t=0 


  return dz
