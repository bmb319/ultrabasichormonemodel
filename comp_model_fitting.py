
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
#z0 = [3.1968, 142.55, 1.4717,  100, 300, 298.31, 0.88, 0.628942699789106, 1.0084, 0.9, 1.3593, 1.0084, 150, 3,	140,	0.7205,1,1, 5.36792170e-02,  3.28039657e-03, 2.47918457e-02,  6.30086107e-04, 2.27306535e-04,  2.42721849e-04, 1]


import numpy as np

#this equation will find the values of x so that dz is 0
#set the parameter that do not need to be fitted to values
#the goal is to find the parameters that we do not know yet lik b5, b.., d3.. these must always be positive
#then z[] values you substitute them for their values in equilibrium, which are found on "main" as the initial conditions


def comp_model_fitting(x): 

  dz = np.zeros(1) 
  z = np.array([1.4717])
  #fitting z2

  z0 = 3.1968
  z1 = 142.55
  z6 = 0.88
  z9 = 0.9
  

  b1 = 15.013  #HA leakage from the cytosol to the extracellular space (fitted)
  b2 = 8.355059385177153 #9.5086 #3.5   #HA release per action potential.
  b3 = 0.05  #HA removal from the extracellular space (fitted) 
  

  gstar_ha_basal =  0.70955  #Equilibrium concentration of g* histamine in H3 receptor. (fitted) 
  gstar_E2_basal =  0.8029051932417229 #0.7221 #Equilibrium concentration of g* E2 in GPER receptor.(fitted) 
  
 
  dz[0] = inhibRHAtoHA(z6, gstar_ha_basal)*activ_E2_to_ha_R_neuron(z9, gstar_E2_basal)*fireha(0)*b2*z1 - VHAT(z[0]) + b1*(z0 - z[0])  - b3*z[0] 


  #si hay una t valriable, ponemos t=0 


  return dz


#tried improving a dampening term, a feedback mechanism...
