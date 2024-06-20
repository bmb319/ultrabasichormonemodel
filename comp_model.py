import numpy as np
from activ_E2_to_ha_R_neuron import *
from activ_E2_to_ha_syn_neuron import *
from BE2in import *
from BPin import *
from fireha import *
from HTin import *
from inhib_cp_to_ERalpha import *
from inhibRHAtoHA import *
from inhibsynHAtoHA import *
from VHAT import *
from VHNMT import *
from VHTDC import *
from VHTL import *
from VMATH import *
from comp_model import *

z = np.zeros(25) 
def comp_model(t,z):
  dz = np.zeros(len(z)) 




  
## ------------------ Histamine neuron model -------------------------------
## Histamine neuron variables
#z[0] Cytosolic histamine (uM)
#z[1] Vesicular histamine (uM)
#z[2] Extracellular histamine (uM)
#z[3] Blood histidine (uM)
#z[4] Cytosolic histidine (uM)
#z[5] Cytosolic histidine pool (uM)
#z[6] Activated g-coupled protein H3 receptor (uM)
#z[7] Activated T protein H3 receptor (uM)
#z[8] Histamine bound to h3 receptor (uM)
#z[9] Activated g-coupled protein GPER (uM)
#z[10] Activated T protein GPER (uM)
#z[11] E2 bound to GPER (uM)

## Histamine neuron constants. 
  b1 = 15.013  #HA leakage from the cytosol to the extracellular space #(fitted)
  b2 =  (0.8)*8.355 #9.5086 #3.5   #HA release per action potential #(fitted) 
  b3 = 0.05  #HA removal from the extracellular space #(fitted) 
  b4 = 0.25  #Strength of stabilization of bHT to bHT0 #(fitted)
  b5 =  2.4875 #From cHT to HTpool. #(fitted)
  b6 = 1.424 #From HTpool to cHT. #(fitted)
  b7 =  1.083 #1.077585 #Other uses of HT remove HT.#(fitted)

  b8 = 100 #Histamine bound to autoreceptors produce G∗ #(fitted)
  b9 = 961.094 #T∗ facilitates the reversion of G∗ to G #(fitted)
  b10 = 8.792311243082075 #G∗ produces T∗ #(fitted)
  b11 = 74 #decay coefficient of T∗ #(fitted)
  b12 = 4.5  #eHA binds to autoreceptors #(fitted)
  b13 = 65.6196645577152 #eHA dissociates from autoreceptors #(fitted)
  g0HH = 10  #Total g-coupled protein for H3 on HA neuron #(fitted)
  t0HH = 12.643209422616819 #Total T protein for H3 on HA neuron #(fitted)
  b0HH = 11  # 10 Total bound H3 receptors on HA neuron #(fitted)
  b14 = 99.9915191655649 #E2 bound to GPER produce G∗ in HA neuron #(fitted)
  b15 = 961.094 #T∗ in GPER facilitates the reversion of G∗ to G in HA neuron #(fitted)
  b16 = 20.002302363640318 #GPER G∗ produces T∗ in HA neuron #(fitted)
  b17 = 66.2992 #decay coefficient of T∗ in GPER #(fitted)
  b18 = 32374.73317631864  #Extracellular E2 binds to GPER receptors. # (fitted) 
  b19 = 65.61789 #Extracellular E2 dissociates from GPER receptors. # (fitted)
  g0GPER = 10  #Total g-coupled protein for GPER in HA neuron. #(fitted)
  t0GPER = 10 #Total T protein for GPER in HA neuron. #(fitted)
  b0GPER = 10  #Total bound E2 in GPER in HA neuron. #(fitted)

#Steady state values.
  gstar_ha_basal =  0.70955  #Equilibrium concentration of g* histamine in H3 receptor. #(fitted)
  gstar_E2_basal =   0.8029 #Equilibrium concentration of g* E2 in GPER receptor #(fitted)
  bht0 = 100 #Steady state value of blood histidine. #(fitted)
  basal_bound_ce2 =  0.955 #(fitted)



#Synapse
  dz[0] = inhibsynHAtoHA(z[6], gstar_ha_basal) * activ_E2_to_ha_syn_neuron(z[24], basal_bound_ce2) * VHTDC(z[4])  - VMATH(z[0], z[1]) -  VHNMT(z[0]) - b1*(z[0] - z[2]) + VHAT(z[2])

  dz[1] = VMATH(z[0], z[1]) - inhibRHAtoHA(z[6], gstar_ha_basal)*activ_E2_to_ha_R_neuron(z[9], gstar_E2_basal)*fireha(t)*b2*z[1]

  
  dz[2] = inhibRHAtoHA(z[6], gstar_ha_basal)*activ_E2_to_ha_R_neuron(z[9], gstar_E2_basal)*fireha(t)*b2*z[1] - VHAT(z[2]) + b1*(z[0] - z[2])  - b3*z[2] 


  
  dz[3] = HTin(t) - VHTL(z[3])  - b4*(z[3] - bht0) #(fitted)

  dz[4] = VHTL(z[3]) - inhibsynHAtoHA(z[6], gstar_ha_basal) * activ_E2_to_ha_syn_neuron(z[24], basal_bound_ce2) * VHTDC(z[4]) - b5*z[4] + b6*z[5] #(fitted)

  dz[5] = b5*z[4] - b6*z[5] - b7*z[5] #(fitted)

#autoreceptors
  dz[6]  = b8*z[8]**2*(g0HH - z[6]) - b9*z[7]*z[6] #(fitted)
  dz[7] = b10*z[6]**2*(t0HH - z[7])  - b11*z[7] #(fitted)
  dz[8] = b12*z[2]*(b0HH - z[8])  - b13*z[8] #(fitted)

#GPER
  dz[9]  = b14*z[11]**2*(g0GPER - z[9]) - b15*z[10]*z[9] #(fitted)
  dz[10] = b16*z[9]**2*(t0GPER - z[10])  - b17*z[10] #(fitted)
  dz[11] =  b18*z[22]*(b0GPER - z[11])  - b19 * z[11] #(fitted)

## ------------------- Estrogen and Progesterone ---------------------------------

#Variables
#z[18] Blood progesterone (uM) 5–20 ng/ml or 3–35 ng/mL, Molar mass = 314.469 g/mol, (17/314.469) = 0.054 uM
#z[19] Extracellular progesterone (uM) 1 ng/g
#z[20] Cytosolic progesterone (uM) from 5 to 10 pg/mg
#z[21] Blood E2 (uM) 109±15 pg/ml (30 to 400 pg/mL)
#z[22] Extracellular E2 (uM) 250 pg/ml?
#z[23] Cytosolic E2 (uM) 6.4±1 pg/mg
#z[24] Bound E2 to ER-alpha receptors in the nucleus of neurons (uM). 10^-10 M


  d1 = 0.50475 #Rate of cytosolic E2 binding to ER-alpha nuclear receptors. (fitted)
  d2 = 0.0001225 #Rate of unbonding E2 from ER-alpha nuclear receptors. (fitted)
  d3 = 0.12375542 #Rate of transport from blood progesterone <-> exytracellular progesterone (h^-1). (fitted)
  d4 = 2.00501032 #Rate of transport from blood progesterone <-> cytosolic progesterone (h^-1).  (fitted)
  d5 = 0.52607241 #Rate of transport from blood E2 <-> extracellular E2 (h^-1). (fitted)
  d6 = 0.97969309 #Rate of transport from blood E2 <-> cytosolic E2 (h^-1). (fitted)
  d7 = 1.54 #Rate of transport from extracelular progesterone <-> cytosolic progesterone (h^-1). (fitted)
  d8 = 1 #Rate of transport from extracelular E2 <-> cytosolic E2 (h^-1). (fitted)
  d9 = 12 #Rate of removal of extracellular progesterone (h^-1). (fitted)
  d10 = 1 #Rate of removal of cytosolic progesterone (h^-1). (fitted)
  d11 = 1 #Rate of removal of extracellular E2 (h^-1). (fitted)
  d12 = 1.5 #Rate of removal of cytosolic E2 (h^-1). (fitted)

  basal_cp = 0.02494 #Basal cytosolic progesterone (uM) #(fitted)
  basal_bp = 0.054 #Basal blood progesterone (uM) # (fitted)
  basal_be2 = 0.000636 #Basal blood E2 (uM) #(fitted)

  dz[18] = BPin(z[18], basal_bp) - d3 * (z[18] - z[19]) - d4 * (z[18] - z[20]) #(fitted)
  dz[19] = d3 * (z[18] - z[19]) + d7 * (z[20] - z[19]) - d9 * z[19] #(fitted)
  dz[20] = d4 * (z[18] - z[20]) - d7 * (z[20] - z[19]) - d10 * z[20] #(fitted)

  dz[21] = BE2in(z[21], basal_be2) - d5 * (z[21] - z[22]) - d6 * (z[21] - z[23]) #(fitted)
  dz[22] = d5 * (z[21] - z[22]) + d8 * (z[23] - z[22]) - d11 * z[22] #(fitted)
  dz[23] = d6 * (z[21] - z[23]) - d8 * (z[23] - z[22]) - d12 * z[23] #(fitted)
  dz[24] = d1 * inhib_cp_to_ERalpha(z[20], basal_cp) * z[23] - d2 * z[24] #(fitted)


  return dz

