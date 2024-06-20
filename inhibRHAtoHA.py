## Function of histamine inhibition of histamine firing.
#DOI:
# Dependent on Serotonin-activated g-protein (b) and the  equilibrium  value
# of  the  activated  G-protein (c).
# Units in uM. 
def inhibRHAtoHA(b, c):
  #b = gstar
  a = 1 - (0.2 * 3.5)*(b - c)     
  return a