## Function of histidine decarboxylase
#DOI:
# Histamine synthesis in cytosol.
# UNITS in uM/h. 
def VHTDC(b):
  # b = cht
  #  c = G*
  km = 270
  vmax = 233 #232
  a =  vmax*(b/(km + b))
  return a