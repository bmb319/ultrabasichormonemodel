## Function of monoamine transporter in vesicles.
#DOI:
# Transport of histamine from cytosol to vesicles
# b = cha
# c = vha
def VMATH(b, c):
  km = 24        
  vmax =  8500 #10552
  a = (vmax*(b/(km + b)) - 5.*c)
  return a