import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Constants
Vmax = 1.45  # μM/s
KM = 0.5 #0.39    # μM
Cp_max = 0.1 # Max release concentration (μM)
release_duration = 10 # Duration of release (s)

# Release term function (mountain-like shape)
def release(t):
    if t < release_duration:
        return Cp_max * (t / release_duration)
    else:
        return Cp_max * (1 - (t - release_duration) / release_duration)

# Differential equation for histamine concentration
def histamine_conc(C, t):
    dCdt = release(t) - (Vmax * C) / (KM + C)
    return dCdt

# Initial histamine concentration
C0 = 0

# Time points
t = np.linspace(0, 20, 1000)  # Define time points for simulation

# Solve differential equation
C = odeint(histamine_conc, C0, t)

# Plot results
plt.plot(t, C, label='Histamine Concentration')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (μM)')
plt.title('Histamine Release and Reuptake Kinetics')
plt.legend()
plt.grid(True)
plt.show()
