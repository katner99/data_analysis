import matplotlib.pyplot as plt  
import numpy as np


# MEP practical 5 diffusion
# Solve diffusion equation using FTCS
# Jahnavi Millan 100371998


# Set model parameters and constants
dt = 240            # time step (s)
n_steps = 24*15     #number of time steps (15 steps = 1 hour)
J = 10              # number of soil levels
D = 1e-7            # Diffusion coefficient (m^2 s^(-1))
soil_depth = 0.1    # soil depth (m)
dz = soil_depth/J   # depth increment (m)
DD = D*dt/dz*dz     # calculate coefficient for FTCS scheme


# set Bounday Conditions
Tair = 40
Trock = 10


# Set up verticle z grid (j index), don't include BCs
z = np.zeros(J)
ja = np.arange(J)
z[ja] = (ja + int(0.5))*dz
print(z)


# set up an array for time grid, note add 1 to account for start time at 0
na = np.arange(n_steps+1)      


# set up T(n,j) array of zeros, where n index is for time, j index is for depth
# n starts at zero and j does not include BCs
# each row is a time step, and columns are soil layers
T = np.zeros((n_steps+1, J))
t = np.zeros(n_steps+1)


# initial condition is all soil layers have temperature Trock
T[0,:] = Trock

for n in range(n_steps):
#for n in range(2):
    for j in range(0,J-1):
        T[n+1,0] = T[n,0] + DD*(T[n,1]-3*T[n,0]+2*Tair)
        T[n+1,J-1] = T[n,J-1] + DD*(2*Trock-3*T[n,J-1]+T[n,J-2])
        T[n+1,j] = T[n,j] + DD*(T[n,j+1]-2*T[n,j]+T[n,j-1])
        t[n+1] = t[n] + dt
    #print(T)


  
# Plot Temperature versus depth through the soil
plt.plot(T[15,:],z,'g-o', label='1 hour')
plt.plot(T[15*2,:],z,'b-o', label='2 hours')
plt.plot(T[15*24,:],z,'r-o', label='24 hours')
plt.xlabel('soil temperature ($^o$ C)')
plt.ylabel('depth (m)')
plt.legend()
plt.axis([9, 10.5, 0, 0.05])
plt.gca().invert_yaxis()
plt.show()

