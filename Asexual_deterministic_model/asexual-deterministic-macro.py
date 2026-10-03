import os
import numpy as np
from os.path import exists
import sys


os.makedirs("data", exist_ok=True)

print("X,Y,dx,dy,dt,T,sigmax,sigmay,Vs,B,c,speed")


# Parameters that can be changed for each simulation run
c0 = 0.65
cmax = 0.65
nb = 1


# Vector of the different climate change speeds that will be used
# to launch simulations.
k_macro = np.linspace(0, nb - 1, nb).astype(int)
c_macro = np.linspace(c0, cmax, num=nb)


# Opening a file to write the results of the simulations in,
# without overwriting previous files.
idata_macrofile = 0
while exists(f"data/data_macro{idata_macrofile}.csv") == True:
    idata_macrofile = idata_macrofile + 1

f = open(f"data/data_macro{idata_macrofile}.csv", "x")
f.write("X,Y,dx,dy,dt,T,sigmax,sigmay,Vs,B,c,speed\n")
f.close()


# Launch of the optimized simulations
for k in k_macro:
    command = (
        "python3 core-asexual-deterministic.py "
        + str(c_macro[k])
        + "  "
        + str(idata_macrofile)
    )

    print(command)
    os.system(command)
