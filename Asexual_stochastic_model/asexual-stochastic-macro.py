import os
import numpy as np
from os.path import exists
import sys


os.makedirs("data", exist_ok=True)

print("X,Y,dx,dy,dt,K,T,sigmax,sigmay,Vs,B,c,speed")


# ============================================================
# Parameters that can be changed for each simulation run
# ============================================================

c0 = 0
cmax = 1
K_macro = 100000
nb = 20


# ============================================================
# Climate-change speeds used for the simulations
# ============================================================

k_macro = np.arange(nb, dtype=int)
c_macro = np.linspace(c0, cmax, num=nb)


# ============================================================
# Creation of the macro output file
# ============================================================

idata_macrofile = 0

while exists(f"data/data_macro{idata_macrofile}.csv") == True:
    idata_macrofile = idata_macrofile + 1

f = open(f"data/data_macro{idata_macrofile}.csv", "x")
f.write("X,Y,dx,dy,dt,K,T,sigmax,sigmay,Vs,B,c,speed\n")
f.close()


# ============================================================
# Launch of the simulations
# ============================================================

for k in k_macro:
    command = (
        "python3 core-asexual-stochastic.py "
        + str(c_macro[k])
        + " "
        + str(K_macro)
        + " "
        + str(idata_macrofile)
    )

    print(command)
    os.system(command)
