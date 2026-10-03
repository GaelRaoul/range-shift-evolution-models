import subprocess
import sys
from pathlib import Path

import numpy as np


print("c,vS,vN,vZS,vZN,Ik1,Ik2,lambda,A,B,C,tfinal")


# Parameters that can be changed for each simulation run
# These defaults reproduce the original sexual sweep:
# c = 1.68, 1.80, ..., 2.40.
c0 = 0.25
cmax = 0.25
nb = 1


# Vector of the different climate change speeds that will be used
# to launch simulations.
k_macro = np.linspace(0, nb - 1, nb).astype(int)
c_macro = np.linspace(c0, cmax, num=nb)


# Opening a file to write the results of the simulations in,
# without overwriting previous files.
data_dir = Path("data")
data_dir.mkdir(parents=True, exist_ok=True)

idata_macrofile = 0
while (data_dir / f"data_macro_sexual{idata_macrofile}.csv").exists():
    idata_macrofile += 1

macro_file = data_dir / f"data_macro_sexual{idata_macrofile}.csv"
with macro_file.open("x") as f:
    f.write("c,vS,vN,vZS,vZN,Ik1,Ik2,lambda,A,B,C,tfinal\n")


# Launch of the simulations
core_script = Path(__file__).with_name("core-sexual-simulations.py")

for k in k_macro:
    command = [
        sys.executable,
        str(core_script),
        str(c_macro[k]),
        str(idata_macrofile),
    ]

    print("Launching:", " ".join(command))
    subprocess.run(command, check=True)


print(f"Macro summary written to {macro_file}")
