import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import os
from os.path import exists
from sklearn.linear_model import LinearRegression
import sys
import csv


# ============================================================
# Parameters that can be changed for each simulation run
# ============================================================

Vs = 1
B = 1.2
X = 100
dx = 0.2
dy = 0.2
dt = 1
T = 40
sigmax = 1
sigmay = 1 * 0.8**2


# ============================================================
# Command-line parameters
# ============================================================

c = float(sys.argv[1])
idata_macrofile = int(float(sys.argv[2]))


# ============================================================
# Definition of the computational grid
# ============================================================

Y = 4 + X * B
Kt = int(T / dt)
Kx = int(X / dx)
Ky = int(Y / dy)

dx = X / Kx
Nx = Kx + 1
dy = Y / Ky
Ny = Ky + 1

x = np.linspace(0, X, Nx)
ix = np.linspace(0, Nx - 1, Nx).astype(int)

y = np.linspace(0, Y, Ny)
iy = np.linspace(0, Ny - 1, Ny).astype(int)

y = y - ((y[1] + y[-1]) / 2 - B * (x[1] + x[-1]) / 2)
x0ini = 0.7 * x[1] + 0.3 * x[-1]
y0ini= B*x0ini

t = np.linspace(0, T, Kt + 1)
it = np.linspace(0, Kt, Kt + 1).astype(int)

x0 = np.copy(x)
y0 = np.copy(x)


# ============================================================
# Definition of the initial population
# ============================================================

# Vectorized version of the original double loop.
Xgrid = x[:, None]
Ygrid = y[None, :]

n00temp = (
    np.exp(
        -(Xgrid - x0ini) ** 2 / 2
        - (Ygrid - B * Xgrid) ** 2 / (2 * 0.74)
    )
    / (np.sqrt(2 * np.pi * 0.74) * 5)
)
n00temp=n00temp/max(np.sum(np.abs(n00temp), axis=1) * dy)

n0 = np.copy(n00temp)
n = np.copy(n0)
ntemp = np.copy(n)
In = np.copy(x)
In0 = np.copy(x)


# ============================================================
# Initiation of recording and data files
# ============================================================

os.makedirs("data", exist_ok=True)

X0rec = np.copy(t)
Xtrec1 = np.copy(t)
Xtrec2 = np.copy(t)
Xtrec3 = np.copy(t)
itrec2 = np.copy(t)
x0rec = np.copy(t)
Nrec = np.copy(t)

Inonzeros0 = np.flatnonzero(In0)
truc = 0


Zn0 = 0 * np.copy(x)
for i in ix[1:len(ix)]:
    for j in iy[1:len(iy)]:
        Zn0[i] = max(
            0,
            min(
                Y,
                Zn0[i] + y[j] * n[i][j] / max(10 ** (-6), In0[i]),
            ),
        )

idatafile = 0
while exists(f"data/results{idatafile}.csv") == True:
    idatafile = idatafile + 1

f = open(f"data/results{idatafile}.csv", "x")
f.write("t,X,Xtheta,I,x1,y1\n")

f2 = open(f"data/coefficients{idatafile}.csv", "x")
f2.write("X,Y,dx,dy,dt,T,sigmax,sigmay,Vs,B,c,speed\n")
f2.write(f"{X},{Y},{dx},{dy},{dt},{T},{sigmax},{sigmay},{Vs},{B},{c},")


# ============================================================
# Definition of the dispersal kernels
# ============================================================

sigmakernelx = sigmax * dt
sigmakernely = sigmay * dt

xconv = np.linspace(-X, X, 2 * Nx - 1)
gamma_x = np.exp(-np.square(xconv) / (2 * sigmakernelx))
gamma_x = gamma_x / sum(gamma_x)

yconv = np.linspace(-Y, Y, 2 * Ny - 1)
gamma_y = np.exp(-np.square(yconv) / (2 * sigmakernely))
gamma_y = gamma_y / sum(gamma_y)


# ============================================================
# Time loop
# ============================================================

for k in it:

    # --------------------------------------------------------
    # 1. Dispersal in trait and space
    #
    # The direct convolutions are intentionally retained here.
    # This preserves the original tail/support behaviour exactly.
    # --------------------------------------------------------

    nconv_x = np.empty_like(n)
    for i in ix:
        nconv_x[i] = np.convolve(n[i], gamma_y, "valid")

    temp = nconv_x.T.copy()
    for j in iy:
        temp[j] = np.convolve(temp[j], gamma_x, "valid")

    ntemp = temp.T

    # --------------------------------------------------------
    # 2. Local density and population update
    # --------------------------------------------------------

    # Same local-density calculation as in the original code.
    In[1:-1] = np.sum(ntemp[1:-1, :], axis=1) * dy

    # Vectorized maladaptation term.
    maladaptation = (
        y[None, 1:-1]
        - B * (x[1:-1, None] - c * t[k])
    )

    death_rate = (
        (1/(2*Vs)) * maladaptation**2
        + In[1:-1, None]
    )

    # The original expression
    #
    #   ntemp*exp(-D*dt)
    #   + ntemp*(exp((1-D)*dt) - exp(-D*dt))
    #
    # simplifies exactly to
    #
    #   ntemp*exp((1-D)*dt).
    #
    # Only the interior entries are updated, exactly as before.
    n[1:-1, 1:-1] = (
        ntemp[1:-1, 1:-1]
        * np.exp((1 - death_rate) * dt)
    )

    # --------------------------------------------------------
    # 3. Boundary conditions
    # --------------------------------------------------------

    # These two formulas are intentionally kept exactly as in
    # the original code (in particular, without a factor dy).
    In[0] = np.sum(ntemp[0])
    In[-1] = np.sum(ntemp[-1])

    n[0, 1:-1] = 0
    n[-1, 1:-1] = 0

    # --------------------------------------------------------
    # 4. Recording
    # --------------------------------------------------------

    Inonzeros = np.flatnonzero(In)
    Igrand = np.flatnonzero(In > np.max(In) / 2)

    pos_front = Igrand[-1]

    Xtrec1[k] = x[Inonzeros[1]]
    Xtrec2[k] = x[Inonzeros[-1]]
    Xtrec3[k] = x[pos_front]
    Nrec[k] = np.sum(In)
    X0rec[k] = x[int(len(ix) / 2)]
    itrec2[k] = Inonzeros[-1]
    x0rec[k] = truc

    f.write(
        f"{t[k]}, {Xtrec2[k]}, {Xtrec3[k]},"
        f"{sum(In) * dx},{B * (x[1] + x[-1]) / 2},"
        f"{(y[1] + y[-1]) / 2}\n"
    )

    # --------------------------------------------------------
    # 5. Shift of the calculation window as the front moves
    # --------------------------------------------------------

    if Xtrec3[k] > (0.3 * x[1] + 0.7 * x[-1]):
        truc = truc + 1
        shift = 5

        x = x + shift * dx

        ntemp = np.copy(n)
        n = np.zeros((Nx, Ny))

        # Original double loop replaced by an array slice.
        n[:-shift, :] = ntemp[shift:, :]

    if Xtrec3[k] < (0.5 * x[1] + 0.5 * x[-1]):
        truc = truc + 1
        shift = 5

        x = x - shift * dx

        ntemp = np.copy(n)
        n = np.zeros((Nx, Ny))

        # Original double loop replaced by an array slice.
        n[shift:, :] = ntemp[:-shift, :]

        # Special refill rule at the left boundary.
        # The range is kept identical to ix[0:shift-1] in the
        # original code.
        for i in ix[0:shift - 1]:
            trait_shift = int(i * B * dx / dy)
            stop = len(iy) - trait_shift

            if trait_shift == 0:
                n[i, :stop] = np.maximum(
                    ntemp[shift, :stop],
                    ntemp[shift, :stop],
                )
            else:
                n[i, :stop] = np.maximum(
                    ntemp[shift, trait_shift:],
                    ntemp[shift, :stop],
                )

    # --------------------------------------------------------
    # 6. Shift of the trait window
    # --------------------------------------------------------

    deltay = (
        B * ((x[1] + x[-1]) / 2 - c * t[k])
        - (y[-1] + y[1]) / 2
    )

    if deltay > 1:
        yshift = 5

        y = y + yshift * dy

        ntemp = np.copy(n)
        n = np.zeros((Nx, Ny))

        # Original double loop replaced by an array slice.
        n[:, :-yshift] = ntemp[:, yshift:]

    if deltay < -1:
        yshift = -5

        y = y + yshift * dy

        ntemp = np.copy(n)
        n = np.zeros((Nx, Ny))

        # For yshift = -5 this means:
        # n[:, 5:] = ntemp[:, :-5]
        n[:, -yshift:] = ntemp[:, :yshift]


# ============================================================
# Calculation and recording of the propagation speed
# ============================================================

iinf = int(len(t) / 2)

tcalc = t[iinf:len(t)]
Xcalc = Xtrec3[iinf:len(t)]

xlin = tcalc.reshape((-1, 1))
ylin = Xcalc

model = LinearRegression().fit(xlin, ylin)
f2.write(f"{model.coef_[0]}")


# ============================================================
# Recording in the macro file
# ============================================================

f3 = open(f"data/data_macro{idata_macrofile}.csv", "a")
f3.write(
    f"{X},{Y},{dx},{dy},{dt},{T},"
    f"{sigmax},{sigmay},{Vs},{B},{c},{model.coef_[0]}\n"
)

print(
    f"{X},{Y},{dx},{dy},{dt},{T},"
    f"{sigmax},{sigmay},{Vs},{B},{c},{model.coef_[0]}"
)

f.close()
f2.close()
f3.close()


# ============================================================
# Final population statistics
# ============================================================

Zn = np.copy(x)
Vn = np.copy(x)

# Same formulas as the original code, vectorized over space.
In[1:-1] = np.sum(n[1:-1, :], axis=1) * dy

Zn[1:-1] = (
    np.sum(n[1:-1, :] * y[None, :], axis=1) * dy
    / (In[1:-1] + 10 ** (-10))
)

Vn[1:-1] = (
    np.sum(
        n[1:-1, :]
        * (y[None, :] - Zn[1:-1, None]) ** 2,
        axis=1,
    )
    * dy
    / (In[1:-1] + 10 ** (-10))
)

small_population = In[1:-1] < 10 ** (-7)
Vn[1:-1][small_population] = 0
Zn[1:-1][small_population] = 0


# ============================================================
# Export of the final profiles
# ============================================================

print(x0ini)
chose=x[1]-x0ini
print(chose)

with open("data/vectorscT40.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["x", "In", "Zn", "Vn"])

    for xvalue, Ivalue, Zvalue, Vvalue in zip(x-x0ini, In, Zn-y0ini, Vn):
        writer.writerow([xvalue, Ivalue, Zvalue, Vvalue])
