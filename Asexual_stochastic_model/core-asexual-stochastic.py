import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import os
from os.path import exists
from sklearn.linear_model import LinearRegression
import sys


# ============================================================
# Parameters that can be changed for each simulation run
# ============================================================

Vs = 1
B = 1.2
X = 25
dx = 0.1
dy = 0.1
dt = 0.2
T = 20
sigmax = 1
sigmay = 0.8**2


# ============================================================
# Command-line parameters
# ============================================================

c = float(sys.argv[1])
K = int(float(sys.argv[2]))
idata_macrofile = int(float(sys.argv[3]))


# ============================================================
# Computational grid
# ============================================================

Y = 4 + X * B
Kt = int(T / dt)
Kx = int(X / dx)
Ky = int(Y / dy)

dx = X / Kx
Nx = Kx + 1
dy = Y / Ky
Ny = Ky + 1
dt = T / Kt

x = np.linspace(0, X, Nx)
ix = np.arange(Nx)

y = np.linspace(0, Y, Ny)
iy = np.arange(Ny)

y = y - ((y[1] + y[-1]) / 2 - B * (x[1] + x[-1]) / 2)
x0ini = (x[1] + x[-1]) / 2

t = np.linspace(0, T, Kt + 1)
it = np.arange(Kt + 1)

x0 = np.copy(x)
y0 = np.copy(x)


# ============================================================
# Initial population
# ============================================================

# Vectorized equivalent of the original double loop.
Xgrid = x[:, None]
Ygrid = y[None, :]

n00 = (
    np.floor(
        np.exp(
            -(Xgrid - x0ini) ** 2 / 10
            - (Ygrid - B * x0ini) ** 2 / 10
        )
        * K
        * dy
    )
    - 50
)
y0ini=B*x0ini

n0 = np.where(n00 > 0, n00, 0)
In0 = np.sum(n0, axis=1)

n = np.copy(n0)
ntemp = np.copy(n)
In = np.copy(x)


# ============================================================
# Recording arrays and output files
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

# Kept from the original code.
Zn0 = 0 * np.copy(x)
for i in ix[1:len(ix)]:
    for j in iy[1:len(iy)]:
        Zn0[i] = max(
            0,
            min(
                Y,
                Zn0[i]
                + y[j] * n[i][j] / max(10 ** (-6), In0[i]),
            ),
        )

idatafile = 0
while exists(f"data/results{idatafile}.csv") == True:
    idatafile = idatafile + 1

f = open(f"data/results{idatafile}.csv", "x")
f.write("t,X,Xtheta,I,x1,y1\n")

f2 = open(f"data/coefficients{idatafile}.csv", "x")
f2.write("X,Y,dx,dy,dt,K,T,sigmax,sigmay,Vs,B,c,speed\n")
f2.write(
    f"{X},{Y},{dx},{dy},{dt},{K},{T},"
    f"{sigmax},{sigmay},{Vs},{B},{c},"
)


# ============================================================
# Dispersal kernels
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
# Precomputation of the two convolution operators
# ============================================================

# For a vector u of length N and a kernel g of length 2N-1,
# np.convolve(u, g, "valid") can be written as G @ u, where
#
#     G[j, i] = g[N - 1 + j - i].
#
# Since the kernels are fixed in time, these matrices only need
# to be built once.

row_y = np.arange(Ny)[:, None]
col_y = np.arange(Ny)[None, :]
gamma_y_matrix = gamma_y[Ny - 1 + row_y - col_y]

row_x = np.arange(Nx)[:, None]
col_x = np.arange(Nx)[None, :]
gamma_x_matrix = gamma_x[Nx - 1 + row_x - col_x]


# ============================================================
# Time loop
# ============================================================

for k in it:

    # --------------------------------------------------------
    # 1. Dispersal in trait and space
    # --------------------------------------------------------

    # Equivalent to:
    #
    #   for i:
    #       nconv_x[i] = np.convolve(n[i], gamma_y, "valid")
    #
    # followed by the analogous convolution in x.
    #
    # Matrix multiplication moves these operations into optimized
    # compiled linear-algebra routines.

    ntemp = n @ gamma_y_matrix.T
    ntemp = gamma_x_matrix @ ntemp

    # --------------------------------------------------------
    # 2. Local population sizes
    # --------------------------------------------------------

    In[1:-1] = np.sum(ntemp[1:-1, :], axis=1)

    # --------------------------------------------------------
    # 3. Stochastic demographic update
    # --------------------------------------------------------

    maladaptation = (
        y[None, 1:-1]
        - B * (x[1:-1, None] - c * t[k])
    )

    rate = (
        (1/(2*Vs)) * maladaptation**2
        + In[1:-1, None] / K
    )

    p_survival = np.exp(-rate * dt)

    p_birth = (
        np.exp((1 - rate) * dt)
        - p_survival
    )

    # np.random.binomial accepts a floating-point scalar number
    # of trials in the original loops by truncating it to an
    # integer.  The vectorized form requires the conversion to
    # be explicit.
    trials = ntemp[1:-1, 1:-1].astype(np.int64)

    n[1:-1, 1:-1] = (
        np.random.binomial(trials, p_survival)
        + np.random.binomial(trials, p_birth)
    )

    # --------------------------------------------------------
    # 4. Boundary conditions
    # --------------------------------------------------------

    In[0] = np.sum(ntemp[0])
    In[-1] = np.sum(ntemp[-1])

    n[0, 1:-1] = 0
    n[-1, 1:-1] = 0

    # --------------------------------------------------------
    # 5. Recording
    # --------------------------------------------------------

    Inonzeros = np.flatnonzero(In)

    Igrandtemp = np.where(In - K / 100 > 0, 1, 0)
    Igrand = np.flatnonzero(Igrandtemp)
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
    # 6. Shift of the spatial calculation window
    # --------------------------------------------------------

    if Xtrec3[k] > (0.3 * x[1] + 0.7 * x[-1]):
        truc = truc + 1
        shift = 5

        x = x + shift * dx

        ntemp = np.copy(n)
        n = np.zeros((Nx, Ny))

        # Vectorized equivalent of the original double loop.
        n[:-shift, :] = ntemp[shift:, :]

    if Xtrec3[k] < (0.5 * x[1] + 0.5 * x[-1]):
        truc = truc + 1
        shift = 5

        x = x - shift * dx

        ntemp = np.copy(n)
        n = np.zeros((Nx, Ny))

        # Vectorized equivalent of the first original double loop.
        n[shift:, :] = ntemp[:-shift, :]

        # Special refill rule at the left boundary.
        # The same range ix[0:shift-1] is retained.
        for i in ix[0:shift - 1]:
            trait_shift = int(i * B * dx / dy)
            stop = len(iy) - trait_shift

            if trait_shift == 0:
                n[i, :stop] = ntemp[shift, :stop]
            else:
                n[i, :stop] = np.maximum(
                    ntemp[shift, trait_shift:],
                    ntemp[shift, :stop],
                )

    # --------------------------------------------------------
    # 7. Shift of the trait calculation window
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

        n[:, :-yshift] = ntemp[:, yshift:]

    if deltay < -1:
        yshift = -5

        y = y + yshift * dy

        ntemp = np.copy(n)
        n = np.zeros((Nx, Ny))

        n[:, -yshift:] = ntemp[:, :yshift]


# ============================================================
# Calculation and recording of propagation speed
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
    f"{X},{Y},{dx},{dy},{dt},{K},{T},"
    f"{sigmax},{sigmay},{Vs},{B},{c},{model.coef_[0]}\n"
)

print(
    f"{X},{Y},{dx},{dy},{dt},{K},{T},"
    f"{sigmax},{sigmay},{Vs},{B},{c},{model.coef_[0]}"
)

f.close()
f2.close()
f3.close()
