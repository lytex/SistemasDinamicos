import sympy as sym
from sympy.functions import sqrt
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np

plt.ion()

x0_e = 8.48528137423857
k, z = sym.symbols('k, z')


Jsym = sym.Matrix([[-10, 10, 0],
[1-k, -1, -sqrt(72)],
[sqrt(72), sqrt(72), sym.Rational(-8, 3)]])
eq = (Jsym - sym.Matrix.eye(3)*z).det()

def LorenzControl(σ=10, β=8/3, ρ=28, k=1):
    def f(x, t=0):
        y = np.empty_like(x)
        y[0] = σ*(x[1] - x[0])
        y[1] = x[0]*(ρ - x[2]) - x[1] - k*(x[0]-x0_e)
        y[2] = x[0]*x[1] - β*x[2]
        return y
    return f

# Valores máximos y mínimos de k
kmin = 0
kmax = 100
kvec = np.linspace(kmin, kmax, 100)

expr = sym.solve(eq, z)
f = sym.lambdify(k, expr, 'numpy')
xv = np.real(f(kvec))
yv = np.imag(f(kvec))

def rloc(kval):
    old_ax = ax1.axis()
    ax1.cla()
    ax1.plot([0, 0], [yv.min(), yv.max()], '-k')
    ax1.plot([xv.min(), xv.max()], [0, 0], '-k')
    ax1.grid(True)
    color = ['r', 'g', 'b'] # Código de color para las raíces
    label = []
    for dim in range(0, 3):
        ax1.plot(xv[dim], yv[dim], f'-{color[dim]}', linewidth=0.5)
        x, y = np.real(f(kval)[dim]), np.imag(f(kval)[dim])
        ax1.plot(x, y, f'{color[dim]}*')
        label.append(f'{x:1.2f}{y:+1.2f}j')
    ax1.set_xlabel('Raíces: '+',  '.join(label))
    if old_ax != (0.0, 1.0, 0.0, 1.0):
        ax1.axis(old_ax)

ax1 = plt.axes([0.075, 0.15, 0.8, 0.8])
ax2 = plt.axes([0.05, 0.04, 0.85, 0.015])
w = Slider(ax2, '$k$', valmin=kmin, valmax=kmax, valinit=kmin, valfmt="%1.6f")
rloc(kmin)
w.on_changed(rloc)
