from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

plt.ion()

x0_e = 8.48528137423857

def LorenzControl(σ=10, β=8/3, ρ=28, k=1):
    def f(x, t=0):
        # x_1(t) ----> x[0], x_2(t) ----> x[1], x_3(t) ----> x[2]
        # x_1(t+dt) -> y[0], x_2(t+dt) -> y[1], x_3(t+dt) -> y[2]
        y = np.empty_like(x)
        y[0] = σ*(x[1] - x[0])
        y[1] = x[0]*(ρ - x[2]) - x[1] - k*(x[0]-x0_e)
        y[2] = x[0]*x[1] - β*x[2]
        return y
    return f


# Valores máximos y mínimos de k
kmin = 0
kmax = 1.25
kvec = np.linspace(kmin, kmax, 100)

t = np.linspace(0.0, 40.0, 2000) # array de tiempo de integración
x0 = np.array([1.0, 1.0, 1.0]) # Condiciones iniciales
def LorenzStability(kval):
    x = odeint(LorenzControl(k=kval), x0, t).T
    ax1.cla()
    ax1.plot(t, x[0], 'b-')
    ax1.plot(t, x[1], 'r-')
    ax1.plot(t, x[2], 'g-')
    ax1.legend(('x', 'y', 'z'), loc='upper right')
    
ax1 = plt.axes([0.075, 0.15, 0.8, 0.8])
ax2 = plt.axes([0.05, 0.06, 0.85, 0.015])
w = Slider(ax2, '$k$', valmin=kmin, valmax=kmax, valinit=kmin, valfmt="%1.6f")
LorenzStability(kmin)
ax2.xaxis.set_visible(True)
ax2.set_xticks(np.concatenate((np.linspace(kmin, kmax, 5), [0.50101])))
w.on_changed(LorenzStability)
