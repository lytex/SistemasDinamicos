# Visualización de cómo afecta el control al espacio de fases
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import numpy as np
from random import choice

plt.ion()

x0_e = 8.48528137423857

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
kmax = 1.25
kvec = np.linspace(kmin, kmax, 100)

t = np.linspace(0.0, 40.0, 4000) # array de tiempo de integración
x0 = np.array([1.0, 1.0, 1.0]) # Condiciones iniciales

plt.ion()

# Variables globales: animación de la bolita
ani = None

def animate(i, *fargs):
    x = fargs[0]
    ln = fargs[1]
    # Pintar la bolita en el instante t[i]
    ln.set_data_3d([x[0][i]], [x[1][i]], [x[2][i]])
    ln.set_color('r')
    ln.set_marker('.')
    ln.set_markersize(10)
    return ln,

def LorenzButterfly(kval):
    global ani
    if ani is not None:
        ani.event_source.stop() # parar ani si no está corriendo
    x = odeint(LorenzControl(k=kval), x0, t).T # integrar el modelo
    ax1.cla()
    fix_ln = ax1.plot(x[0], x[1], x[2], 'b-')[0]
    ln = ax1.plot([], [], [])[0]
    fig = plt.gcf()
    ani = FuncAnimation(fig, animate, frames=np.arange(t.size), repeat=False, fargs=(x,ln), interval=0, blit=True)
    return ani

plt.figure(figsize=(10,8))
ax1 = plt.axes([0.075, 0.15, 0.8, 0.8], projection='3d')
ax2 = plt.axes([0.05, 0.06, 0.85, 0.015])
# Resultan ser unos valores muy buenos para el control:
kinit = choice([0.926, 1.1603])
w = Slider(ax2, '$k$', valmin=kmin, valmax=kmax, valinit=kinit, valfmt="%1.6f")
ax2.xaxis.set_visible(True)
ax2.set_xticks(np.concatenate((np.linspace(kmin, kmax, 5), [0.50101])))

ani = LorenzButterfly(kinit)

plt.connect('button_press_event', lambda event: ani.event_source.stop())
plt.connect('button_release_event', lambda event: ani.event_source.start())
w.on_changed(LorenzButterfly)
