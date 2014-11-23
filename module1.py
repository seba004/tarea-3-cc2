import numpy as np
#import plot

def wave_equation_forward_differences(P, Nx, Nt):
    """
    Solves the heat equation using forward differences
    """
    x = np.linspace(P["xmin"], P["xmax"], Nx)
    t = np.linspace(P["tmin"], P["tmax"], Nt)
    dx = x[1]-x[0]
    dt = t[1]-t[0]
    S2 = (P["c"]*dt/dx)**2
    print "CFL condition: (c*dt/dx)^2 = %.4f <= 1.0 ?" %S2
    # Storage
    u = np.zeros((Nx, Nt))
    # Time Loop
    i = np.arange(1,Nx-1)
    for j, tj in enumerate(t):
        if j==0:
            u[:,0] = P["f"](x)
            u[ 0,0] = P["l"](tj)
            u[-1,0] = P["r"](tj)
        elif j==1:
            u[i, j] = .5*S2*u[i+1,j-1] + (1-S2)*u[i,j-1] + .5*S2*u[i-1,j-1] - dt*P["g"](x[i])
            u[ 0,1] = P["l"](tj)
            u[-1,1] = P["r"](tj)
        else:
            u[i, j] = S2*u[i+1,j-1] + (2-2*S2)*u[i,j-1] + S2*u[i-1,j-1] - u[i,j-2]
            u[ 0,j] = P["l"](tj)
            u[-1,j] = P["r"](tj)
    return x, t, u





shift = lambda x : np.where(x<=0.4, 0., 1.)*np.where(x>=0.6, 0., 1.)

f1 = lambda x: np.sin(5*2*np.pi*x) * shift(x)
g1 = lambda x: 0
l1 = lambda t: 0*t
r1 = lambda t: 0*t
P1 = {"xmin":0, "xmax":1, "tmin":0, "tmax":1, "c":1.0,
      "f":f1, "g":g1, "l":l1, "r":r1}

f2 = lambda x: np.sin(5*2*np.pi*x)
g2 = lambda x: np.cos(5*2*np.pi*x)
l2 = lambda t: 0*t
r2 = lambda t: 0*t
P2 = {"xmin":0, "xmax":1, "tmin":0, "tmax":1, "c":1.0,
      "f":f2, "g":g2, "l":l2, "r":r2}

f3 = lambda x: np.sin(2*np.pi*x)
g3 = lambda x: 2*np.pi*np.cos(2*np.pi*x)
l3 = lambda t: 0*t
r3 = lambda t: 0*t
P3 = {"xmin":0, "xmax":1, "tmin":0, "tmax":1, "c":1.,
      "f":f3, "g":g3, "l":l3, "r":r3}

f4 = lambda x: np.exp(-(x-.5)**2/0.01)
g4 = lambda x: 0*x
l4 = lambda t: t
r4 = lambda t: t
P4 = {"xmin":0, "xmax":1, "tmin":0, "tmax":1, "c":1.,
      "f":f4, "g":g4, "l":l4, "r":r4}





P = P4
#x, t, u = wave_equation_forward_differences(P, 40, 20) # Unstable
#x, t, u = wave_equation_forward_differences(P, 40, 40) # Stable
x, t, u = wave_equation_forward_differences(P, 100, 100) # Stable

#plot.show(x,t,u)
