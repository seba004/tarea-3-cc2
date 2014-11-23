#from matplotlib import pyplot as plt
#import matplotlib.animation as animation
import  numpy as np
import math
from scipy import integrate
import math

#edptype si es 1 es nomal si es 0 es de cond peridocas
#soltype si es 1 es explicita si es 0 es implicita
def diferenciacion(P,h,k,edptype,soltype,alpha,beta):
    if (edptype==1 and soltype ==1):

        x,t,u=normal_edo(P,alpha,beta,h,k)
        return x,t,u

def normal_edo(P,alpha,beta,h,k):
    #se calcula Nx y Nt
    Nx=int((P["xmax"]-P["xmin"])/(float(h)))
    Nt=int((P["tmax"])/(float(k)))
    #se obtien el "mallado"
    x = np.linspace(P["xmin"], P["xmax"], Nx)
    t = np.linspace(P["tmin"], P["tmax"], Nt)
    dx = h
    dt = k
    print("hola")
    S2 = lambda x: ((P["c"](x))*dt/dx)**2
    print (S2)
    #print "CFL condition: (c*dt/dx)^2 = %.4f <= 1.0 ?" %S2
    # Storage
    u = np.zeros((Nx, Nt))
    # Time Loop

    if(alpha==1 and beta==1):
        for j, tj in enumerate(t):
            if j==0:
                for i in range(1,Nx-1):
                    u[i,0] = P["f"](x[i])
                u[ 0,0] = P["l"](tj)
                u[-1,0] = P["r"](tj)
            elif j==1:
                for i in range(1,Nx-1):
                    u[i, j] = .5*(S2(x[i]))*u[i+1,j-1] + (1-(S2(x[i])))*u[i,j-1] + .5*(S2(x[i]))*u[i-1,j-1] - dt*P["g"](x[i])
                u[ 0,1] = P["l"](tj)
                u[-1,1] = P["r"](tj)
            else:
                for i in range(1,Nx-1):
                    u[i, j] = (S2(x[i]))*u[i+1,j-1] + (2-2*(S2(x[i])))*u[i,j-1] + (S2(x[i]))*u[i-1,j-1] - u[i,j-2]
                u[ 0,j] = P["l"](tj)
                u[-1,j] = P["r"](tj)
        return x, t, u

    if(alpha==0 and beta==0):#cond de nueman
        for j,tj in enumerate(t):
            if(j==0):
                for i in range(1,Nx-1):
                    u[i,0] = P["f"](x[i])
                u[ 0,0] = u[1,0]+dx*P["l"](tj)# Foward
                u[-1,0] = u[-2,0]+dx*P["r"](tj)# Backward
            elif(j==1):
                for i in range(1,Nx-1):
                    u[i, j] = .5*(S2(x[i]))*u[i+1,j-1] + (1-(S2(x[i])))*u[i,j-1] + .5*(S2(x[i]))*u[i-1,j-1] - dt*P["g"](x[i])
                u[ 0,1] = u[1,1]+dx*P["l"](tj)
                u[-1,1] = u[-2,1]+dx*P["r"](tj)
            else:
                for i in range(1,Nx-1):
                    u[i, j] = (S2(x[i]))*u[i+1,j-1] + (2-2*(S2(x[i])))*u[i,j-1] + (S2(x[i]))*u[i-1,j-1] - u[i,j-2]
                u[ 0,j] = u[1,j]+dx*P["l"](tj)
                u[-1,j] = u[-2,j]+dx*P["r"](tj)
        return x, t, u
    #mixto con alpha 1  y beta 0
    elif(alpha==1 and beta==0):
        for j,tj in enumerate(t):
            if(j==0):
                for i in range(1,Nx-1):
                    u[i,0] = P["f"](x[i])
                u[ 0,0] = P["l"](tj)
                u[-1,0] = u[-2,0]+dx*P["r"](tj)# Backward
            elif(j==1):
                for i in range(1,Nx-1):
                    u[i, j] = .5*(S2(x[i]))*u[i+1,j-1] + (1-(S2(x[i])))*u[i,j-1] + .5*(S2(x[i]))*u[i-1,j-1] - dt*P["g"](x[i])
                u[ 0,1] = P["l"](tj)
                u[-1,1] = u[-2,1]+dx*P["r"](tj)#Backward
            else:
                for i in range(1,Nx-1):
                    u[i, j] = (S2(x[i]))*u[i+1,j-1] + (2-2*(S2(x[i])))*u[i,j-1] + (S2(x[i]))*u[i-1,j-1] - u[i,j-2]
                u[ 0,j] = P["l"](tj)
                u[-1,j] = u[-2,j]+dx*P["r"](tj)# Backward
        return x, t, u
    #mixto con alpha0 y beta 1
    elif(alpha==0 and beta ==1):
        for j,tj in enumerate(t):
            if(j==0):
                for i in range(1,Nx-1):
                    u[i,0] = P["f"](x[i])
                u[ 0,0] = u[1,0]+dx*P["l"](tj)#Usando Foward Difference
                u[-1,0] = P["r"](tj)
            elif(j==1):
                for i in range(1,Nx-1):
                    u[i, j] = .5*(S2(x[i]))*u[i+1,j-1] + (1-(S2(x[i])))*u[i,j-1] + .5*(S2(x[i]))*u[i-1,j-1] - dt*P["g"](x[i])
                u[ 0,1] = u[1,1]+dx*P["l"](tj)#Usando Foward Difference
                u[-1,1] = P["r"](tj)
            else:
                for i in range(1,Nx-1):
                    u[i, j] = (S2(x[i]))*u[i+1,j-1] + (2-2*(S2(x[i])))*u[i,j-1] + (S2(x[i]))*u[i-1,j-1] - u[i,j-2]
                u[ 0,j] = u[1,j]+dx*P["l"](tj)#Usando Foward Difference
                u[-1,j] = P["r"](tj)
        return x, t, u
    else:#robin
        for j,tj in enumerate(t):
            if(j==0):
                for i in range(1,Nx-1):
                    u[i,0] = P["f"](x[i])
                u[ 0,0] = (float(dx)/(alpha*dx+alpha-1))*P["l"](tj) - u[1,0]*((1-alpha)/(dx*P["alpha"]+P["alpha"]-1))#Usando Foward Difference
                u[-1,0] = P["r"](tj)*(dx/(P["betha"]*dx-P["betha"]+1))+ u[-2,0]*((1-P["betha"])/(P["betha"]*dx-P["betha"]+1))
            elif(j==1):
                for i in range(1,Nx-1):
                    u[i, j] = .5*(S2(x[i]))*u[i+1,j-1] + (1-(S2(x[i])))*u[i,j-1] + .5*(S2(x[i]))*u[i-1,j-1] - dt*P["g"](x[i])
                u[ 0,1] = (float(dx)/(P["alpha"]*dx+P["alpha"]-1))*P["l"](tj) - u[1,1]*((1-P["alpha"])/(dx*P["alpha"]+P["alpha"]-1))
                u[-1,1] = P["r"](tj)*(dx/(P["betha"]*dx-P["betha"]+1))+ u[-2,1]*((1-P["betha"])/(P["betha"]*dx-P["betha"]+1))
            else:
                for i in range(1,Nx-1):
                    u[i, j] = (S2(x[i]))*u[i+1,j-1] + (2-2*(S2(x[i])))*u[i,j-1] + (S2(x[i]))*u[i-1,j-1] - u[i,j-2]
                u[ 0,j] = (float(dx)/(P["alpha"]*dx+P["alpha"]-1))*P["l"](tj) - u[1,j]*((1-P["alpha"])/(dx*P["alpha"]+P["alpha"]-1))
                u[-1,j] = P["r"](tj)*(dx/(P["betha"]*dx-P["betha"]+1))+ u[-2,j]*((1-P["betha"])/(P["betha"]*dx-P["betha"]+1))
        return x, t, u

def seccion1():
    f1 = lambda x: np.exp(-200*x**2)
    g1 = lambda x: 400*x*np.exp(-200*x**2)
    l1 = lambda t: t*0
    r1 = lambda t: t*0
    c1=lambda x: 1
    t_max=100
    alpha=0
    beta=0
    edptype=1
    soltype=1
    h=0.002
    k=0.001
    P1 = {"xmin":0, "xmax":1, "tmin":0, "tmax":t_max  , "c":c1,"f":f1,"g":g1, "l":l1, "r":r1}
    x,t,u=diferenciacion(P1,h,k,edptype,soltype,alpha,beta)
    return x,t,u
    #show(xa,ta,ua)


if __name__ == '__main__':
    x,t,u=seccion1()