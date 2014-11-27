#from matplotlib import pyplot as plt
#import matplotlib.animation as animation
import  numpy as np
import math
from scipy import integrate
import time
import plot

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
                u[-1,1] = u[-2,1]+dx*P["r"](tj)
            else:
                for i in range(1,Nx-1):
                    u[i, j] = (S2(x[i]))*u[i+1,j-1] + (2-2*(S2(x[i])))*u[i,j-1] + (S2(x[i]))*u[i-1,j-1] - u[i,j-2]
                u[ 0,j] = P["l"](tj)
                u[-1,j] = u[-2,j]+dx*P["r"](tj)
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
                u[ 0,1] = u[1,1]+dx*P["l"](tj)
                u[-1,1] = P["r"](tj)
            else:
                for i in range(1,Nx-1):
                    u[i, j] = (S2(x[i]))*u[i+1,j-1] + (2-2*(S2(x[i])))*u[i,j-1] + (S2(x[i]))*u[i-1,j-1] - u[i,j-2]
                u[ 0,j] = u[1,j]+dx*P["l"](tj)
                u[-1,j] = P["r"](tj)
        return x, t, u
    else:#robin
        for j,tj in enumerate(t):
            if(j==0):
                for i in range(1,Nx-1):
                    u[i,0] = P["f"](x[i])
                u[ 0,0] = (float(dx)/(alpha*dx+alpha-1))*P["l"](tj) - u[1,0]*((1-alpha)/(dx*alpha+alpha-1))#Usando Foward Difference
                u[-1,0] = P["r"](tj)*(dx/(beta*dx-beta+1))+ u[-2,0]*((1-beta)/(beta*dx-beta+1))
            elif(j==1):
                for i in range(1,Nx-1):
                    u[i, j] = .5*(S2(x[i]))*u[i+1,j-1] + (1-(S2(x[i])))*u[i,j-1] + .5*(S2(x[i]))*u[i-1,j-1] - dt*P["g"](x[i])
                u[ 0,1] = (float(dx)/(alpha*dx+alpha-1))*P["l"](tj) - u[1,1]*((1-alpha)/(dx*alpha+alpha-1))
                u[-1,1] = P["r"](tj)*(dx/(beta*dx-beta+1))+ u[-2,1]*((1-beta)/(beta*dx-beta+1))
            else:
                for i in range(1,Nx-1):
                    u[i, j] = (S2(x[i]))*u[i+1,j-1] + (2-2*(S2(x[i])))*u[i,j-1] + (S2(x[i]))*u[i-1,j-1] - u[i,j-2]
                u[ 0,j] = (float(dx)/(alpha*dx+alpha-1))*P["l"](tj) - u[1,j]*((1-alpha)/(dx*alpha+alpha-1))
                u[-1,j] = P["r"](tj)*(dx/(beta*dx-beta+1))+ u[-2,j]*((1-beta)/(beta*dx-beta+1))
        return x, t, u

def seccion1():
    f1 = lambda x: np.exp(-200*x**2)
    g1 = lambda x: 0#400*x*np.exp(-200*x**2)
    l1 = lambda t: t*0
    r1 = lambda t: t*0
    c1=lambda x: 1
    t_max=100
    alpha=1
    beta=1
    edptype=1
    soltype=1
    h=0.002
    k=0.001
    P1 = {"xmin":-1, "xmax":1, "tmin":0, "tmax":t_max  , "c":c1,"f":f1,"g":g1, "l":l1, "r":r1}
    x,t,u=diferenciacion(P1,h,k,edptype,soltype,alpha,beta)
    plot.show(x,t,u)
def seccion2():
    f1 = lambda x: np.exp(-200*x**2)
    g1 = lambda x: 400*x*np.exp(-200*x**2)
    l1 = lambda t: t*0
    r1 = lambda t: t*0
    c1=lambda x: 1
    t_max=1
    alpha=1
    beta=1
    edptype=1
    soltype=1
    h1=0.02
    k1=0.001
    h2=0.004
    k2=0.001
    h3=0.002
    k3=0.001

    P1 = {"xmin":-1, "xmax":1, "tmin":0, "tmax":t_max  , "c":c1,"f":f1,"g":g1, "l":l1, "r":r1}
    tiempo1 = time.clock()
    x1,t1,u1=diferenciacion(P1,h1,k1,edptype,soltype,alpha,beta)
    tiempo2 = time.clock()
    tiempo100=tiempo2-tiempo1

    tiempo3 = time.clock()
    x2,t2,u2=diferenciacion(P1,h2,k2,edptype,soltype,alpha,beta)
    tiempo4 = time.clock()
    tiempo500=tiempo4-tiempo3

    tiempo5 = time.clock()
    x3,t3,u3=diferenciacion(P1,h3,k3,edptype,soltype,alpha,beta)
    tiempo6 = time.clock()
    tiempo1000=tiempo6-tiempo5

    print("el tiempo de 100 es:",tiempo100)
    print("el tiempo de 500 es:",tiempo500)
    print("el tiempo de 1000 es:",tiempo1000)
    plot.show(x1,t1,u1)
    plot.show(x2,t2,u2)
    plot.show(x3,t3,u3)

def seccion3():
    f1 = lambda x: np.exp(-200*x**2)
    g1 = lambda x: 400*x*np.exp(-200*x**2)
    l1 = lambda t: t*0
    r1 = lambda t: t*0
    c1=lambda x: 1
    t_max=1
    alpha=1
    beta=1
    edptype=1
    soltype=1
    h1=0.004
    k1=0.0009
    h2=0.004
    k2=0.001
    h3=0.004
    k3=0.01

    P1 = {"xmin":-1, "xmax":1, "tmin":0, "tmax":t_max  , "c":c1,"f":f1,"g":g1, "l":l1, "r":r1}
    x1,t1,u1=diferenciacion(P1,h1,k1,edptype,soltype,alpha,beta)
    x2,t2,u2=diferenciacion(P1,h2,k2,edptype,soltype,alpha,beta)
    x3,t3,u3=diferenciacion(P1,h3,k3,edptype,soltype,alpha,beta)
    plot.show(x1,t1,u1)#funciona
    plot.show(x2,t2,u2)#funciona
    plot.show(x3,t3,u3)#muere

def seccion4():
    f1 = lambda x: x*0
    g1 = lambda x: math.sin(math.pi*(x+1))
    l1 = lambda t: t*0
    r1 = lambda t: t*0
    c1=lambda x: 2
    t_max=1
    alpha1=1
    alpha0=0
    beta1=1
    beta0=0
    edptype=1
    soltype=1
    h=0.004
    k=0.001
    P1 = {"xmin":-1, "xmax":1, "tmin":0, "tmax":t_max  , "c":c1,"f":f1,"g":g1, "l":l1, "r":r1}
    x1,t1,u1=diferenciacion(P1,h,k,edptype,soltype,alpha1,beta1)#dir
    x2,t2,u2=diferenciacion(P1,h,k,edptype,soltype,alpha0,beta0)#neu
    x3,t3,u3=diferenciacion(P1,h,k,edptype,soltype,alpha1,beta0)#mix1
    x4,t4,u4=diferenciacion(P1,h,k,edptype,soltype,alpha0,beta1)#mix2
    x5,t5,u5=diferenciacion(P1,h,k,edptype,soltype,0.5,0.5)#robin
    plot.show(x1,t1,u1)
    plot.show(x2,t2,u2)
    plot.show(x3,t3,u3)
    plot.show(x4,t4,u4)
    plot.show(x5,t5,u5)

def seccion5():
    f1 = lambda x: np.exp(-200*x**2)+((x+1)*x)/2.0
    g1 = lambda x: (1.0/2.0)+x-400*np.exp(-200*x**2)*x
    l1 = lambda t: math.sin(t)
    r1 = lambda t: ((math.sin(t))*(math.cos(t)))/t
    c1=lambda x: (1.0/5.0)+((math.sin(x-1))**2)
    t_max=20
    alpha=1
    beta=1
    edptype=1
    soltype=1
    k=0.001
    h=0.004
    P1 = {"xmin":-1, "xmax":1, "tmin":0, "tmax":t_max  , "c":c1,"f":f1,"g":g1, "l":l1, "r":r1}
    x,t,u=diferenciacion(P1,h,k,edptype,soltype,alpha,beta)
    plot.show(x,t,u)

if __name__ == '__main__':
    seccion5()
