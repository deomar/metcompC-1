import numpy as np
from numpy.linalg import *
        
def init():
    N=25
    L=int(np.sqrt(N))
    dt=0.05
    tmax=200
    v0 = 1
    phi = 2*np.pi*np.random.random(N)
    x = []
    y= []
    for j in range(L):
        for i in range(L):
            x.append(0.5+i)
            y.append(0.5+j)
    vx = v0*np.cos(phi)
    vy = v0*np.sin(phi)
    theta = 2*np.pi*np.random.random(N)
    itmax=int(tmax/dt)
    tau = 1
    Req = 5/6
    R0 = 1
    Fad = 0.75
    Frep = 30
    mi = 1
    nse_lim = 0.6/(2*np.sqrt(dt))

    return L,N,dt,tmax,x,y,vx,vy,itmax,tau,v0,R0, Req, Fad,Frep,mi,theta


class particle():
    def __init__(self,x,y,vx,vy,ident,theta,N):
        self.r = np.array([x,y])
        self.v = np.array([vx,vy])
        self.T = theta
        self.ni = np.array([np.cos(self.T),np.sin(self.T)])
        self.ident = ident
        self.Force = np.zeros(2)
        self.u = self.v/np.linalg.norm(self.v)
    def evol_r(self,dt):
        self.r=self.r+self.v*dt
        return

    def evol_v(self):
        self.v = v0*self.ni + mi*self.Force
        self.u = self.v/np.linalg.norm(self.v)
        return
        
    def evol_theta(self,dt):
        xi = np.random.uniform(-nse_lim,nse_lim)
        vnorm = np.linalg.norm(self.v)
        p = self.v/vnorm
        prod = (self.ni[0]*p[1] - self.ni[1]*p[0])
        self.T += xi*np.sqrt(dt) + np.arcsin(prod)*dt
        return
        
    def evol_ni(self):
        self.ni[0] = np.cos(self.T)
        self.ni[1] = np.sin(self.T)
        
    def reflex(self,L):
        if self.r[0] >= L : self.r[0]= self.r[0]-L
        if self.r[0] <= 0 : self.r[0]= self.r[0]+L
        if self.r[1] >= L : self.r[1]= self.r[1]-L
        if self.r[1] <= 0 : self.r[1]= self.r[1]+L
        
    def force(self,dr):
        d=np.linalg.norm(dr)
        eij = dr/d
        if d < Req :
            f = Frep*((d-Req)/Req)
        elif d <= R0:
            f = Fad*((d-Req)/(R0-Req))
        else :
            f= 0
        return eij*f
        
    def force(self,dr):
        d=np.linalg.norm(dr)
        eij = dr/d
        if d < Req :
            f = Frep*((d-Req)/Req)
        elif d <= R0:
            f = Fad*((d-Req)/(R0-Req))
        else :
            f= 0
        return eij*f
        
    def force(self,dr):
        d=np.linalg.norm(dr)
        eij = dr/d
        if d < Req :
            f = Frep*((d-Req)/Req)
        elif d <= R0:
            f = Fad*((d-Req)/(R0-Req))
        else :
            f= 0
        return eij*f
        
    def forces_between_particles(self,N):
        for i in range(self.ident+1,N):
            dr = part[i].r-self.r
            forc=self.force(dr)
            self.Force+=forc
            part[i].Force-=forc     
        return

    def zero_forces(self,N): 
        self.Force=0.0
        return
          
L,N,dt,tmax,x,y,vx,vy,itmax,tau,v0,R0, Req, Fad,Frep,mi,theta = init()
part = list(particle(x[i],y[i],vx[i],vy[i],i,theta[i],N) for i in range(N))
t=0

for it in range(itmax):
    t += dt
    vx=[]
    vy=[]
    list(map(lambda i:i.zero_forces(N), part))
    list(map(lambda i:i.forces_between_particles(N), part))
    list(map(lambda i:i.evol_v(), part))
    list(map(lambda i:i.evol_r(dt), part))
    list(map(lambda i:i.reflex(L), part))
    list(map(lambda i:i.evol_theta(dt), part))
    list(map(lambda i:i.evol_ni(), part))
    list(map(lambda i:i.zero_forces(N), part))
        
    # Temos dois vetores onde são as coordenadas X e Y das velocidades dividas pelo seus módulos 
    list(map(lambda i:vx.append(i.u[0]), part))
    list(map(lambda i:vy.append(i.u[1]), part))
        
    #Printamos o valor de V por t
    print(t,np.sqrt(sum(vx)**2 + sum(vy)**2)/N)
        
        
