import numpy as np
from numpy.linalg import *


def init():
    N=100 # Número de células
    L=int(np.sqrt(N)) 
    dt=0.005
    tmax=30
    v0 = 0.5
    # Inicialização das posições e das velocidades das células
    phi = 2*np.pi*np.random.random(N)
    theta = 2*np.pi*np.random.random(N)
    x = []
    y= []
    for j in range(L):
        for i in range(L):
            x.append(0.5+i)
            y.append(0.5+j)
    vx = v0*np.cos(phi)
    vy = v0*np.sin(phi)
    itmax=int(tmax/dt)
    #Valores para os termos conforme o artigo PHYSICAL REVIEW E 74, 061908 (2006)
    tau = 1
    R0 = 1
    Req = 5/6
    Fad = 0.75
    Frep = 30
    mi = 1
    # Intensidade do ruído 
    eta = 5.3
    nse_lim = eta/(2*np.sqrt(dt))
    
    return L,N,dt,tmax,x,y,vx,vy,itmax,tau,v0,R0, Req, Fad,Frep,mi,theta, nse_lim,eta


class particle():
    def __init__(self,x,y,vx,vy,i,theta,N):
        self.r = np.array([x,y])
        self.v = np.array([vx,vy])
        self.T = theta
        self.ni = np.array([np.cos(self.T),np.sin(self.T)])
        self.ident = i
        self.Force = np.zeros(2)
    
    #Concede a distância entre cada partícula para o cálculo da força e atribui para as partículas interagentes. 
    def forces_between_particles(self,N):
        for i in range(self.ident+1,N):
            dr = part[i].r - self.r
            forc=self.force(dr) #self.force(dr) contém as forças atrativas-repulsivas
            self.Force+=forc
            part[i].Force-=forc    
        return
    
    # Cálculo da força dado a distancia entre as partículas
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

    #Evolução da velocidade - Integração do tipo Euler
    def evol_v(self):
        self.v = v0*self.ni + mi*self.Force
        return
        
    #Evolução da posição - Integração do tipo Euler
    def evol_r(self,dt):
        self.r=self.r+self.v*dt
        return
        
    # Condição de contorno periódica
    def Cont_periodico(self,L):  
        if self.r[0] >= L : self.r[0]= self.r[0]-L
        if self.r[0] <= 0 : self.r[0]= self.r[0]+L
        if self.r[1] >= L : self.r[1]= self.r[1]-L
        if self.r[1] <= 0 : self.r[1]= self.r[1]+L        
       
    #Evolução do ângulo theta - Integração do tipo Euler   
    def evol_theta(self,dt):
        xi = np.random.uniform(-nse_lim,nse_lim) #Ruído branco com magnitude nse_lim
        vnorm = np.linalg.norm(self.v)
        u = self.v/vnorm
        prod = (self.ni[0]*u[1] - self.ni[1]*u[0])
        #Como estamos trantando com ruído devemos tirar a raiz de dt para a integração
        self.T += xi*np.sqrt(dt) + np.arcsin(prod)*dt
        return
    
    # Atualização do vetor n da autoproupulsão    
    def evol_ni(self):
        self.ni[0] = np.cos(self.T)
        self.ni[1] = np.sin(self.T)
    
    #Método que zera as forças
    def zero_forces(self,N): 
        self.Force=0.0
        return
    def gnuplot(self):
        print(self.r[0], self.r[1], self.v[0], self.v[1])


L,N,dt,tmax,x,y,vx,vy,itmax,tau,v0,R0, Req, Fad,Frep,mi,theta,nse_lim,eta = init()
part = list(particle(x[i],y[i],vx[i],vy[i],i,theta[i],N) for i in range(N))
t=0

#Script para formatação do gif
print("set terminal gif animate delay 2\n")
print("set output 'video_part.gif'\n")
for it in range(itmax):
    t+=dt
    if it%10==0:
        #Script para formatação do gif
        print("unset ke\n")
        print("set title 'N={}, eta = {} v_0 = {}' font ',14' \n".format(N,eta,v0))
        print("set xrange [0:{}]\n".format(L))
        print("set yrange [0:{}]\n".format(L))
        print("set size square \n")
        print("plot '-' u 1:2:(0.5) w circles fs solid 0.5 lc 2, '-' u 1:2:3:4 w vector filled head lw 2 lc 8 \n")
        list(map(lambda i:i.gnuplot(),part))
        print("e\n")
        
    list(map(lambda i:i.zero_forces(N), part))
    list(map(lambda i:i.forces_between_particles(N), part))
    list(map(lambda i:i.forces_between_particles(N), part))
    list(map(lambda i:i.evol_v(), part))
    list(map(lambda i:i.evol_r(dt), part))
    list(map(lambda i:i.Cont_periodico(L), part))
    list(map(lambda i:i.evol_theta(dt), part))
    list(map(lambda i:i.evol_ni(), part))
    list(map(lambda i:i.zero_forces(N), part))



