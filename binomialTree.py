'''
Created on 22/05/2014

@author: joseprubiol
'''

import math

def trigeorgisTree(K,T,S,sig,r,N,PorC):

    dt=T/N;
    nu=r-0.5*sig*sig;
    dxu=math.sqrt(sig*sig*dt+(nu*dt)**2);
    dxd=-dxu;
    pu=0.5+0.5*(nu*dt/dxu);
    pd=1-pu;
    disc=math.exp(-r*dt);

    St = [0] * (N+1)
    C = [0] * (N+1)
    
    St[0]=S*math.exp(N*dxd);
    
    for j in xrange(1, N+1): 
        St[j] = St[j-1] * math.exp(dxu-dxd);
    
    for j in xrange(1, N+1):
        if PorC == 'p':
            C[j] = max(K-St[j],0);
        elif PorC == 'c':
            C[j] = max(St[j]-K,0);
    
    for i in xrange(N, 0, -1):
        for j in xrange(0, i):
            C[j] = disc*(pu*C[j+1]+pd*C[j]);
            
    return C[0]

def CRRTree(K,T,S,sig,r,N,PorC):
    
    dt=T/N;
    dxu=math.exp(sig*math.sqrt(dt));
    dxd=math.exp(-sig*math.sqrt(dt));
    pu=((math.exp(r*dt))-dxd)/(dxu-dxd);
    pd=1-pu;
    disc=math.exp(-r*dt);

    St = [0] * (N+1)
    C = [0] * (N+1)
    
    St[0]=S*dxd**N;
    
    for j in xrange(1, N+1): 
        St[j] = St[j-1] * math.exp(dxu-dxd);
    
    for j in xrange(1, N+1):
        if PorC == 'p':
            C[j] = max(K-St[j],0);
        elif PorC == 'c':
            C[j] = max(St[j]-K,0);
    
    for i in xrange(N, 0, -1):
        for j in xrange(0, i):
            C[j] = disc*(pu*C[j+1]+pd*C[j]);
            
    return C[0]

def JRTree(K,T,S,sig,r,N,PorC):
    
    dt=T/N;
    dxu=math.exp((r-(sig**2/2))*dt+sig*math.sqrt(dt));
    dxd=math.exp((r-(sig**2/2))*dt-sig*math.sqrt(dt));
    pu=0.5;
    pd=1-pu;
    disc=math.exp(-r*dt);

    St = [0] * (N+1)
    C = [0] * (N+1)
    
    St[0]=S*dxd**N;
    
    for j in xrange(1, N+1): 
        St[j] = St[j-1] * math.exp(dxu-dxd);
    
    for j in xrange(1, N+1):
        if PorC == 'p':
            C[j] = max(K-St[j],0);
        elif PorC == 'c':
            C[j] = max(St[j]-K,0);
    
    for i in xrange(N, 0, -1):
        for j in xrange(0, i):
            C[j] = disc*(pu*C[j+1]+pd*C[j]);
            
    return C[0]
