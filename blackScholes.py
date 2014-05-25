'''
Created on 22/05/2014

@author: joseprubiol
'''

import math

# Abramowitz and Stegun 7.1.26 approximation
def erf(x):
    
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)

    return y

# custom erf implementation
def fnor1(x): 
    
    y=0.5*(1+erf(x/math.sqrt(2)));
    
    return y

# math erf implementation
def fnor2(x):

    y=0.5*(1+math.erf(x/math.sqrt(2)));
    
    return y

def BS(t,St,K,T,r,sig,PorC):

    Tmt=T-t;
    ATmt=sig*math.sqrt(Tmt);
    logo=math.log(St/K);
    Ap=(logo+(r+0.5*sig**2)*Tmt)/ATmt;
    An=Ap-ATmt;
    
    if PorC == 'c':
        p=St*fnor2(Ap)-K*math.exp(-r*Tmt)*fnor2(An);
    elif PorC == 'p':
        p=K*math.exp(-r*Tmt)*fnor2(-An)-St*fnor2(-Ap);
    
    return p
