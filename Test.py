'''
Created on 25/05/2014

@author: joseprubiol
'''

from binomialTree import *
from blackScholes import *

preu = BS(0,10,10,2.5,0.02,0.3,'p')
print 'Black Scholes put price',preu

preu = BS(0,10,10,2.5,0.02,0.3,'c')
print 'Black Scholes call price',preu

preu = trigeorgisTree(10,2.5,10,0.3,0.02,1000,'p')
print 'Trigeorgis put price',preu

preu = trigeorgisTree(10,2.5,10,0.3,0.02,1000,'c')
print 'Trigeorgis call price',preu

preu = CRRTree(10,2.5,10,0.3,0.02,1000,'p')
print 'Cox-Ross-Rubinstein put price',preu

preu = CRRTree(10,2.5,10,0.3,0.02,1000,'c')
print 'Cox-Ross-Rubinstein call price',preu

preu = JRTree(10,2.5,10,0.3,0.02,1000,'p')
print 'Jarrow-Rudd put price',preu

preu = JRTree(10,2.5,10,0.3,0.02,1000,'c')
print 'Jarrow-Rudd call price',preu