import numpy as np
import matplotlib.pyplot as plt

from binomialTree import *
from blackScholes import *

price = BS(0,10,10,2.5,0.02,0.3,'p')
print 'Black Scholes put price',price

price = BS(0,10,10,2.5,0.02,0.3,'c')
print 'Black Scholes call price',price

price = trigeorgisTree(10,2.5,10,0.3,0.02,1000,'p')
print 'Trigeorgis put price',price

price = trigeorgisTree(10,2.5,10,0.3,0.02,1000,'c')
print 'Trigeorgis call price',price

price = CRRTree(10,2.5,10,0.3,0.02,1000,'p')
print 'Cox-Ross-Rubinstein put price',price

price = CRRTree(10,2.5,10,0.3,0.02,1000,'c')
print 'Cox-Ross-Rubinstein call price',price

price = JRTree(10,2.5,10,0.3,0.02,1000,'p')
print 'Jarrow-Rudd put price',price

price = JRTree(10,2.5,10,0.3,0.02,1000,'c')
print 'Jarrow-Rudd call price',price

runs = list(range(500,10500,500))
trigeorgis = []
CRR = []
JR = []

BSPrice = BS(0,10,10,2.5,0.02,0.3,'c')
BS = [BSPrice,BSPrice]
BSRuns = [500,10000]

for i in runs:
    print i
    trigeorgis.append(trigeorgisTree(10,2.5,10,0.3,0.02,i,'c'))
    CRR.append(CRRTree(10,2.5,10,0.3,0.02,i,'c'))
    JR.append(JRTree(10,2.5,10,0.3,0.02,i,'c'))

plt.plot(BSRuns, BS, label='Black Scholes')
plt.plot(runs, trigeorgis, label='Trigeorgis')
plt.plot(runs, CRR, label='CRR')
plt.plot(runs, JR, label='JR')
plt.legend(loc='upper right')
plt.show()