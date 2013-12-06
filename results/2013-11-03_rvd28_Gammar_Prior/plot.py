#!/usr/bin/env python


import scipy.stats as ss
import numpy as np
import matplotlib.pyplot as plt

alpha, loc, beta=2.5e4, 1, 5
data=ss.gamma.rvs(alpha,loc=loc,scale=beta,size=5000)
myHist = plt.hist(data, 100, normed=True)
rv = ss.gamma(alpha,loc,beta)
x = np.linspace(0,100) 
h = plt.plot(x, rv.pdf(x), lw=2)
plt.show()