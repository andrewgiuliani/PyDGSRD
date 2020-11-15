import numpy as np
import ipdb
def LegendreVandermonde(x, p):
    V = np.zeros( (x.size, p+1) )
    dV = np.zeros( (x.size, p+1) )
   
#    ipdb.set_trace(context=21)
    odd = 2*np.arange(11) + 1
    normalization = np.sqrt(odd)/odd

    #normalization = np.array([np.sqrt(2.), (1./3.)*np.sqrt(6.), (1./5.)*np.sqrt(10.), (1./7.)*np.sqrt(14.), (1./3.)*np.sqrt(2.), (1./11.)*np.sqrt(22.), (1./13.)*np.sqrt(26.), (1./15.)*np.sqrt(30.), (1./17.)*np.sqrt(34.), (1./19.)*np.sqrt(38.), (1./21.)*np.sqrt(42.)])
    for i in range(p+1):
        V[:,i] = Legendre(x,i)/normalization[i]
        dV[:,i] = 2*dLegendre(x,i)/normalization[i]
        #dV[:,i] = dLegendre(x,i)/normalization[i]
    return V,dV



def Legendre(x,p) : 
    if p==0 :
        y = np.ones(x.size)
    elif p==1:
        y = x
    elif p==2:
        y = (3*x**2-1)/2
    elif p==3:
        y = (5*x**3 - 3*x)/2
    elif p==4:
        y = (35*x**4 - 30*x**2 + 3)/8
    elif p==5:
        y = (63*x**5 - 70*x**3 + 15*x)/8
    elif p==6:
        y = (231*x**6 - 315*x**4 + 105*x**2 -5)/16
    elif p==7:
        y = (429*x**7 - 693*x**5 + 315*x**3 - 35*x)/16
    elif p==8:
        y = (6435*x**8 - 12012*x**6 + 6930*x**4 - 1260*x**2 + 35)/128
    elif p==9:
        y = (12155*x**9 - 25740*x**7 + 18018*x**5 - 4620*x**3 + 315*x)/128
    elif p==10:
        y = (46189*x**10 - 109395*x**8 + 90090*x**6 - 30030*x**4 + 3465*x**2 - 63)/256
    return y

def dLegendre(x,p) :
    if p==0:
        y = np.zeros(x.size)
    elif p==1:
        y = np.ones(x.size)
    elif p==2:
        y = 3*x
    elif p==3:
        y = (15*x**2 - 3)/2
    elif p==4:
        y = (140*x**3 - 60*x)/8
    elif p==5:
        y = (315*x**4 - 210*x**2 + 15)/8
    elif p==6:
        y = (1386*x**5 - 1260*x**3 + 210*x)/16
    elif p==7:
        y = (3003*x**6 - 3465*x**4 + 945*x**2 - 35)/16
    elif p==8:
        y = (51480*x**7 - 72072*x**5 + 27720*x**3 - 2520*x)/128
    elif p==9:
        y = (109395*x**8 - 180180*x**6 + 90090*x**4 - 13860*x**2 + 315)/128
    elif p==10:
        y = (461890*x**9 - 875160*x**7 + 540540*x**5 - 120120*x**3 + 6930*x)/256
    return y
