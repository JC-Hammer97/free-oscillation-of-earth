# -*- coding: utf-8 -*-

import numpy
import pylab



def sigma(x): 
    p=0                    
    r = x*(6371.0)    
    if 6371.0 > r >= 6368.0:                    #creates values of rho, p, for each layer within the earth 
        p+=1.030                                #which can be used to create values of sigma
    elif r >= 6357.0:
        p+=2.802
    elif r >= 6352.0:
        p+=2.902
    elif r >= 5951.0:
        p+=7.15855 - 3.85999*(x)
    elif r >= 5701.0:
        p+=11.11978 - 7.87054*(x)
    elif r >= 3485.7:
        p+=6.81430 - 1.66273*(x) - 1.18531*(x**2)
    elif r >= 1217.1:
        p+=12.58416 - 1.69929*(x) - 1.94128*(x**2) - 7.11215*(x**3)
    elif r >= 0:
        p+=13.01219 - 8.45292*(x**2)

    p=p*(10**3)
    p0=5.55e3  
    sig=p/p0
    
    return p, sig
"""
j = []                                          #create empty list
for i in numpy.arange(0,1,0.001):               
    j.append(sigma(i)[1])                       #adds sigma to list for each increment value
pylab.plot(numpy.arange(0,1,0.001)*6371,j)      #sigma(i)[1] is the density, sigma(i)[0] is the value of sigma
pylab.xlabel('Radius (km)')
pylab.ylabel('Sigma')                                                #plots sigma against x, x is the radius but can be scaled 0-1

"""

def m(x,p):
    Vs=0
    r = x*(6371.0)
    if 6371.0 > r >= 6368.0:
        Vs+=0   
    elif r >= 6357.0:
        Vs+=3.55000
    elif r >= 6352.0:
        Vs+=3.75000
    elif r >= 6291.0:
        Vs+=4.65400
    elif r >= 6151.0:
        Vs+=4.34060
    elif r >= 5951.0:
        Vs+=15.09536 - 11.01544*(x)
    elif r >= 5701.0:
        Vs+=15.04371 - 10.69726*(x)
    elif r >= 3485.7:
        Vs+=9.20501 - 6.85512*(x) + 9.39892*(x**2) - 6.25575*(x**3)
    elif r >= 1217.1:
        Vs+=0
    elif r >= 0:
        Vs+=3.56454 - 3.45241*(x**2)
        
    Vs=Vs*0.001
    mu=(Vs**2)*p*1e12
    mu0=2.911e11
    em=mu/mu0

    return em

"""
q = []                                          
for w in numpy.arange(0,1,0.001):               
    q.append(m(w,sigma(w)[0]))                
pylab.plot(numpy.arange(0,1,0.001)*6371,q)      
pylab.xlabel('Radius (km)')
pylab.ylabel('m')


"""                                        

#print(m(0.58076,sigma(0.58076)[0]))
 

"""
def runkut(n,x,y,h):
    "Advances the solution of diff eqn defined by derivs from x to x+h"
    y0=y[:]
    k1=dy(n,x,y)
    for i in range(n):
        y[i]=y0[i]+0.5*h*k1[i]
    k2=dy(n,x+0.5*h,y)
    for i in range(n):
        y[i]=y0[i]+h*(0.2071067811*k1[i]+0.2928932188*k2[i])
    k3=dy(n,x+0.5*h,y)
    for i in range(n):
        y[i]=y0[i]-h*(0.7071067811* k2[i]-1.7071067811*k3[i])
    k4=dy(n,x+h,y)
    for i in range(n):
        a=k1[i]+0.5857864376*k2[i]+3.4142135623*k3[i]+k4[i]
        y[i]=y0[i]+0.16666666667*h*a
        
    x+=h
    return(x,y)
    

def dy(n,x,y):
    "Defines the differential eqns"
    s=sigma(x)[1]
    M=m(x,sigma(x)[0])
    #M=s=1
    node=3
    y1 = y[0]
    y2 = y[1]
    u = [0,0]
    u[0] = y1/x + y2/M
    u[1] = ((M*((node**2)+node-2)/(x**2))-((s)*h))*y1 - 3*y2/x
    return u

y2list=[]

for h in numpy.arange(5.0,20.0,0.001):
    N=452
    x=0.548; y=[1.0,0.0] 
    #print(h)
    for j in range(0,N):
        (x,y)=runkut(2,x,y,0.001)
        #print(x,y[0],y[1])
    #print(x,y[1])
    y2list.append(y[1])
    
    if -0.0001 < y2list[-1] < 0.0001:
        print (h)

pylab.plot(numpy.arange(5.0, 20.0, 0.001), y2list)
pylab.xlabel('lambda')
pylab.ylabel('y2 at Earths surface (m)')


"""

#5.8   4.406
#14.3  10.523
#25.3  17.948 with variable sigma

p0=5.55e3 
mu0=2.911e11
h1=4.406
omega=((h1*mu0)/(p0*((6371*1000)**2)))**0.5
#print(omega)
period=(2*numpy.pi)/omega
#print(period)
mins=period/60
print(mins)
h2=10.523
omega=((h2*mu0)/(p0*((6371*1000)**2)))**0.5
#print(omega)
period=(2*numpy.pi)/omega
#print(period)
mins=period/60
print(mins)
h3=17.948
omega=((h3*mu0)/(p0*((6371*1000)**2)))**0.5
#print(omega)
period=(2*numpy.pi)/omega
#print(period)
mins=period/60
print(mins)

