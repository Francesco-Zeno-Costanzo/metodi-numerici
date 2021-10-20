import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt

N=10000 	#grandezza del campionamento
M=(2**31)-1 	#periodo del generatore

#funzione che genera numeri casuali fra 0 ed M
def GEN(r0,n=1, a=16807,c=0): 
    if n==1:
        return np.mod(a*r0+c,M)
    else:
        r=np.zeros(n)
        r[0]=r0
        for i in range(1,n):
            r[i]=np.mod(a*r[i-1]+c,M)
    return r 

#funzione da campionare  
def f(x, g):
	return np.sqrt(1-x**2)*np.exp(-g*x)

#funzione che sappiamo campionare analiticamente
def h(x, g):
	return np.exp(-g*x)

#genero il campione uniforme e lo normalizzo
l=-GEN(97, N)
l1=l/M
u=np.array([])#array che conterrà l'accettanza al variare di g
q=np.arange(0.1, 11, 0.1)

x=np.linspace(-1., 1., 1000000)#array per i plot

for g in q:
	#campionamento analitico		
	y=np.log((1/np.sinh(g))/(2*l1+np.exp(g)/np.sinh(g)))/g
	#campionamento di f
	naccept=0
	counter=0
	t=np.array([])#array che conterrà il campionamento
	while naccept<N:
		for i in range(len(y)):
		    	counter+=1
		    	#numero casuale in [0, c*h(x)] nella fattispecie va bene c=1
		    	z=((GEN(-l[i]))/M)*h(y[i], g)
		    	if z<f(y[i], g): #test di accettanza
		    		t=np.insert(t, len(t), y[i])
		    		naccept+=1
	#grafici del campionamento vs distribuzione per alcuni g
	#sia campionamento che distribuzione sono normalizzati
	if g==1.0 or g==2.0 or g==4.0 or g==6.0:
		plt.figure(1)
		plt.suptitle("$\sqrt{1-x^2}e^{-\gamma x}$")
		if g==1.0:
			Norm=si.simps(f(x, g), x, dx=1/len(x), even='first')
			plt.subplot(221)
			plt.title('$\gamma=1$')
			plt.hist(t, 100,  histtype = 'step', density=True)
			plt.plot(x, f(x, g)/Norm)
			plt.grid()
		if g==2.0:
			Norm=si.simps(f(x, g), x, dx=1/len(x), even='first' )
			plt.subplot(222)
			plt.title('$\gamma=2$')
			plt.hist(t, 100,  histtype = 'step', density=True)
			plt.plot(x, f(x, g)/Norm)	
			plt.grid()
		if g==4.0:
			Norm=si.simps(f(x, g), x, dx=1/len(x), even='first' )
			plt.subplot(223)
			plt.title('$\gamma=4$')
			plt.hist(t, 100,  histtype = 'step', density=True)
			plt.plot(x, f(x, g)/Norm)	
			plt.grid()	
		if g==6.0:
			Norm=si.simps(f(x, g), x, dx=1/len(x), even='first' )
			plt.subplot(224)
			plt.title('$\gamma=6$')
			plt.hist(t, 100,  histtype = 'step', density=True)
			plt.plot(x, f(x, g)/Norm)	
			plt.grid()
	#calcolo accettanza		
	u=np.insert(u, len(u), naccept/counter)

#grafico l'accettanza
plt.figure(2)
plt.title("Accettanza in funzione di $\gamma$")
plt.plot(q, u)
plt.xscale('log')
plt.grid()
plt.show()
