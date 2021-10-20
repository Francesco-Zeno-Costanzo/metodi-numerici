import numpy as np
import random as rn
import scipy.integrate as si
import matplotlib.pyplot as plt
#campionamento di una gaussiana

s=1		#varianza della distribuzione da campionare
aver=5		#media della distribuzione da campionare

steps, x, acc=np.loadtxt('metrogauss.txt',skiprows=2, unpack=True)
corr=np.loadtxt('metrogausscorr.txt',skiprows=0, unpack=True)
errdb, bloc=np.loadtxt('metrogausserrdb.txt',skiprows=0, unpack=True)
errbt, bloc1=np.loadtxt('metrogausserr.txt',skiprows=0, unpack=True)

N=len(steps)

def p(x):
    return np.exp(-((x-aver)**2)/(2*s))


#array per il plot e nrmalizzazione
t=np.linspace(0, 10, N+1)
Norm=si.simps(p(t), t, dx=1/len(t), even='first' )

print(acc[-1]/N)

#elimino parte del campione per essere sicuro ceh la catena abbia termalizzato
x1=x[10000:]
m1=np.mean(x1)

dm3= np.std(x1)/np.sqrt(len(x1)-1)


fig1 = plt.figure(1)
plt.suptitle('Algoritmo Metropolis, campionamento gaussiana', fontsize=20)

plt.subplot(211)
plt.title('$\overline{x}=%f \pm (\sigma_{blocking} = %f)  (\sigma_{bootstrap} = %f)  (\sigma_{naiv} = %f)$'%(m1, errdb[len(errdb)-2], errbt[len(errbt)-3], dm3), fontsize=20)
plt.grid()
plt.hist(x1, 100,  histtype = 'step', density=True)
plt.plot(t, p(t)/Norm)


#controlliamo la termalizzazione della catena
plt.subplot(212)
plt.title('Evoluzione della catena', fontsize=15)
plt.xlabel('iterazioni', fontsize=15)
plt.plot(steps, x)
plt.grid()

tint=sum(corr)
plt.figure(2)
plt.title(r'autocorrelazione a due punti, $\tau_{int}$=%f' %tint, fontsize=15)
plt.ylabel('C(k)')
plt.xlabel('k')
plt.grid()
plt.plot(corr)

plt.figure(3)
plt.title("errore con datablocking \n i numeri sul grafico sono l'ordine di grandezza del numero  dei blocchi", fontsize=15)
plt.ylabel('errore')
plt.xlabel('k')
x=np.linspace(1, len(errdb), len(errdb))
plt.grid()
plt.yscale('log')
plt.errorbar(x, errdb, fmt='o')
plt.xticks([k for k in range(1,len(errdb)+1)])
for i in range(len(errdb)):
	plt.annotate("%.0e" %int(bloc[i]), (x[i], errdb[i]))

plt.figure(4)
plt.title("errore con bootstrap con blocking \n  i numeri sul grafico sono l'ordine di grandezza del numero  dei blocchi", fontsize=15)
plt.ylabel('errore')
plt.xlabel('k')
x=np.linspace(1, len(errbt), len(errbt))
plt.grid()
plt.yscale('log')
plt.errorbar(x, errbt, fmt='o')
plt.xticks([k for k in range(1,len(errbt)+1)])
for i in range(len(errbt)):
	plt.annotate("%.0e" %int(bloc1[i]), (x[i], errbt[i]))
	
plt.show()
