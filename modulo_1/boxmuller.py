import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt
#campionamento di una gaussiana
M=(2**31)-1	#periodo generatore
N=30000		#grandezza del campione
s=1		#varianza della distribuzione da campionare
m=5		#media della distribuzione da campionare

def f(x):
    return np.exp(-((x-m)**2)/(2*s))

#array per il plot e nomalizzazione di f
x=np.linspace(0, 10, 10000)
Norm=si.simps(f(x), x, dx=1/len(x), even='first' )

#funzione che genera numeri casuali fra 0 ed M
def GEN(r0, n=1, a=16807, c=0):
    """
    generatore conguenziale lineare
    Parametri
    ---------
    r0 : int
        seed della generazione
    n : int, opzionale
        dimensione array da generare, di default è 1
    a : int, opzionale
        moltiplicatore del generatore, di default è 16807
    c : int, opzionale
        incremento del generatore, di default è 0

    Returns
    ---------
    r : array
        array con numeri distribuiti casualmente
    """
    if n==1:
        return np.mod(a*r0+c,M)
    else:
        r=np.zeros(n)
        r[0]=r0
        for i in range(1,n):
            r[i]=np.mod(a*r[i-1]+c,M)
    return r

#camipnamento tramite doppio cambio di variabile
x1 = GEN(56, N)/M
x2 = GEN(81, N)/M

phi = 2*np.pi*x1
z = -np.log(1-x2)

y1 = np.sqrt(2*z*s**2)*np.cos(phi) + m
y2 = np.sqrt(2*z*s**2)*np.sin(phi) + m

#media ed errore con estimatore naiv
m1 = np.mean(y1)
dm = np.std(y1)/np.sqrt(len(y1)-1)

#media ed errore tramite resampling, bootstrap
K=100
X=np.zeros((K, N))
X[1,:]=y1
for j in range(K):
    for l in range(N):
        xx=int((GEN(l)/M)*(N))
        X[j, l]=X[1, xx]

h=np.array([])
for t in range(K):
    h=np.insert(h, len(h), np.mean(X[t, :]))

m2 = np.mean(h)
dm2 = np.std(h, ddof=1)

print('media= %f +- %f'%(m1, dm))
print('media= %f +- %f'%(m2, dm2))

#grafico del campionamento, sia istogramma che distribuzione normalizzati
plt.figure(1)
plt.title('Campionamento gaussiana tramite box-muller', fontsize=15)
plt.xlabel('$\mu = %f \pm %f$ naiv \n $\mu = %f \pm %f$ bootstrap' %(m1, dm, m2, dm2), fontsize=15)
plt.hist(y1, 100,  histtype = 'step', density=True)
plt.plot(x, f(x)/Norm)
plt.grid()

plt.show()
