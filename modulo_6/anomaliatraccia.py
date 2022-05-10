import numpy as np
import random as rn
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline

massa = 0.05

nt, o1, o2, o3, do1, do2, do3 = np.loadtxt('misure.dat', unpack=True)
nt, e, de = np.loadtxt('dens.dat', unpack=True)

#densità rinormalizzata
dens = 0.5*nt**2*(e - e[-1])
d_dens = 0.5*nt**2*np.sqrt(de**2 + de[-1]**2)
#temperatura per unità di massa
T = 1/(nt*massa)

#traccia tensore energia impulso
tra = nt**2*(o1 - o1[-1])
dtra = nt**2*np.sqrt(do1**2 + do1[-1]**2)

def errore(x, y, dy):
    '''
    La funzione esegue un boostrap ricampionando la curva
    estrando un nuovo punto con una distibuzione gaussiana
    centrata nel valore centrale del punto e con deviazione
    standard l'errore sul punto; la curva viene interpolata
    e se calcola l'integrale, l'errore sarà deviazione standard
    dei vari integrali calcolati
    '''
    N = 100 # numero ricampionamenti
    P = np.zeros((len(x), N))
    for i in range(N):
        T = x
        m = y
        s = dy
        #estrazione con box-muller
        phi = 2*np.pi*np.random.rand(len(T))
        z = -np.log(1-np.random.rand(len(T)))
        inte = np.sqrt(2*z*s**2)*np.cos(phi) + m
        #interpolazione e calolo dell'integrale
        s2 = InterpolatedUnivariateSpline(T, inte, k=2)
        pr = np.array([s2.integral(T[0], i) for i in T])
        P[:, i] = pr
    
    err = np.std(P, ddof=1, axis=1)
    return err

def pressione(x, y, z, dy):
    '''
    la funzion esegue il calcolo dell'integrale dell'anomalia
    di traccia intermpolando l'integranda con un polinomio di
    grado due e poi calcolando l'integrale vero e proprio
    '''
    T = x[::-1]
    nt = z[::-1]
    tra = y[::-1]
    dtra = dy[::-1]
    #interpolo per calcolare l'integrale
    inte = tra*massa*nt
    dinte = dtra*massa*nt
    
    s2 = InterpolatedUnivariateSpline(T, inte, k=2)
    pr = np.array([s2.integral(T[0], i) for i in T])
    err = errore(T, inte, dinte)
    
    pres = pr[::-1]
    inte = inte[::-1]
    dint = dinte[::-1]
    perr = err[::-1]
    return pres, inte, perr, dint, s2

press, inte, dpress, dinte, spline = pressione(T, tra, nt, dtra)
dens2 = tra + press
ddens2 = np.sqrt(dtra**2 + dpress**2)

t = np.linspace(np.min(T), np.max(T), 1000)
plt.figure(4)
plt.suptitle("Metodo dell'anomalia di traccia")
plt.subplot(211)
plt.title('Traccia e pressione')
plt.errorbar(T, tra, dtra, fmt='.', color='b')
plt.errorbar(T, inte, dinte, fmt='*', color='k')
plt.errorbar(T, press, dpress, fmt='^', color='r')
plt.plot(t, spline(t), 'y')
plt.grid()

plt.subplot(212)
plt.title('Densità di energia normalizata')
plt.ylabel('$\epsilon / T^2$')
plt.xlabel('T/m')
k = np.linspace(np.min(T), np.max(T),1000)
o = np.ones(len(k))
plt.errorbar(T, dens, d_dens, fmt='.', color='b', label='primo metodo')
plt.errorbar(T, dens2, ddens2, fmt='.', color='k', label='metodo anomalia')
plt.plot(k, np.pi/6*o, label='valore asintotico')
plt.legend(loc='best')
plt.grid()


plt.show()
