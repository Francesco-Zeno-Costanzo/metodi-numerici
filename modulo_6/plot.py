import numpy as np
import random as rn
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline

massa = 0.05

##Densità di energia 
##primo modo di calcolarla

nt, e, de = np.loadtxt('dens.dat', unpack=True)
#densità rinormalizzata
dens = 0.5*nt**2*(e - e[-1])
d_dens = 0.5*nt**2*np.sqrt(de**2 + de[-1]**2)
#temperatura per unità di massa
T = 1/(nt*massa)

plt.figure(1, figsize=(14,6))
plt.subplot(121)
plt.title('Densità di energia non normalizata')
plt.xlabel('nt')
plt.ylabel('$\epsilon/T^2$')
plt.errorbar(nt, e, de, fmt='.')
plt.xscale('log')
plt.yscale('log')
plt.grid()

plt.subplot(122)
plt.title('Densità di energia normalizata')
plt.ylabel('$\epsilon / T^2$')
plt.xlabel('T/m')
k = np.linspace(0, np.max(T),1000)
o = np.ones(len(k))
plt.errorbar(T, dens, d_dens, fmt='.')
plt.plot(k, np.pi/6*o, label='valore asintotico')
plt.legend(loc='best')
plt.grid()

"""
##====================================================================================
"""
## LIMITE AL CONTINUO
## OVER-RELAXATION
nt, eco, deco = np.loadtxt('densm.dat', unpack=True)

dens1 = 0.5*nt**2*eco
d_dens1 = 0.5*nt**2*deco

m = 1/(10*nt)
x = m[2:]
y = dens1[2:]
dy = d_dens1[2:]

def f(x, a, b):
    '''funzione di fit
    '''
    return a + b*x**2
    

#Eseguiamo il fit e stampiamo i risultati:
print('metodo con over-relaxation')
print(f'a_t = {np.pi/6}')
pars, covm = curve_fit(f, x, y, sigma=dy, absolute_sigma=False)
print('a_f = %.5f +- %.5f ' % (pars[0], np.sqrt(covm.diagonal()[0])))
print('b_f = %.5f +- %.5f ' % (pars[1], np.sqrt(covm.diagonal()[1])))


#Calcoliamo il chi quadro,indice ,per quanto possibile, della bontà del fit:
chisq = sum(((y - f(x, *pars))/dy)**2.)
ndof = len(y) - len(pars)
print('chi quadro = %.3f (%d dof)' % (chisq, ndof))


#Definiamo un matrice di zeri che divverà la matrice di correlazione:
c=np.zeros((len(pars),len(pars)))
#Calcoliamo le correlazioni e le inseriamo nella matrice:
for i in range(0, len(pars)):
    for j in range(0, len(pars)):
       c[i][j]=(covm[i][j])/(np.sqrt(covm.diagonal()[i])*np.sqrt(covm.diagonal()[j]))
print(c) #matrice di correlazione


#Grafichiamo il risultato
fig1 = plt.figure(3, figsize=(14,6))
plt.suptitle('limite al continuo 2D',fontsize=15)
#Parte superiore contenetnte il fit:
frame1=fig1.add_axes((.1,.35,.4,.6))
#frame1=fig1.add_axes((trasla lateralmente, trasla verticamente, larghezza, altezza))
frame1.set_title('con over-relaxation',fontsize=10)
plt.ylabel('densità di energia',fontsize=10)
plt.grid()


plt.errorbar(x, y, dy, fmt='.', color='black', label='dati') #grafico i punti
plt.errorbar(m[:2], dens1[:2], d_dens1[:2], fmt='.', color='red', label='outliers')
t=np.linspace(0,np.max(m), 10000)

plt.plot(t, f(t, *pars), color='blue', alpha=0.5, label='best fit') #grafico del best fit
plt.legend(loc='best')#inserisce la legenda nel posto migliorte


#Parte inferiore contenente i residui
frame2=fig1.add_axes((.1,.1,.4,.2))

#Calcolo i residui normalizzari
ff=(y-f(x, *pars))/dy
frame2.set_ylabel('Residui Normalizzati')
plt.xlabel('massa',fontsize=10)
#plt.ticklabel_format(axis = 'both', style = 'sci', scilimits = (0,0))

plt.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
plt.plot(x, ff, '.', color='black') #grafico i residui normalizzati
plt.grid()

##=============================
## SENZA OVER-RELAXATION

nt, eno, deno = np.loadtxt('densmo.dat', unpack=True)

dens1 = 0.5*nt**2*eno
d_dens1 = 0.5*nt**2*deno

m = 1/(10*nt)
x = m[1:]
y = dens1[1:]
dy = d_dens1[1:]


#Eseguiamo il fit e stampiamo i risultati:
print('\n')
print('metodo senza over-relaxation')
print(f'a_t = {np.pi/6}')
pars, covm = curve_fit(f, x, y, sigma=dy, absolute_sigma=False)
print('a_f = %.5f +- %.5f ' % (pars[0], np.sqrt(covm.diagonal()[0])))
print('b_f = %.5f +- %.5f ' % (pars[1], np.sqrt(covm.diagonal()[1])))


#Calcoliamo il chi quadro,indice ,per quanto possibile, della bontà del fit:
chisq = sum(((y - f(x, *pars))/dy)**2.)
ndof = len(y) - len(pars)
print('chi quadro = %.3f (%d dof)' % (chisq, ndof))


#Definiamo un matrice di zeri che divverà la matrice di correlazione:
c=np.zeros((len(pars),len(pars)))
#Calcoliamo le correlazioni e le inseriamo nella matrice:
for i in range(0, len(pars)):
    for j in range(0, len(pars)):
       c[i][j]=(covm[i][j])/(np.sqrt(covm.diagonal()[i])*np.sqrt(covm.diagonal()[j]))
print(c) #matrice di correlazione


#Grafichiamo il risultato
fig1 = plt.figure(3, figsize=(14,6))
#Parte superiore contenetnte il fit:
frame1=fig1.add_axes((.55,.35,.4,.6))
#frame1=fig1.add_axes((trasla lateralmente, trasla verticamente, larghezza, altezza))
frame1.set_title('senza over-relaxation',fontsize=10)
plt.ylabel('densità di energia',fontsize=10)
plt.grid()


plt.errorbar(x, y, dy, fmt='.', color='black', label='dati') #grafico i punti
plt.errorbar(m[:1], dens1[:1], d_dens1[:1], fmt='.', color='red', label='outliers')
t=np.linspace(0,np.max(m), 10000)

plt.plot(t, f(t, *pars), color='blue', alpha=0.5, label='best fit') #grafico del best fit
plt.legend(loc='best')#inserisce la legenda nel posto migliorte


#Parte inferiore contenente i residui
frame2=fig1.add_axes((.55,.1,.4,.2))

#Calcolo i residui normalizzari
ff=(y-f(x, *pars))/dy
frame2.set_ylabel('Residui Normalizzati')
plt.xlabel('massa',fontsize=10)
#plt.ticklabel_format(axis = 'both', style = 'sci', scilimits = (0,0))

plt.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
plt.plot(x, ff, '.', color='black') #grafico i residui normalizzati
plt.grid()

"""
##=========================================================================================
"""
##Metodo dell'anomalia di traccia

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
