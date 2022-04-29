import numpy as np
import scipy.special as ssp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

## Acthung:
# nella sezione dove viene calcolato lo splitting dei livelli
# viene letto il file di input della simulazione il che potrebbe
# dare problemi se la correllazione di x e di x^2 dossero fatte
# con parametri diversi, quindi a meno di non creare due file
# risulta opportuno eseguire il codice a blocchi
##


##ENERGIA CINETICA
y2, dy2, Dy2, dDy2 = np.loadtxt(r'datieneplot.dat', unpack=True)

eta = np.array([1/l for l in range(1, len(y2)+1)])

def f(x, a, b):
    ''' funzione di fit
    '''
    return a + b*x**2

y = y2
x = eta
dy = dy2 

#Eseguiamo il fit e stampiamo i risultati:
pars, covm = curve_fit(f, x, y, sigma=dy, absolute_sigma=False)
print('a = %.5f +- %.5f ' % (pars[0], np.sqrt(covm.diagonal()[0])))
print('b = %.5f +- %.5f ' % (pars[1], np.sqrt(covm.diagonal()[1])))


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
fig1 = plt.figure(1)
#Parte superiore contenetnte il fit:
frame1=fig1.add_axes((.1,.35,.8,.6))
#frame1=fig1.add_axes((trasla lateralmente, trasla verticamente, larghezza, altezza))
frame1.set_title('Fit del valore medio di $y^2$',fontsize=20)
plt.ylabel('$<y^2>$',fontsize=10)
#plt.ticklabel_format(axis = 'both', style = 'sci', scilimits = (0,0))#notazione scientifica sugliassi
plt.grid()


plt.errorbar(x, y, dy, fmt='.', color='black', label='dati') #grafico i punti
t=np.linspace(0,np.max(x), 10000)

plt.plot(t, f(t, *pars), color='blue', alpha=0.5, label='best fit') #grafico del best fit
plt.legend(loc='best')#inserisce la legenda nel posto migliorte


#Parte inferiore contenente i residui
frame2=fig1.add_axes((.1,.1,.8,.2))

#Calcolo i residui normalizzari
ff=(y-f(x, *pars))/dy
frame2.set_ylabel('Residui Normalizzati')
plt.xlabel('$\eta$',fontsize=10)
#plt.ticklabel_format(axis = 'both', style = 'sci', scilimits = (0,0))

x1=np.linspace(0,np.max(x), 1000)
plt.plot(x1, 0*x1, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
plt.plot(eta, ff, '.', color='black') #grafico i residui normalizzati
plt.grid()
"""
##========================================================================
"""
##ENERGIA POTENZIALE

plt.figure(2)
plt.subplot(121)
plt.title('Termine cinetico non normalizzato')
plt.xlabel(r'$\eta$')
plt.errorbar(eta, -Dy2/(2*eta**2), dDy2/(2*eta**2), fmt='.', color='black')
plt.grid()

plt.subplot(122)
plt.title('Termine cinetico normalizzato')
plt.xlabel(r'$\eta$')
plt.errorbar(eta, 1/(2*eta) - Dy2/(2*eta**2), dDy2/(2*eta**2), fmt='.', color='black')
plt.grid()
"""
##========================================================================
"""
##FUNZIONA D'ONDA
#i cammini vengono letti come una matrice
cammini01 = np.loadtxt(r'cammini1.dat')
cammini15 = np.loadtxt(r'cammini15.dat')
N01 = 20*1
N15 = 20*15
#trasformo la matrice in un vettore dato che nono importa l'ordine
f01 = np.reshape(cammini01, int(N01*1e5))
f15 = np.reshape(cammini15, int(N15*1e5))

cammini01 = np.loadtxt(r'cammino1.dat')
cammini15 = np.loadtxt(r'cammino15.dat')
N01 = 3*1
N15 = 3*15
g01 = np.reshape(cammini01, int(N01*1e5))
g15 = np.reshape(cammini15, int(N15*1e5))

def G(x, m):
    '''soluzione analitica oscillatore armonico
    '''
    return (1/(np.pi)**(1/4))*(1/np.sqrt((2**m)*ssp.gamma(m+1)))*ssp.eval_hermite(m, x)*np.exp(-(x**2)/2)

plt.figure(3)
plt.suptitle("Funzione d'onda stato fondamentale", fontsize=13)
x = np.linspace(-3,3, 10000)
plt.subplot(221)
plt.hist(f01, bins=60, density=True, label=r'$\beta \omega$ = 20; $\eta$=1')
plt.plot(x, abs(G(x, 0))**2, 'k', label='Soluzione analitica')
plt.legend(loc='best')
plt.grid()

plt.subplot(222)
plt.hist(f15, bins=60, density=True, label=r'$\beta \omega$ = 20; $\eta$=0.07')
plt.plot(x, abs(G(x, 0))**2, 'k', label='Soluzione analitica')
plt.legend(loc='best')
plt.grid()

plt.subplot(223)
plt.hist(g01, bins=60, density=True, label=r'$\beta \omega$ = 3; $\eta$=1')
plt.plot(x, abs(G(x, 0))**2, 'k', label='Soluzione analitica')
plt.legend(loc='best')
plt.grid()

plt.subplot(224)
plt.hist(g15, bins=60, density=True, label=r'$\beta \omega$ = 3; $\eta$=0.07')
plt.plot(x, abs(G(x, 0))**2, 'k', label='Soluzione analitica')
plt.legend(loc='best')
plt.grid()

"""
##========================================================================

"""
#ENERGIA INTERNA

y2, dy2, Dy2, dDy2 = np.loadtxt(r'datieneUplot.dat', unpack=True)
eta = 0.1
#energia corettemente rinormalizzata
U = 1/(2*eta) - 1/(2*eta**2)*Dy2 + y2/2
dU = dy2/2 + dDy2/(2*eta**2)
t = np.array([l for l in range(1, len(y2)+1)])
t = 1/t

def E(x, a, b):
    ''' funzione di fit
    '''
    return a + b/(np.exp(1/x)-1)

pars, covm = curve_fit(E, t, U, sigma=dU, absolute_sigma=False)
print('a = %.5f +- %.5f ' % (pars[0], np.sqrt(covm.diagonal()[0])))
print('b = %.5f +- %.5f ' % (pars[1], np.sqrt(covm.diagonal()[1])))

chisq = sum(((U - E(t, *pars))/dU)**2.)
ndof = len(U) - len(pars)
print('chi quadro = %.3f (%d dof)' % (chisq, ndof))


#Definiamo un matrice di zeri che divverà la matrice di correlazione:
c=np.zeros((len(pars),len(pars)))
#Calcoliamo le correlazioni e le inseriamo nella matrice:
for i in range(0, len(pars)):
    for j in range(0, len(pars)):
       c[i][j]=(covm[i][j])/(np.sqrt(covm.diagonal()[i])*np.sqrt(covm.diagonal()[j]))
print(c) #matrice di correlazione

fig8 = plt.figure(4)
frame8 = fig8.add_axes((.1,.35,.8,.6))
frame8.set_title('Energia interna in funzione della temperatura')
frame8.set_ylabel('Energia interna')
frame8.grid()
frame8.errorbar(t, U, dU, fmt='.', color='black', label='dati')
p = np.linspace(np.min(t), np.max(t), 1000)
frame8.plot(p, E(p, *pars), 'b', label='best fit')   
frame8.legend(loc='best')

frame9 = fig8.add_axes((.1,.1,.8,.2))
ff = (U - E(t, *pars))/dU
frame9.set_ylabel('Residui Normalizzati')
frame9.set_xlabel(r'$1/(\beta \omega)$')
frame9.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
frame9.plot(t, ff, '.', color='black') #grafico i residui normalizzati
frame9.grid()


"""
##========================================================================
"""
##CORRELAZIONE DI X

A = np.loadtxt(r'input.txt', unpack=True)
_,_,P,_,_,nret,_ , shift = A # leggo i dati necessari
# sono interi ma vengono letti come reali
P = int(P); nret = int(nret); shift = int(shift)

#vettori che serviranno poi
e = np.zeros(nret)
de = np.zeros(nret)
eta = np.array([1/i for i in range(1+shift, nret+shift+1)])
colors = plt.cm.jet(np.linspace(0, 1, nret))


def f(x, A, de):
    '''funzione di fit
    '''
    return A*np.exp(-x*de)
        
plt.figure(5)
plt.title('Andamento della correlazione di x a vari eta')
plt.ylabel('C(k)')
plt.xlabel('k')
plt.grid()

# dato che le simulazioni sono fatte a N*eta costante
# la lungezza dei vettori a vari eta cambia ed è quindi
# necessario leggere un correlazione per volta

for j in range(nret):
    #seleziono una singloa righa del file, i.e. una sola correlazione
    a = np.loadtxt(r'adatiycplot.dat', skiprows=j, max_rows=1)
    k = P*((j+1) + shift)//2 #lunghezza correlazione per ogni cammino
    
    #seleziono solo meta curva
    x = np.linspace(0, k//2 - 1,k//2)
    y = a[0:k//2]
    dy = a[k:k+k//2] 
    
    plt.errorbar(x, y, dy, fmt='.', color=colors[j], label=fr'$\eta$ = {eta[j]:.2f}')
    plt.legend(loc='best')
   
    pars, covm = curve_fit(f, x, y, sigma=dy)
    chisq = sum(((y - f(x, *pars))/dy)**2.)
    ndof = len(y) - len(pars)
    print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
    
    t = np.linspace(np.min(x), np.max(x), 1000)
    plt.plot(t, f(t, *pars), color=colors[j])
    #salvo i valori delle energie e rispettivi errori
    e[j] = pars[1]
    de[j] = np.sqrt(covm.diagonal()[1])
    

def g(x, m, q):
    '''funziine di fit
    '''
    return m*x + q
 
#parametrizzo la curva in eta^2 
#probabblmente dovuto agli errori sistematci
x = eta**2
#avendo fittato in funzione della lunghezza dell'array
#per avere le energie fisiche devo dividere per eta
y = e/eta
dy = de/eta

#sono outlier nella specifica simulazione
x=np.delete(x,3);y=np.delete(y,3);dy=np.delete(dy,3)
x=np.delete(x,3);y=np.delete(y,3);dy=np.delete(dy,3)
x=np.delete(x,4);y=np.delete(y,4);dy=np.delete(dy,4)

pars, covm = curve_fit(g, x, y, sigma=dy)
print('m = %.5f +- %.5f ' % (pars[0], np.sqrt(covm.diagonal()[0])))
print('q = %.5f +- %.5f ' % (pars[1], np.sqrt(covm.diagonal()[1])))

chisq = sum(((y - g(x, *pars))/dy)**2.)
ndof = len(y) - len(pars)
print('chi quadro = %.3f (%d dof)' % (chisq, ndof))


c=np.zeros((len(pars),len(pars)))

for i in range(0, len(pars)):
    for j in range(0, len(pars)):
       c[i][j]=(covm[i][j])/(np.sqrt(covm.diagonal()[i])*np.sqrt(covm.diagonal()[j]))
print(c) #matrice di correlazione

  
fig3 = plt.figure(6)
frame3 = fig3.add_axes((.1,.35,.8,.6))
frame3.set_title('Splitting tra fondamentale e primo eccitato')
frame3.set_ylabel('$(E_1 - E_0)/\hbar \omega$')
frame3.grid()


frame3.errorbar(eta**2, e/eta, de/eta, fmt='.', color='red', label='outliers')
frame3.errorbar(x, y, dy, fmt='.', color='blue', label='dati')
t = np.linspace(np.min(eta**2), np.max(eta**2), 1000)
frame3.plot(t, g(t, *pars), 'k', label='best fit')   
frame3.legend(loc='best')
frame4 = fig3.add_axes((.1,.1,.8,.2))

ff = (y - g(x, *pars))/dy
frame4.set_ylabel('Residui Normalizzati')
frame4.set_xlabel('$\eta^2$')

frame4.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
frame4.plot(x, ff, '.', color='black') #grafico i residui normalizzati
frame4.grid()

"""
##=======================================================================
"""

##CORRELAZIONE DI X^2

A = np.loadtxt(r'input.txt', unpack=True)
_,_,P,_,_,nret,_ , shift = A # leggo i dati necessari
# sono interi ma vengono letti come reali
P = int(P); nret = int(nret); shift = int(shift)

#vettori che serviranno poi
e = np.zeros(nret)
de = np.zeros(nret)
eta = np.array([1/i for i in range(1+shift, nret+shift+1)])
colors = plt.cm.jet(np.linspace(0, 1, nret))


def f(x, A, de):
    '''funzione di fit
    '''
    return A*np.exp(-x*de)
        
plt.figure(7)
plt.title('Andamento della correlazione di x a vari eta')
plt.ylabel('C(k)')
plt.xlabel('k')
plt.grid()

# dato che le simulazioni sono fatte a N*eta costante
# la lungezza dei vettori a vari eta cambia ed è quindi
# necessario leggere un correlazione per volta

for j in range(nret):
    #seleziono una singloa righa del file, i.e. una sola correlazione
    a = np.loadtxt(r'adatiy2cplot.dat', skiprows=j, max_rows=1)
    k = P*((j+1) + shift)//2 #lunghezza correlazione per ogni cammino
    
    #seleziono solo meta curva
    x = np.linspace(0, k//2 - 1,k//2)
    y = a[0:k//2]
    dy = a[k:k+k//2] 
    
    plt.errorbar(x, y, dy, fmt='.', color=colors[j], label=fr'$\eta$ = {eta[j]:.2f}')
    plt.legend(loc='best')
   
    pars, covm = curve_fit(f, x, y, sigma=dy)
    chisq = sum(((y - f(x, *pars))/dy)**2.)
    ndof = len(y) - len(pars)
    print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
    
    t = np.linspace(np.min(x), np.max(x), 1000)
    plt.plot(t, f(t, *pars), color=colors[j])
    #salvo i valori delle energie e rispettivi errori
    e[j] = pars[1]
    de[j] = np.sqrt(covm.diagonal()[1])
    

def g(x, m, q):
    '''funziine di fit
    '''
    return m*x + q
 
#parametrizzo la curva in eta^2 
#probabblmente dovuto agli errori sistematci
x = eta**2
#avendo fittato in funzione della lunghezza dell'array
#per avere le energie fisiche devo dividere per eta
y = e/eta
dy = de/eta

#sono outlier nella specifica simulazione
x=np.delete(x,-1);y=np.delete(y,-1);dy=np.delete(dy,-1)

pars, covm = curve_fit(g, x, y, sigma=dy)
print('m = %.5f +- %.5f ' % (pars[0], np.sqrt(covm.diagonal()[0])))
print('q = %.5f +- %.5f ' % (pars[1], np.sqrt(covm.diagonal()[1])))

chisq = sum(((y - g(x, *pars))/dy)**2.)
ndof = len(y) - len(pars)
print('chi quadro = %.3f (%d dof)' % (chisq, ndof))


c=np.zeros((len(pars),len(pars)))

for i in range(0, len(pars)):
    for j in range(0, len(pars)):
       c[i][j]=(covm[i][j])/(np.sqrt(covm.diagonal()[i])*np.sqrt(covm.diagonal()[j]))
print(c) #matrice di correlazione



fig6 = plt.figure(8)
frame6 = fig6.add_axes((.1,.35,.8,.6))
frame6.set_title('Splitting tra fondamentale e secondo eccitato')
frame6.set_ylabel('$(E_2 - E_0)/\hbar \omega$')
frame6.grid()

frame6.errorbar(eta**2, e/eta, de/eta, fmt='.', color='red', label='outliers')
frame6.errorbar(x, y, dy, fmt='.', color='blue', label='dati')
t = np.linspace(np.min(x), np.max(x), 1000)
frame6.plot(t, g(t, *pars), 'k', label='best fit')   
frame6.legend(loc='best')
 
frame7 = fig6.add_axes((.1,.1,.8,.2))

ff = (y - g(x, *pars))/dy
frame7.set_ylabel('Residui Normalizzati')
frame7.set_xlabel('$\eta^2$')

frame7.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
frame7.plot(x, ff, '.', color='black') #grafico i residui normalizzati
frame7.grid()



plt.show()
