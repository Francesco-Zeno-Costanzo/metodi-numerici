import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

## Acthung:
# Per calcolare lo splitting dei livelli
# viene letto il file di input della simulazione il che potrebbe
# dare problemi se la correllazione di x e di x^2 dossero fatte
# con parametri diversi, quindi a meno di non creare due file
# risulta opportuno eseguire il codice a blocchi
##      

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
        
plt.figure(1)
plt.title('Andamento della correlazione di x a vari eta')
plt.ylabel('C(k)')
plt.xlabel('k')
plt.grid()

# dato che le simulazioni sono fatte a N*eta costante
# la lungezza dei vettori a vari eta cambia ed è quindi
# necessario leggere un correlazione per volta

for j in range(nret):
    #seleziono una singloa righa del file, i.e. una sola correlazione
    a = np.loadtxt(r'datiplot/datiyc.dat', skiprows=j, max_rows=1)
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

  
fig = plt.figure(2)
frame1 = fig.add_axes((.1,.35,.8,.6))
frame1.set_title('Splitting tra fondamentale e primo eccitato')
frame1.set_ylabel('$(E_1 - E_0)/\hbar \omega$')
frame1.grid()


frame1.errorbar(eta**2, e/eta, de/eta, fmt='.', color='red', label='outliers')
frame1.errorbar(x, y, dy, fmt='.', color='blue', label='dati')
t = np.linspace(np.min(eta**2), np.max(eta**2), 1000)
frame1.plot(t, g(t, *pars), 'k', label='best fit')   
frame1.legend(loc='best')

frame2 = fig.add_axes((.1,.1,.8,.2))

ff = (y - g(x, *pars))/dy
frame2.set_ylabel('Residui Normalizzati')
frame2.set_xlabel('$\eta^2$')

frame2.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
frame2.plot(x, ff, '.', color='black') #grafico i residui normalizzati
frame2.grid()

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
        
plt.figure(3)
plt.title('Andamento della correlazione di x a vari eta')
plt.ylabel('C(k)')
plt.xlabel('k')
plt.grid()

# dato che le simulazioni sono fatte a N*eta costante
# la lungezza dei vettori a vari eta cambia ed è quindi
# necessario leggere un correlazione per volta

for j in range(nret):
    #seleziono una singloa righa del file, i.e. una sola correlazione
    a = np.loadtxt(r'datiplot/datiy2c.dat', skiprows=j, max_rows=1)
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



fig2 = plt.figure(4)
frame3 = fig2.add_axes((.1,.35,.8,.6))
frame3.set_title('Splitting tra fondamentale e secondo eccitato')
frame3.set_ylabel('$(E_2 - E_0)/\hbar \omega$')
frame3.grid()

frame3.errorbar(eta**2, e/eta, de/eta, fmt='.', color='red', label='outliers')
frame3.errorbar(x, y, dy, fmt='.', color='blue', label='dati')
t = np.linspace(np.min(x), np.max(x), 1000)
frame3.plot(t, g(t, *pars), 'k', label='best fit')   
frame3.legend(loc='best')
 
frame4 = fig2.add_axes((.1,.1,.8,.2))

ff = (y - g(x, *pars))/dy
frame4.set_ylabel('Residui Normalizzati')
frame4.set_xlabel('$\eta^2$')

frame4.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
frame4.plot(x, ff, '.', color='black') #grafico i residui normalizzati
frame4.grid()


plt.show()
