import numpy as np
import scipy.special as ssp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

y2, dy2, Dy2, dDy2 = np.loadtxt(r'datiy2plot.dat', unpack=True)

eta = np.array([1/l for l in range(1, len(y2)+1)])

def f(x, a, b):
    return a + b*x**2
    

#Eseguiamo il fit e stampiamo i risultati:
pars, covm = curve_fit(f, eta, y2, sigma=dy2, absolute_sigma=False)
print('a = %.5f +- %.5f ' % (pars[0], np.sqrt(covm.diagonal()[0])))
print('b = %.5f +- %.5f ' % (pars[1], np.sqrt(covm.diagonal()[1])))


#Calcoliamo il chi quadro,indice ,per quanto possibile, della bontà del fit:
chisq = sum(((y2 - f(eta, *pars))/dy2)**2.)
ndof = len(y2) - len(pars)
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


plt.errorbar(eta, y2, dy2, fmt='.', color='black', label='dati') #grafico i punti
t=np.linspace(0,np.max(eta), 10000)

plt.plot(t, f(t, *pars), color='blue', alpha=0.5, label='best fit') #grafico del best fit
plt.legend(loc='best')#inserisce la legenda nel posto migliorte


#Parte inferiore contenente i residui
frame2=fig1.add_axes((.1,.1,.8,.2))

#Calcolo i residui normalizzari
ff=(y2-f(eta, *pars))/dy2
frame2.set_ylabel('Residui Normalizzati')
plt.xlabel('$\eta$',fontsize=10)
#plt.ticklabel_format(axis = 'both', style = 'sci', scilimits = (0,0))

x1=np.linspace(0,np.max(eta), 1000)
plt.plot(x1, 0*x1, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
plt.plot(eta, ff, '.', color='black') #grafico i residui normalizzati
plt.grid()

##========================================================================

a = np.loadtxt(r'datiycplot.dat')
K = len(a[0,:])//2
W = len(a[:,0])
e = np.zeros(W)
de = np.zeros(W)
x = np.linspace(0, K-1, K)
eta = np.array([1/i for i in range(1, W+1)])

colors = plt.cm.jet(np.linspace(0, 1, W))

def f(x, A, de):
        return A*np.exp(-x*de)
        
plt.figure(2)
plt.title('Andamento della correlazione di x a vari eta')
plt.ylabel('C(k)')
plt.xlabel('k')
plt.grid()

for j in range(W):
    
    y = a[j, 0:K]
    dy = a[j, K:] 
    
    plt.errorbar(x, y, dy, fmt='.', ecolor=colors[j], label=fr'$\eta$ = {eta[j]:.2f}')
    plt.legend(loc='best')
   
    pars, covm = curve_fit(f, x, y, sigma=dy)
    chisq = sum(((y - f(x, *pars))/dy)**2.)
    ndof = len(y) - len(pars)
    print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
    
    t = np.linspace(np.min(x), np.max(x), 1000)
    plt.plot(t, f(t, *pars), color=colors[j])
    e[j] = pars[1]
    de[j] = np.sqrt(covm.diagonal()[1])

def g(x, m, q):
    return m*x + q
 
x = eta**2
y = e/eta
dy = de/eta  
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

  
fig3 = plt.figure(3)
frame3 = fig3.add_axes((.1,.35,.8,.6))
frame3.set_title('Splitting tra fondamentale e primo eccitato')
frame3.set_ylabel('$(E_1 - E_0)/\hbar \omega$')
frame3.grid()

frame3.errorbar(x, y, dy, fmt='.')
t = np.linspace(np.min(x), np.max(x), 1000)
frame3.plot(t, g(t, *pars))   
 
frame4 = fig3.add_axes((.1,.1,.8,.2))

ff = (y - g(x, *pars))/dy
frame4.set_ylabel('Residui Normalizzati')
frame4.set_xlabel('$\eta^2$')

frame4.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
frame4.plot(x, ff, '.', color='black') #grafico i residui normalizzati
frame4.grid()

##=======================================================================

cammini = open(r'cammino.dat', "r").read()
cammini = cammini.split('\n')

N = int(cammini[0])
P = int(float(cammini[1]))
eta = np.array([1/i for i in range(1, N+1)])

cammini = cammini[2:N+2]

plt.figure(4)
plt.suptitle('Cammini a vari eta')
i = 1
f = np.array([])
for C in cammini:
    plt.subplot(3, 2, i)
    plt.grid()
    K = C.split()
    
    cammino = np.array([])
    for c in K:
        cammino = np.insert(cammino, len(cammino), float(c))
        f = np.insert(f, len(f), float(c))
    x = np.linspace(0, P, len(cammino))
    plt.plot(x, cammino, label=fr'$\eta$={eta[i-1]:.2f}')
    plt.legend(loc='best')
    i += 1


def G(x, m):
    return (1/(np.pi)**(1/4))*(1/np.sqrt((2**m)*ssp.gamma(m+1)))*ssp.eval_hermite(m, x)*np.exp(-(x**2)/2)
plt.figure(5)

plt.hist(f, bins=20, density=True)

x=np.linspace(-2,2, 10000)
plt.plot(x, abs(G(x, 0))**2)
plt.grid()


##=======================================================================
a = np.loadtxt(r'datiy2cplot.dat')
K = len(a[0,:])//2
W = len(a[:,0])
e = np.zeros(W)
de = np.zeros(W)
x = np.linspace(0, K-1, K)
eta = np.array([1/i for i in range(1, W+1)])

colors = plt.cm.jet(np.linspace(0, 1, W))

def f(x, A, de):
        return A*np.exp(-x*de)
plt.figure(6)
plt.title('Andamento della correlazione di x^2 a vari eta')
plt.ylabel('C(k)')
plt.xlabel('k')
plt.grid()
 
for j in range(W):
    
    y = a[j, 0:K]
    dy = a[j, K:] 
    
    plt.errorbar(x, y, dy, fmt='.', ecolor=colors[j], label=fr'$\eta$ = {eta[j]:.2f}')
    plt.legend(loc='best')
	
    pars, covm = curve_fit(f, x, y, sigma=dy)
    chisq = sum(((y - f(x, *pars))/dy)**2.)
    ndof = len(y) - len(pars)
    print('chi quadro = %.3f (%d dof)' % (chisq, ndof))
    
    t = np.linspace(np.min(x), np.max(x), 1000)
    plt.plot(t, f(t, *pars), color=colors[j])
    e[j] = pars[1]
    de[j] = np.sqrt(covm.diagonal()[1])

def g(x, m, q):
    return m*x + q
 
x = eta**2
y = e/eta
dy = de/eta  
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



fig6 = plt.figure(7)
frame6 = fig6.add_axes((.1,.35,.8,.6))
frame6.set_title('Splitting tra fondamentale e secondo eccitato')
frame6.set_ylabel('$(E_2 - E_0)/\hbar \omega$')
frame6.grid()

frame6.errorbar(x, y, dy, fmt='.')
t = np.linspace(np.min(x), np.max(x), 1000)
frame6.plot(t, g(t, *pars))   
 
frame7 = fig6.add_axes((.1,.1,.8,.2))

ff = (y - g(x, *pars))/dy
frame7.set_ylabel('Residui Normalizzati')
frame7.set_xlabel('$\eta^2$')

frame7.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
frame7.plot(x, ff, '.', color='black') #grafico i residui normalizzati
frame7.grid()


plt.show()
