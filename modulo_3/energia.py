import numpy as np
import scipy.special as ssp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import grafici


##ENERGIA CINETICA
y2, dy2, Dy2, dDy2 = np.loadtxt(r'datiplot/datiene.dat', unpack=True)
eta = np.array([1/l for l in range(1, len(y2)+1)])

def f(x, a, b):
    ''' funzione di fit
    '''
    return a + b*x**2

y = y2
x = eta
dy = dy2
init = np.array([1, 1])

print('Valor medio di y^2')
Title = 'Fit del valore medio di $y^2$'
ylabel = '$<y^2>$'
xlabel = '$\eta$'
grafici.fit(f, x, y, dy, init, 1, Title, xlabel, ylabel)
print('\n')

##========================================================================

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

##========================================================================

##FUNZIONA D'ONDA
#i cammini vengono letti come una matrice
cammini01 = np.loadtxt(r'datiplot/cammini1.dat')
cammini15 = np.loadtxt(r'datiplot/cammini15.dat')
N01 = 20*1
N15 = 20*15
#trasformo la matrice in un vettore dato che nono importa l'ordine
f01 = np.reshape(cammini01, int(N01*1e5))
f15 = np.reshape(cammini15, int(N15*1e5))

cammini01 = np.loadtxt(r'datiplot/cammino1.dat')
cammini15 = np.loadtxt(r'datiplot/cammino15.dat')
N01 = 3*1
N15 = 3*15
g01 = np.reshape(cammini01, int(N01*1e5))
g15 = np.reshape(cammini15, int(N15*1e5))

def G(x, m):
    '''soluzione analitica oscillatore armonico
    '''
    psi = (1/(np.pi)**(1/4))*(1/np.sqrt((2**m)*ssp.gamma(m+1)))*ssp.eval_hermite(m, x)*np.exp(-(x**2)/2)
    return psi
    
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


##========================================================================


#ENERGIA INTERNA

print('fit per il limite al continuo:')

B = np.array([3, 4, 5, 6, 7, 8, 10, 15, 20, 30, 40, 60], dtype=int)
U = np.zeros(len(B))
dU = np.zeros(len(B))

def f(x, a, b):
    ''' funzione di fit
    '''
    return a + b*x**2

Y = np.zeros((len(B), 10))
X = np.zeros((len(B), 10))
P = np.zeros((len(B), 2))
dY = np.zeros((len(B), 10))
   
for i, b in enumerate(B):
    y2, dy2, Dy2, dDy2 = np.loadtxt(fr'datiplot/ene{b}.dat', unpack=True)
    eta = np.array([1/l for l in range(1, len(y2)+1)])

    y = 1/(2*eta) - Dy2/(2*eta**2) +y2/2 ; Y[i,:] = y
    x = eta; X[i,:] = x
    dy = dy2/2 + dDy2/(2*eta**2); dY[i,:] = dy
    init = np.array([1, 1])

    Title = "Fit del valore medio dell'energia"
    ylabel = '$U$'
    xlabel = '$\eta$'
    pars, dpars = grafici.fit(f, x, y, dy, init, 1, Title, xlabel, ylabel, False)
    P[i,:] = pars
    U[i] = pars[0]
    dU[i] = dpars[0]
    print('\n')

print('fit energia in funzione del tempo')    
Title = "Fit del valore medio dell'energia"
ylabel = '$U$'
xlabel = '$\eta^2$'
grafici.plotfit(f, P, X, Y, dY, B, 4, Title, xlabel, ylabel)

def E(x, a, b):
    ''' funzione di fit
    '''
    return a + b/(np.exp(1/x)-1)

Title = 'Energia interna in funzione della temperatura'
ylabel = 'Energia interna'
xlabel = r'$1/(\beta \omega)$'
init = np.array([0.5, 1])    
pars, dpars = grafici.fit(E, 1/B, U, dU, init, 5, Title, xlabel, ylabel)   

plt.show()
