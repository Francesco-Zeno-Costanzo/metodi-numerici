import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#il nuero finale identifica la grandezza del reticolo
E15, M15, C15, X15, cb15, dE15, dM15, dC15, dX15, dcb15 = np.loadtxt(r'datiplot15.dat', unpack=True)
E20, M20, C20, X20, cb20, dE20, dM20, dC20, dX20, dcb20 = np.loadtxt(r'datiplot20.dat', unpack=True)
E25, M25, C25, X25, cb25, dE25, dM25, dC25, dX25, dcb25 = np.loadtxt(r'datiplot25.dat', unpack=True)
E30, M30, C30, X30, cb30, dE30, dM30, dC30, dX30, dcb30 = np.loadtxt(r'datiplot30.dat', unpack=True)
E35, M35, C35, X35, cb35, dE35, dM35, dC35, dX35, dcb35 = np.loadtxt(r'datiplot35.dat', unpack=True)
E40, M40, C40, X40, cb40, dE40, dM40, dC40, dX40, dcb40 = np.loadtxt(r'datiplot40.dat', unpack=True)
E45, M45, C45, X45, cb45, dE45, dM45, dC45, dX45, dcb45 = np.loadtxt(r'datiplot45.dat', unpack=True)
E50, M50, C50, X50, cb50, dE50, dM50, dC50, dX50, dcb50 = np.loadtxt(r'datiplot50.dat', unpack=True)


H = 0
#valori esatti per gli indici e la temperatura critica
#n = 1
#b = 1/8
#g = 7/4
#a = 0
#bc = 0.4406868

##Alcuni plot a titolo espositivo

B = np.zeros(len(E15))
for i in range(len(B)):
    B[i] = 0.35 + i*(0.5-0.35)/len(B)

plt.figure(1)
plt.title(f"Simulazione del modello di Ising 2D tramite Metropolis \n Campo magnetico esterno B={H}")
plt.xlabel(r"$\beta$ [u.a.]", fontsize=15)
plt.ylabel("Energia [u.a.]", fontsize=15)
plt.grid()
plt.errorbar(B, E20, dE20, fmt='.', color='black', label='L=20')
plt.errorbar(B, E30, dE30, fmt='v', color='black', label='L=30')
plt.errorbar(B, E40, dE40, fmt='*', color='black', label='L=40')
plt.errorbar(B, E50, dE50, fmt='^', color='black', label='L=50')
plt.legend(loc='best')


plt.figure(2)
plt.title(f"Simulazione del modello di Ising 2D tramite Metropolis \n Campo magnetico esterno B={H}")
plt.xlabel(r"$\beta$ [u.a.]", fontsize=15)
plt.ylabel("Magnetizzazione [u.a.]", fontsize=15)
plt.grid()
plt.errorbar(B, M20, dM20, fmt='.', color='black', label='L=20')
plt.errorbar(B, M30, dM30, fmt='v', color='black', label='L=30')
plt.errorbar(B, M40, dM40, fmt='*', color='black', label='L=40')
plt.errorbar(B, M50, dM50, fmt='^', color='black', label='L=50')
plt.legend(loc='best')


plt.figure(3)
plt.title(f"Simulazione del modello di Ising 2D tramite Metropolis \n Campo magnetico esterno B={H}")
plt.xlabel(r"$\beta$ [u.a.]", fontsize=15)
plt.ylabel("Calore specifico [u.a.]", fontsize=15)
plt.grid()
plt.errorbar(B, C20, dC20, fmt='.', color='black', label='L=20')
plt.errorbar(B, C30, dC30, fmt='v', color='black', label='L=30')
plt.errorbar(B, C40, dC40, fmt='*', color='black', label='L=40')
plt.errorbar(B, C50, dC50, fmt='^', color='black', label='L=50')
plt.legend(loc='best')


plt.figure(4)
plt.title(f"Simulazione del modello di Ising 2D tramite Metropolis \n Campo magnetico esterno B={H}")
plt.xlabel(r"$\beta$ [u.a.]", fontsize=15)
plt.ylabel("Suscettività [u.a.]", fontsize=15)
plt.grid()
plt.errorbar(B, X20, dX20, fmt='.', color='black', label='L=20')
plt.errorbar(B, X30, dX30, fmt='v', color='black', label='L=30')
plt.errorbar(B, X40, dX40, fmt='*', color='black', label='L=40')
plt.errorbar(B, X50, dX50, fmt='^', color='black', label='L=50')
plt.legend(loc='best')


##plot del cumulante di binder da cui si può stimare il beta critico

fig, main_ax = plt.subplots() #creo le variabili per il garfico
main_ax.set_title(f"Simulazione del modello di Ising 2D tramite Metropolis \n Campo magnetico esterno B={H}")
main_ax.set_xlabel(r"$\beta$ [u.a.]", fontsize=15)
main_ax.set_ylabel("Cumulante di binder", fontsize=15)
main_ax.errorbar(B, cb20, dcb20, fmt='.', color='black', label='L=20')
main_ax.errorbar(B, cb30, dcb30, fmt='v', color='black', label='L=30')
main_ax.errorbar(B, cb40, dcb40, fmt='*', color='black', label='L=40')
main_ax.errorbar(B, cb50, dcb50, fmt='^', color='black', label='L=50')
main_ax.legend(loc='best')
main_ax.grid()
#aggiungo un sotto-grafico in alto a destra che servirà per osservare meglio il punto critico
right_inset_ax = fig.add_axes([.55, .5, .3, .3])
right_inset_ax.errorbar(B, cb20, dcb20, fmt='.', linestyle='--', color='black', label='L=20')
right_inset_ax.errorbar(B, cb30, dcb30, fmt='v', linestyle='--', color='blue', label='L=30')
right_inset_ax.errorbar(B, cb40, dcb40, fmt='*', linestyle='--', color='red', label='L=40')
right_inset_ax.errorbar(B, cb50, dcb50, fmt='^', linestyle='--', color='green', label='L=50')
right_inset_ax.legend(loc='best')
right_inset_ax.grid()
right_inset_ax.set_xlim(0.436, 0.444)
right_inset_ax.set_ylim(1.05, 1.25)

#valore di beta ricavato
bc = 0.440
dbc = 0.003
print(f"beta critico = {bc:.3f} +- {dbc:.3f}")

##grafico del comportamento della magnetizzazione nelle due fais della transizione

mag, ene = np.loadtxt(r'dati20hist.dat', skiprows=3, unpack=True)
x = mag[0:40000]/(20*20)        #fase disordinata
z = mag[380000:410000]/(20*20)  #fase ordinata

plt.figure(6)
plt.suptitle('Comportamneto della magnetizzazione reticolo 20X20', fontsize=20)
plt.subplot(221)
plt.title('Fase disordinata', fontsize=15)
plt.xlabel('Magnetizzazione', fontsize=15)
plt.ylabel('P(M)', fontsize=15)
plt.grid()
plt.hist(x, 100, histtype='step', density=True)
plt.subplot(222)
plt.title('Fase ordinata', fontsize=15)
plt.xlabel('Magnetizzazione', fontsize=15)
plt.ylabel('P(M)', fontsize=15)
plt.grid()
plt.hist(z, 200, histtype='step', density=True)
plt.subplot(223)
plt.ylabel('Magnetizzazione', fontsize=15)
plt.xlabel('Passi della catena', fontsize=15)
plt.grid()
plt.plot(x)
plt.subplot(224)
plt.ylabel('Magnetizzazione', fontsize=15)
plt.xlabel('Passi della catena', fontsize=15)
plt.grid()
t = np.linspace(380000, 410000, 30000)
plt.plot(t, z)


def FX(x, m, g, q):
    '''funzione modello per i fit successivi, tranne per il calore specifico
    '''
    return m*x**g + q

def ChiRes(y, dy, x, p, f, k=0):
    ''' 
    Funzione per il calcolo del chi quadro e residui
    se viene passato k uguale ad uno restituisce  i residui
    '''
    res = (y - f(x, *p))/dy
    chisq = sum(res**2.) #chi quadro
    ndof = len(y) - len(p) #gradi di libertà
    if k == 1:
        return res
    else:
        return chisq, ndof

def Corr(covm, n):
    ''' funzione per il calcolo della matrice di correlazione
    '''
    corr = np.zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            corr[i][j] = (covm[i][j])/(np.sqrt(covm.diagonal()[i])*np.sqrt(covm.diagonal()[j]))
    return corr

##stima di gamma/nu
print('stima di gamma/nu')

#funzione per estrarre il massimo della suscettività e relativo errore
def MS(*arg):
    '''*arg permette alla funzione di prendere in input infiniti argomenti.
       La funzione è scritta affinche la prima metà degli argomenti passati siano
       gli array dei valori centrali mentre la seconda parte quelli degli errori.
       Restitusce due array contenenti il massimo dei valori centrali e relativi errori
    '''
    Q = len(arg)//2 #lunghezza degli array
    l = np.zeros(Q)
    dl = np.zeros(Q)
    k = np.zeros(Q, dtype=np.int64)#verra usato come indice quindi deve essere intero

    for i in range(len(arg)):

        if i < Q:
            l[i] = np.max(arg[i]) #trovo il massimo
            k[i] = np.where(arg[i] == l[i])[0][0] #conservo il rispettivo indice

        else:
            #trovo il valore dell'errore del punto usando l'indice
            dl[i-Q] = arg[i][k[i-Q]]

    return l, dl

MX, dMX = MS(X15, X20, X25, X30, X35, X40, X45, X50, dX15, dX20, dX25, dX30, dX35, dX40, dX45, dX50)
N = np.arange(15, 51, 5)

pars, covm = curve_fit(FX, N, MX, sigma=dMX, absolute_sigma=False)
print(f'm = {pars[0]:.5f} +- {np.sqrt(covm.diagonal()[0]):.5f}')
print(f'g = {pars[1]:.5f} +- {np.sqrt(covm.diagonal()[1]):.5f}')
print(f'q = {pars[2]:.5f} +- {np.sqrt(covm.diagonal()[2]):.5f}')


chi1, Nd1 = ChiRes(MX, dMX, N, pars, FX)
print(f'chi quadro = {chi1:.3f}  ({Nd1} dof)')

c1 = Corr(covm, len(pars))
print(c1) 

fig1 = plt.figure(7)

frame1 = fig1.add_axes((.1, .35, .8, .6))
#frame1=fig1.add_axes((trasla lateralmente, trasla verticamente, larghezza, altezza))
frame1.set_title('Massimo della suscettività al variare di L', fontsize=20)
frame1.set_ylabel('$\chi_{max} [a.u.]$', fontsize=15)
plt.grid()

plt.errorbar(N, MX, dMX, fmt='.', color='black', label='dati') 
t = np.linspace(np.min(N), np.max(N), 10000)

plt.plot(t, FX(t, *pars), color='blue', label='best fit') 
plt.legend(loc='best')


#Parte inferiore contenente i residui
frame2 = fig1.add_axes((.1, .1, .8, .2))
frame2.set_ylabel('Residui Normalizzati')
frame2.set_xlabel('L [a.u.]', fontsize=15)

ff = ChiRes(MX, dMX, N, pars, FX, 1)
x1 = np.linspace(np.min(N), np.max(N), 1000)
plt.plot(x1, 0*x1, color='red', linestyle='--', alpha=0.5)
plt.plot(N, ff, '.', color='black') 
plt.grid()

print('\n')
'''
#I seguenti dati sono presi su un reticolo 40x40 
#per beta da 0.444 fino ad 1, quindi dopo la temperatura critica
#per semplicità non è stato calcolato il cumulante di binder
'''

e40, m40, c40, x40, de40, dm40, dc40, dx40 = np.loadtxt(r'datiplot40dopopicco.dat', unpack=True)

##stima di gamma
print('stima di gamma')

#dato che prima abbiamo stimato gamma/nu con questa stima possiamo stimare sia gamma che nu
b = np.zeros(len(e40))
for i in range(len(b)):
    b[i] = 0.444 + i*(1-0.444)/len(b)

#seleziono solo un range dei dati
x = b[5:35]-bc
y = x40[5:35]
dy = dx40[5:35]
#valori che mi aspetto per i parametri ottimali
init = np.array([0.00352, -7/4, -0.013]) #aiutano la convergenza del fit

pars1, covm1 = curve_fit(FX, x, y, init, sigma=dy, absolute_sigma=False)
print(f'm = {pars1[0]:.5f} +- {np.sqrt(covm1.diagonal()[0]):.5f}')
print(f'g = {pars1[1]:.5f} +- {np.sqrt(covm1.diagonal()[1]):.5f}')
print(f'q = {pars1[2]:.5f} +- {np.sqrt(covm1.diagonal()[2]):.5f}')


chi2, Nd2 = ChiRes(y, dy, x, pars1, FX)
print(f'chi quadro = {chi2:.3f}  ({Nd2} dof)')

c2 = Corr(covm1, len(pars1))
print(c2) 


fig1 = plt.figure(8)

frame1 = fig1.add_axes((.1, .35, .8, .6))

frame1.set_title('Suscettivita reticolo 40X40 dopo il punto critico', fontsize=20)
frame1.set_ylabel('$\chi$[u.a.]', fontsize=15)
plt.grid()

plt.errorbar(x, y, dy, fmt='.', color='black', label='dati')
t = np.linspace(np.min(x), np.max(x), 10000)

plt.plot(t, FX(t, *pars1), color='blue', label='best fit') 
plt.legend(loc='best')

frame2 = fig1.add_axes((.1, .1, .8, .2))
frame2.set_ylabel('Residui Normalizzati')
frame2.set_xlabel(r'$\beta - \beta_c$ [u.a.]', fontsize=15)

ff = ChiRes(y, dy, x, pars1, FX, 1)
x1 = np.linspace(np.min(x), np.max(x), 1000)

plt.plot(x1, 0*x1, color='red', linestyle='--', alpha=0.5)
plt.plot(x, ff, '.', color='black') 
plt.grid()

print('\n')

##stima di beta
print('stima di beta')

#seleziono solo un range dei dati
x = b[:14] - bc
y = m40[:14]
dy = dm40[:14]

#valori che mi aspetto per i parametri ottimali
init = np.array([0.93, 1/8, 0.64]) #aiutano la convergenza del fit

pars2, covm2 = curve_fit(FX, x, y, init, sigma=dy, absolute_sigma=False)
print(f'm = {pars2[0]:.5f} +- {np.sqrt(covm2.diagonal()[0]):.5f}')
print(f'b = {pars2[1]:.5f} +- {np.sqrt(covm2.diagonal()[1]):.5f}')
print(f'q = {pars2[2]:.5f} +- {np.sqrt(covm2.diagonal()[2]):.5f}')


chi3, Nd3 = ChiRes(y, dy, x, pars2, FX)
print(f'chi quadro = {chi3:.3f}  ({Nd3} dof)')

c3 = Corr(covm2, len(pars2))
print(c3)


#Grafichiamo il risultato
fig1 = plt.figure(9)

frame1 = fig1.add_axes((.1, .35, .8, .6))

frame1.set_title('Magnetizzazione reticolo 40X40 dopo il punto critico', fontsize=20)
frame1.set_ylabel('M [u.a.]', fontsize=15)
plt.grid()


plt.errorbar(x, y, dy, fmt='.', color='black', label='dati')
t = np.linspace(np.min(x), np.max(x), 10000)

plt.plot(t, FX(t, *pars2), color='blue', label='best fit') 
plt.legend(loc='best')


frame2 = fig1.add_axes((.1, .1, .8, .2))
frame2.set_ylabel('Residui Normalizzati')
frame2.set_xlabel(r'$\beta - \beta_c$ [u.a.]', fontsize=15)

ff = ChiRes(y, dy, x, pars2, FX, 1)
x1 = np.linspace(np.min(x), np.max(x), 1000)

plt.plot(x1, 0*x1, color='red', linestyle='--', alpha=0.5)
plt.plot(x, ff, '.', color='black') 
plt.grid()


##stima di alpha
print('\n')
print('stima di alpha')

def Fa(x, m, g):
    '''funzione modello per il calore specifico
    '''
    return m*x**g 
    
CM, dCM = MS(C15, C20, C25, C30, C35, C40, C45, C50, dC15, dC20, dC25, dC30, dC35, dC40, dC45, dC50)



#valori che mi aspetto per i parametri ottimali
init = np.array([0.1, 0]) #aiutano la convergenza del fit

pars3, covm3 = curve_fit(Fa, N, CM, sigma=dCM, absolute_sigma=False)
print(f'm = {pars3[0]:.5f} +- {np.sqrt(covm3.diagonal()[0]):.5f}')
print(f'a = {pars3[1]:.5f} +- {np.sqrt(covm3.diagonal()[1]):.5f}')



chi4, Nd4 = ChiRes(CM, dCM, N, pars3, Fa)
print(f'chi quadro = {chi4:.3f}  ({Nd4} dof)')

c4 = Corr(covm3, len(pars3))
print(c4) 

fig1 = plt.figure(10)

frame1 = fig1.add_axes((.1, .35, .8, .6))
#frame1=fig1.add_axes((trasla lateralmente, trasla verticamente, larghezza, altezza))
frame1.set_title('Massimo del calore specifico al variare di L', fontsize=20)
frame1.set_ylabel('Calore specifico [a.u.]', fontsize=15)
plt.grid()

plt.errorbar(N, CM, dCM, fmt='.', color='black', label='dati') 
t=np.linspace(np.min(N), np.max(N), 10000)

plt.plot(t, Fa(t, *pars3), color='blue', label='best fit') 
plt.legend(loc='best')


frame2 = fig1.add_axes((.1, .1, .8, .2))
frame2.set_ylabel('Residui Normalizzati')
frame2.set_xlabel(r'$\beta - \beta_c$ [a.u.]', fontsize=15)

ff = ChiRes(CM, dCM, N, pars3, Fa, 1)
x1 = np.linspace(np.min(N), np.max(N), 1000)

plt.plot(x1, 0*x1, color='red', linestyle='--', alpha=0.5)
plt.plot(N, ff, '.', color='black') 
plt.grid()

#caclolo deigli indici critici ed associati errori
b = pars2[1]
db = np.sqrt(covm2.diagonal()[1])

g = -pars1[1]
dg = np.sqrt(covm1.diagonal()[1])

n = -pars1[1]/pars[1]
dn = np.sqrt((np.sqrt(covm.diagonal()[1])/pars[1])**2 + (np.sqrt(covm1.diagonal()[1])/pars1[1])**2)*n

a = pars3[1]
da = np.sqrt(covm3.diagonal()[1])


print('\n')
print(f"gamma = {g:.5f} +- {dg:.5f}")
print(f"beta  = {b:.5f} +- {db:.5f}")
print(f"nu    = {n:.5f} +- {dn:.5f}")
print(f"alpha = {a:.5f} +- {da:.5f}")

print('dg/nu = %.5f' %(pars[1]-7/4))
print('db = %.5f' %(pars2[1]-1/8))
print('dg = %.5f' %(pars1[1]+7/4))
print('da = %.5f' %(pars3[1]-0))

##grafici del finite size scaling

plt.figure(11)
plt.title("Finite size scaling della magnetizzazione", fontsize=20)
plt.xlabel(r"$(\beta-\beta_c)L^{1/ \nu}$", fontsize=15)
plt.ylabel(r"$|M|/L^{-b/ \nu}$", fontsize=15)
plt.grid()
plt.errorbar((B-bc)*20**(1/n), M20/(20**(-b/n)), fmt='.', color='black', label='L=20')
plt.errorbar((B-bc)*30**(1/n), M30/(30**(-b/n)), fmt='v', color='black', label='L=30')
plt.errorbar((B-bc)*40**(1/n), M40/(40**(-b/n)), fmt='*', color='black', label='L=40')
plt.errorbar((B-bc)*50**(1/n), M50/(50**(-b/n)), fmt='^', color='black', label='L=50')
plt.legend(loc='best')

plt.figure(12)
plt.title("Finite size scaling della suscettività", fontsize=20)
plt.xlabel(r"$(\beta-\beta_c)L^{1/ \nu}$", fontsize=15)
plt.ylabel(r"$ \chi /L^{\gamma/ \nu}$", fontsize=15)
plt.grid()
plt.errorbar((B-bc)*20**(1/n), X20/(20**(g/n)), fmt='.', color='black', label='L=20')
plt.errorbar((B-bc)*30**(1/n), X30/(30**(g/n)), fmt='v', color='black', label='L=30')
plt.errorbar((B-bc)*40**(1/n), X40/(40**(g/n)), fmt='*', color='black', label='L=40')
plt.errorbar((B-bc)*50**(1/n), X50/(50**(g/n)), fmt='^', color='black', label='L=50')

plt.legend(loc='best')

plt.show()
