import numpy as np
import matplotlib.pyplot as plt

import grafici

#il nuero finale identifica la grandezza del reticolo
E15, M15, C15, X15, cb15, dE15, dM15, dC15, dX15, dcb15 = np.loadtxt(r'datiplot/15.dat', unpack=True)
E20, M20, C20, X20, cb20, dE20, dM20, dC20, dX20, dcb20 = np.loadtxt(r'datiplot/20.dat', unpack=True)
E25, M25, C25, X25, cb25, dE25, dM25, dC25, dX25, dcb25 = np.loadtxt(r'datiplot/25.dat', unpack=True)
E30, M30, C30, X30, cb30, dE30, dM30, dC30, dX30, dcb30 = np.loadtxt(r'datiplot/30.dat', unpack=True)
E35, M35, C35, X35, cb35, dE35, dM35, dC35, dX35, dcb35 = np.loadtxt(r'datiplot/35.dat', unpack=True)
E40, M40, C40, X40, cb40, dE40, dM40, dC40, dX40, dcb40 = np.loadtxt(r'datiplot/40.dat', unpack=True)
E45, M45, C45, X45, cb45, dE45, dM45, dC45, dX45, dcb45 = np.loadtxt(r'datiplot/45.dat', unpack=True)
E50, M50, C50, X50, cb50, dE50, dM50, dC50, dX50, dcb50 = np.loadtxt(r'datiplot/50.dat', unpack=True)


H = 0
#valori esatti per gli indici e la temperatura critica
#n = 1
#b = 1/8
#g = 7/4
#a = 0
#bc = 0.4406868

##Alcuni plot a titolo espositivo

param = np.loadtxt('input.txt', unpack=True)
b_min = param[3]
b_max = param[4]
N = int(param[5])

B = np.linspace(b_min, b_max, N)


Title = f'Simulazione del modello di Ising 2D tramite Metropolis \n Campo magnetico esterno B={H}'
xlabel = r'$\beta$ [u.a.]'   
grafici.plot(B, [E20, E30, E40, E50], 
                [dE20, dE30, dE40, dE50],
                [20, 30, 40, 50], 1, Title, xlabel, "Energia [u.a.]")
                
grafici.plot(B, [M20, M30, M40, M50],
                [dM20, dM30, dM40, dM50],
                [20, 30, 40, 50], 2, Title, xlabel, "Magetizzazione [u.a.]")
                
grafici.plot(B, [C20, C30, C40, C50], 
                [dC20, dC30, dC40, dC50],
                [20, 30, 40, 50], 3, Title, xlabel, "Calore specifico [u.a.]")
                 
grafici.plot(B, [X20, X30, X40, X50], 
                [dX20, dX30, dX40, dX50], 
                [20, 30, 40, 50], 4, Title, xlabel, "Suscettività [u.a.]") 
 
grafici.plotbinder(B, [cb20, cb30, cb40, cb50], 
                      [dcb20, dcb30, dcb40, dcb50],
                      [20, 30, 40, 50], 5, Title, xlabel, "Cumulante di binder") 

##plot del cumulante di binder da cui si può stimare il beta critico


#valore di beta ricavato
bc = 0.440
dbc = 0.003
print(f"beta critico = {bc:.3f} +- {dbc:.3f}")

##grafico del comportamento della magnetizzazione nelle due fais della transizione

mag, ene = np.loadtxt(r'datiplot/20hist.dat', skiprows=3, unpack=True)
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


def F(x, m, g, q):
    '''funzione modello per i fit successivi
    '''
    return m*x**g + q

##stima di gamma/nu
print('stima di gamma/nu')

#funzione per estrarre il massimo della suscettività e relativo errore
def MS(*arg):
    '''
    *arg permette alla funzione di prendere in input infiniti argomenti.
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
Title = 'Massimo della suscettività al variare di L'
ylabel = r'$\chi_{max} [a.u.]$'
xlabel = r'L [a.u.]$'
init = np.array([1, 1, 1])
  
pars, dpars = grafici.fit(F, N, MX, dMX, init, 7, Title, xlabel, ylabel)

print('\n')
'''
#I seguenti dati sono presi su un reticolo 40x40 
#per beta da 0.444 fino ad 1, quindi dopo la temperatura critica
#per semplicità non è stato calcolato il cumulante di binder
'''

e40, m40, c40, x40, de40, dm40, dc40, dx40 = np.loadtxt(r'datiplot/40dopopicco.dat', unpack=True)

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

Title = 'Suscettivita reticolo 40X40 dopo il punto critico'
ylabel = r'$chi [a.u.]$'
xlabel = r'$\beta - \beta_c$ [u.a.]$'

pars1, dpars1 = grafici.fit(F, x, y, dy, init, 8, Title, xlabel, ylabel)

print('\n')

##stima di beta
print('stima di beta')

#seleziono solo un range dei dati
x = b[:14] - bc
y = m40[:14]
dy = dm40[:14]

#valori che mi aspetto per i parametri ottimali
init = np.array([0.93, 1/8, 0.64]) #aiutano la convergenza del fit

Title = 'Magnetizzazione reticolo 40X40 dopo il punto critico'
ylabel = r'M [a.u.]'
xlabel = r'$\beta - \beta_c$ [u.a.]$'

pars2, dpars2 = grafici.fit(F, x, y, dy, init, 9, Title, xlabel, ylabel)


##stima di alpha
print('\n')
print('stima di alpha')

def Fa(x, m, g):
    '''funzione modello per i fit successivi
    '''
    return m*x**g 
    
CM, dCM = MS(C15, C20, C25, C30, C35, C40, C45, C50, dC15, dC20, dC25, dC30, dC35, dC40, dC45, dC50)
N = np.arange(15, 51, 5)


#valori che mi aspetto per i parametri ottimali
init = np.array([0.1, 0]) #aiutano la convergenza del fit

Title = 'Massimo del calore specifico al variare di L'
ylabel = r'Calore specifico [a.u.]'
xlabel = r'$\beta - \beta_c$ [u.a.]$'

pars3, dpars3 = grafici.fit(Fa, N, CM, dCM, init, 10, Title, xlabel, ylabel)

#caclolo deigli indici critici ed associati errori


g = -pars1[1]
dg = dpars[1]

n = -pars1[1]/pars[1]
dn = np.sqrt((dpars[1]/pars[1])**2 + (dpars[1]/pars1[1])**2)*n

b = pars2[1]
db = dpars2[1]

a = pars3[1]
da = dpars3[1]


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

Title = 'Finite size scaling della magnetizzazione'
xlabel = r'$(\beta-\beta_c)L^{1/ \nu}$'
ylabel = r'$|M|/L^{-b/ \nu}$'
    
grafici.FSS(B, 1/n, -b/n, bc, 
            [M20, M30, M40, M50],
            [dM20, dM30, dM40, dM50], 
            [20, 30, 40, 50], 11, Title, xlabel, ylabel)

Title = 'Finite size scaling della suscettività'
xlabel = r'$(\beta-\beta_c)L^{1/ \nu}$'
ylabel = r'$ \chi /L^{\gamma/ \nu}$'

grafici.FSS(B, 1/n, g/n, bc, 
            [X20, X30, X40, X50],
            [dX20, dX30, dX40, dX50],
            [20, 30, 40, 50], 12, Title, xlabel, ylabel) 


##estrapolazione di beta continuo
print('\n')
print('stima di beta critico')

def taglio(y, z, b, q):
    '''
    funzione per selezionare un certo range di dati
    '''
    Y = np.zeros((len(y), int(2*q)))
    Z = np.zeros((len(z), int(2*q)))
    B = np.zeros((len(z), int(2*q)))
    
    for j in range(len(y)):
        m = np.max(y[j])
        i = int(np.where(y[j] == m)[0][0]) #conservo il rispettivo indice
        Y[j,:] = y[j][i-q:i+q]
        Z[j,:] = z[j][i-q:i+q]
        B[j,:] = b[i-q:i+q]
    return Y, Z, B
    

X, dX, B = taglio([X15, X20, X25, X30, X35, X40, X45, X50], 
           [dX15, dX20, dX25, dX30, dX35, dX40, dX45, dX50], B, 5)
           
L = np.array([15, 20, 25, 30, 35, 40, 45, 50])
def G(x, a, b, c):
    return a - b*(x - c)**2

Pars =  np.zeros((len(X[:,0]), 3))
dPars =  np.zeros((len(X[:,0]), 3))

print('fit massimo della suscettività')

for x, y, dy, j, l in zip(B, X, dX, range(len(X[:,0])), L):
    
    init = np.array([np.max(y), l ,x[5]])
    Pars[j,:], dPars[j,:] = grafici.fit(G, x, y, dy, init, 0, Title, xlabel, ylabel, False)

Title = f'Suscettività intorno al punto critico'
xlabel = r'$\beta$ [u.a.]'
ylabel = '$\chi [u.a.]$'

grafici.plotfit(G, Pars, B, X, dX, L, 13, Title, xlabel, ylabel)

print('\n')

def Bl(x, m, q):
    return x*m + q
    
betacl = Pars[:, 2]
dbetacl = dPars[:, 2]
init=np.array([1, 0.4406868])
 
Title = 'Estrapolazione di beta al continuo'
xlabel = r'1/L'
ylabel = r'$ \beta_{c}(L)$'  
pars4, dpars4 = grafici.fit(Bl, 1/L, betacl, dbetacl, init, 14, Title, xlabel, ylabel)


print('\n')
print(f'beta critico: {pars4[1]:.5f} -+ {dpars4[1]:.5f}')

plt.show()
