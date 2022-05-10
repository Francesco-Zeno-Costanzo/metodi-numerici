import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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


##====================================================================================

## LIMITE AL CONTINUO
## OVER-RELAXATION
nt, eco, deco = np.loadtxt('densmo.dat', unpack=True)

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

nt, eno, deno = np.loadtxt('densm.dat', unpack=True)

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

plt.show()
