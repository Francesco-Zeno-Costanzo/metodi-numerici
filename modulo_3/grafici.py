import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

    
def fit(F, x, y, dy, init, k, Title, xlabel, ylabel, plot=True):
    '''
    Funzione che esegui il fi e fa il plot associato
    F è la funzione di fit, x l'array dei dati sulle x
    y l'array di dati sulle y e dy il relativo errore
    init è l'array dei parametri inizali e k il numero
    della figura ceh può non essere visualizzata se
    plot assume valore False, di default è True
    '''
    pars, covm = curve_fit(F, x, y, init, sigma=dy)
    err = np.sqrt(covm.diagonal())
    for p, dp, i in zip(pars, err, range(len(pars))):
        print(f"pars{i} = {p:.5f} +- {dp:.5f}")
    

    chisq = sum(((y - F(x, *pars))/dy)**2.)
    ndof = len(y) - len(pars)
    print('chi quadro = %.3f (%d dof)' % (chisq, ndof))


    c = np.zeros((len(pars), len(pars)))

    for i in range(0, len(pars)):
        for j in range(0, len(pars)):
            c[i][j]=(covm[i][j])/(np.sqrt(covm.diagonal()[i])*np.sqrt(covm.diagonal()[j]))
    print(c) #matrice di correlazione

    if plot :
        fig = plt.figure(k)
        frame1 = fig.add_axes((.1,.35,.8,.6))
        frame1.set_title(Title, fontsize=20)
        frame1.set_ylabel(ylabel, fontsize=15)
        frame1.grid()

        frame1.errorbar(x, y, dy, fmt='.', color='black', label='dati')
        t = np.linspace(np.min(x), np.max(x), 1000)
        frame1.plot(t, F(t, *pars), color='blue', label='best fit')   
        frame1.legend(loc='best')
        frame2 = fig.add_axes((.1,.1,.8,.2))

        ff = (y - F(x, *pars))/dy
        frame2.set_ylabel('Residui Normalizzati')
        frame2.set_xlabel(xlabel, fontsize=15)

        frame2.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
        frame2.plot(x, ff, '.', color='black') #grafico i residui normalizzati
        frame2.grid()
    
    return pars, err
    
    
def plotfit(F, Pars, X, args, d_args, L, k, title, xlabel, ylabel):
    '''
    funzione per plottare più fit
    che devono già essere stati esguiti
    su un unico grafico
    '''
    fig = plt.figure(k)
    
    frame1 = fig.add_axes((.1,.35,.8,.6))
    frame1.set_title(title, fontsize=20)
    frame1.set_ylabel(ylabel, fontsize=15)
    frame1.grid()
    
    frame2 = fig.add_axes((.1,.1,.8,.2))
    frame2.set_ylabel('Residui Normalizzati')
    frame2.set_xlabel(xlabel, fontsize=15)
    frame2.grid()
    
    colors = plt.cm.jet(np.linspace(0, 1, len(args)))
    
    for x, y, dy, pars, i in zip(X, args, d_args, Pars, range(len(args))):
      
        frame1.errorbar(x, y, dy, fmt='.', color=colors[i], label=fr'dati per $\beta \omega$={L[i]}')
        t = np.linspace(np.min(x), np.max(x), 1000)
        frame1.plot(t, F(t, *pars), color=colors[i])   
        frame1.legend(loc='best')
        

        ff = (y - F(x, *pars))/dy
        

        frame2.plot(t, 0*t, color='red', linestyle='--', alpha=0.5) #grafico la retta costantemente zero
        frame2.plot(x, ff, '.', color=colors[i]) #grafico i residui normalizzati
        

