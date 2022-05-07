The purpose of this module is to study the thermodynamics of the quantum system via path integral. In fact, it is possible to write the average value of an operator in terms of path integral with due care and therefore adopt montecarlo techniques

<img src="http://latex.codecogs.com/svg.latex?Z&space;=&space;\mathcal{N}&space;\int&space;\mathcal{D}x(\tau)&space;\exp\Biggl(&space;\frac{-S_E[x(\tau)]}{\hbar}&space;\Biggr)&space;\hspace{10&space;mm}&space;<O>&space;=&space;\frac{\int&space;\mathcal{D}x(\tau)&space;\exp\Biggl(&space;\frac{-S_E[x(\tau)]}{\hbar}&space;\Biggr)&space;O[x(\tau)]}{\int&space;\mathcal{D}x(\tau)&space;\exp\Biggl(&space;\frac{-S_E[x(\tau)]}{\hbar}&space;\Biggr)}" title="http://latex.codecogs.com/svg.latex?Z = \mathcal{N} \int \mathcal{D}x(\tau) \exp\Biggl( \frac{-S_E[x(\tau)]}{\hbar} \Biggr) \hspace{10 mm} <O> = \frac{\int \mathcal{D}x(\tau) \exp\Biggl( \frac{-S_E[x(\tau)]}{\hbar} \Biggr) O[x(\tau)]}{\int \mathcal{D}x(\tau) \exp\Biggl( \frac{-S_E[x(\tau)]}{\hbar} \Biggr)}" />

The expression of the mean value of O has a probabilistic interpretation that will be the basis of the simulation using Monte Carlo techniques; we have the mean value of a function of the path evaluated on a probability distribution for the paths which is:

<img src="http://latex.codecogs.com/svg.latex?P[x(\tau)]&space;=&space;\frac{\exp\Biggl(&space;\frac{-S_E[x(\tau)]}{\hbar}&space;\Biggr)}{\int&space;\mathcal{D}x(\tau)&space;\exp\Biggl(&space;\frac{-S_E[x(\tau)]}{\hbar}&space;\Biggr)}" title="http://latex.codecogs.com/svg.latex?P[x(\tau)] = \frac{\exp\Biggl( \frac{-S_E[x(\tau)]}{\hbar} \Biggr)}{\int \mathcal{D}x(\tau) \exp\Biggl( \frac{-S_E[x(\tau)]}{\hbar} \Biggr)}" />

The simulated system is the harmonic oscillator.

The oscarm.f code is the one that simulates the system and makes the necessary measurements, which are analyzed in the other two fortran codes.

In enrgia.py are show the results for the internal energy of sistem.
grafici.py is a code that energia.py uses to plot the results so that the latter is more ordered and shorter.
correlazioni.py compute the the splitting of energy level E_1 - E_0 and E_2 - E_0


