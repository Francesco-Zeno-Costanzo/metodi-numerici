This code are some exampe of distribution sampling with metropolis and accept and reject method.
accrig.py use the last method to sample

<a href="https://www.codecogs.com/eqnedit.php?latex=\sqrt{1-x^2}e^{-\gamma&space;x}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sqrt{1-x^2}e^{-\gamma&space;x}" title="\sqrt{1-x^2}e^{-\gamma x}" /></a>

metrogauss.f use metropolis to sample a gaussian distributione with average = 5 and variance = 1
The Metropolis Algorithm, on the other hand, is a monte carlo method for sampling an unknown distribution starting from one proportional to it, p(x).
We have to generate a sequence of extractions x0, . . . xN , starting from a given x0. We generate xx uniformly in the interval [xk - d, xk + d];
if the probability ratio between the new variable and the old one is greater than a randomly generated number in [0, 1] the move is accepted, otherwise rejected.

<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;xx&space;\in&space;[x_k&space;-&space;\delta,&space;x_k&space;&plus;&space;\delta&space;]\\&space;r&space;\in&space;[0,&space;1]\\&space;\frac{p(xx)}{p(x_k)}&space;<&space;r&space;\rightarrow&space;x_{k&plus;1}=xx\\&space;\frac{p(xx)}{p(x_k)}&space;>&space;r&space;\rightarrow&space;x_{k&plus;1}=x_k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;xx&space;\in&space;[x_k&space;-&space;\delta,&space;x_k&space;&plus;&space;\delta&space;]\\&space;r&space;\in&space;[0,&space;1]\\&space;\frac{p(xx)}{p(x_k)}&space;<&space;r&space;\rightarrow&space;x_{k&plus;1}=xx\\&space;\frac{p(xx)}{p(x_k)}&space;>&space;r&space;\rightarrow&space;x_{k&plus;1}=x_k" title="\\ xx \in [x_k - \delta, x_k + \delta ]\\ r \in [0, 1]\\ \frac{p(xx)}{p(x_k)} < r \rightarrow x_{k+1}=xx\\ \frac{p(xx)}{p(x_k)} > r \rightarrow x_{k+1}=x_k" /></a>

block.f use the data-blocking technique to compute the error of metropolis sample because the data are correlated
bootstrap.f is binned bootstrap to do the same thing
corr.f compute the autocorrelation function of metropolis sample
