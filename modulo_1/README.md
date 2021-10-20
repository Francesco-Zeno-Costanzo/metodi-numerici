These codes are an example of distribution sampling

accrig.py uses accept-reject method to sample:

<a href="https://www.codecogs.com/eqnedit.php?latex=\sqrt{1-x^2}e^{-\gamma&space;x}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sqrt{1-x^2}e^{-\gamma&space;x}" title="\sqrt{1-x^2}e^{-\gamma x}" /></a>


The accept-reject method allows to sample a distribution p (x) with a known g (x), such that c * g (x)> p (x) for all x, where c is a constant.
We extract a variable x distributed according to g (x), then a second variable y between [0, c * g (x)] uniformly distributed; if y <p (x) then I accept the variable x. Finally the new set of accepted variables x will be distributed according to p (x)


metrogauss.f uses metropolis algorithm to sample a gaussian distribution with average = 5 and variance = 1


The Metropolis Algorithm, is a Monte Carlo method for sampling an unknown distribution starting from one proportional to it, p(x).
We have to generate a sequence of draws x0, . . . xN , starting from a given x0. We generate xx uniformly in the interval [xk - d, xk + d];
if the probability ratio between the new variable and the old one is greater than a randomly generated number in [0, 1] the move is accepted, otherwise rejected.

<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;xx&space;\in&space;[x_k&space;-&space;\delta,&space;x_k&space;&plus;&space;\delta&space;]\\&space;r&space;\in&space;[0,&space;1]\\&space;\frac{p(xx)}{p(x_k)}&space;<&space;r&space;\rightarrow&space;x_{k&plus;1}=xx\\&space;\frac{p(xx)}{p(x_k)}&space;>&space;r&space;\rightarrow&space;x_{k&plus;1}=x_k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;xx&space;\in&space;[x_k&space;-&space;\delta,&space;x_k&space;&plus;&space;\delta&space;]\\&space;r&space;\in&space;[0,&space;1]\\&space;\frac{p(xx)}{p(x_k)}&space;<&space;r&space;\rightarrow&space;x_{k&plus;1}=xx\\&space;\frac{p(xx)}{p(x_k)}&space;>&space;r&space;\rightarrow&space;x_{k&plus;1}=x_k" title="\\ xx \in [x_k - \delta, x_k + \delta ]\\ r \in [0, 1]\\ \frac{p(xx)}{p(x_k)} < r \rightarrow x_{k+1}=xx\\ \frac{p(xx)}{p(x_k)} > r \rightarrow x_{k+1}=x_k" /></a>

block.f uses the data-blocking technique to compute the error of the metropolis sample because the data are correlated


bootstrap.f is binned bootstrap and it is used to do the same thing


corr.f computes the autocorrelation function of the metropolis sample

plotmetro.py plots the data returned by the previous codes


boxmuller.py uses boxmuller algorithm to sample the same gaussian distribution, but in this case the error is computed with the classical standard deviation and with the classical bootstrap.
The boxmuller algorithm allows to sample the gaussian distubution by making variable changes:


<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;P(y)&space;=&space;e^{-y^2}&space;dy&space;\rightarrow&space;P(y_1,&space;y_2)&space;=&space;e^{-y_1^2}&space;e^{-y_2^2}&space;dy_1&space;dy_2\hspace{5&space;mm}&space;\text{variable&space;exchange}\\&space;y_1&space;=&space;\rho&space;\cos(\vartheta)\\&space;y_2&space;=&space;\rho&space;\sin(\vartheta)\\&space;P(\rho,&space;\vartheta)&space;=&space;\rho&space;e^{-\rho^2}&space;d\rho&space;d\vartheta\hspace{10&space;mm}&space;\text{variable&space;change&space;again}\\&space;z=\rho^2&space;\rightarrow&space;P(z,&space;\vartheta)=\frac{1}{2}e^{-z}&space;dz&space;d\vartheta\\&space;x_1,&space;x_2&space;\in&space;[0,&space;1)&space;\rightarrow&space;\vartheta=2\pi&space;x_1,&space;z=-\ln(1-x_2)&space;\rightarrow\\&space;\rightarrow&space;y_1=\sqrt{z}\cos(\vartheta),&space;y_2=\sqrt{z}\sin(\vartheta)&space;\\&space;\text{now&space;}&space;y_1&space;\text{&space;and&space;}&space;y_2&space;\text{&space;are&space;gaussian&space;distributed}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;P(y)&space;=&space;e^{-y^2}&space;dy&space;\rightarrow&space;P(y_1,&space;y_2)&space;=&space;e^{-y_1^2}&space;e^{-y_2^2}&space;dy_1&space;dy_2\hspace{5&space;mm}&space;\text{variable&space;exchange}\\&space;y_1&space;=&space;\rho&space;\cos(\vartheta)\\&space;y_2&space;=&space;\rho&space;\sin(\vartheta)\\&space;P(\rho,&space;\vartheta)&space;=&space;\rho&space;e^{-\rho^2}&space;d\rho&space;d\vartheta\hspace{10&space;mm}&space;\text{variable&space;change&space;again}\\&space;z=\rho^2&space;\rightarrow&space;P(z,&space;\vartheta)=\frac{1}{2}e^{-z}&space;dz&space;d\vartheta\\&space;x_1,&space;x_2&space;\in&space;[0,&space;1)&space;\rightarrow&space;\vartheta=2\pi&space;x_1,&space;z=-\ln(1-x_2)&space;\rightarrow\\&space;\rightarrow&space;y_1=\sqrt{z}\cos(\vartheta),&space;y_2=\sqrt{z}\sin(\vartheta)&space;\\&space;\text{now&space;}&space;y_1&space;\text{&space;and&space;}&space;y_2&space;\text{&space;are&space;gaussian&space;distributed}" title="\\ P(y) = e^{-y^2} dy \rightarrow P(y_1, y_2) = e^{-y_1^2} e^{-y_2^2} dy_1 dy_2\hspace{5 mm} \text{variable exchange}\\ y_1 = \rho \cos(\vartheta)\\ y_2 = \rho \sin(\vartheta)\\ P(\rho, \vartheta) = \rho e^{-\rho^2} d\rho d\vartheta\hspace{10 mm} \text{variable change again}\\ z=\rho^2 \rightarrow P(z, \vartheta)=\frac{1}{2}e^{-z} dz d\vartheta\\ x_1, x_2 \in [0, 1) \rightarrow \vartheta=2\pi x_1, z=-\ln(1-x_2) \rightarrow\\ \rightarrow y_1=\sqrt{z}\cos(\vartheta), y_2=\sqrt{z}\sin(\vartheta) \\ \text{now } y_1 \text{ and } y_2 \text{ are gaussian distributed}" /></a>

