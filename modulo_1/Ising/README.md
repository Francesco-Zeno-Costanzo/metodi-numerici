In this folder there are the code to simulate and analize data for ising's model. 
The code that performs the simulation is ising2d.f; term.f calculates the thermodynamic quantities and relative errors.

The ising model is an attempt to simulate the structure of a physical ferromagnetic substance, or more accurately, to simulate a domain in a ferromagnetic substance (or anti-ferromagnetic).
We will consider a 2-dimensional periodic lattice. Associated with each lattice site is a spin variable which is a number that is eitheir +1 or -1.
The hamiltonian of the system is:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{H}&space;=&space;-J\sum_{<ij>}&space;\sigma_i&space;\sigma_j&space;-&space;B\sum_{i=0}^N&space;\sigma_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{H}&space;=&space;-J\sum_{<ij>}&space;\sigma_i&space;\sigma_j&space;-&space;B\sum_{i=0}^N&space;\sigma_i" title="\mathcal{H} = -J\sum_{<ij>} \sigma_i \sigma_j - B\sum_{i=0}^N \sigma_i" /></a>

J>0 corresponds to the ferromagnetic case, J<0 anti-ferromagnetic case. From now on, J = 1 and B = 0 will be assumed.

The thermodynamic quantities are calculated in the following way:

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;E&=\frac{<\varepsilon>}{V}&space;\\&space;C&=\frac{(<\varepsilon^2>-<\varepsilon>^2)}{V}\\&space;M&=\frac{<m>}{V}\\&space;\chi&=\frac{(<m^2>-<m>^2)}{V}&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{align*}&space;E&=\frac{<\varepsilon>}{V}&space;\\&space;C&=\frac{(<\varepsilon^2>-<\varepsilon>^2)}{V}\\&space;M&=\frac{<m>}{V}\\&space;\chi&=\frac{(<m^2>-<m>^2)}{V}&space;\end{align*}" title="\begin{align*} E&=\frac{<\varepsilon>}{V} \\ C&=\frac{(<\varepsilon^2>-<\varepsilon>^2)}{V}\\ M&=\frac{<m>}{V}\\ \chi&=\frac{(<m^2>-<m>^2)}{V} \end{align*}" /></a>

where m is the sum of the spins and epsilon the energy given by the Hamiltonian:

<a href="https://www.codecogs.com/eqnedit.php?latex=m=\sum_i&space;\sigma_i&space;\hspace{10&space;mm}&space;\varepsilon&space;=&space;\frac{\mathcal{H}}{2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?m=\sum_i&space;\sigma_i&space;\hspace{10&space;mm}&space;\varepsilon&space;=&space;\frac{\mathcal{H}}{2}" title="m=\sum_i \sigma_i \hspace{10 mm} \varepsilon = \frac{\mathcal{H}}{2}" /></a>

another important value that will be calculated is the binder cumulant, the kurtosis of the order parameter, defined as:

<a href="https://www.codecogs.com/eqnedit.php?latex=U=\frac{<m^4>}{<m^2>^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U=\frac{<m^4>}{<m^2>^2}" title="U=\frac{<m^4>}{<m^2>^2}" /></a>

More information can be found in the pdf of the final report. Also in the montecarlo repository you can find a code for the simulation of ising written in C.
