# SFA OAM Linear
Using the SFA to compute the orbital angular momentum of a photoelectron ionized via a strong linear sin^2 pulse laser field. This code is an edited verions of the circular code based on the results of this [preprint arXiv paper](https://arxiv.org/abs/2102.07453)

To build on linux use the build.sh script, run:
`./build.sh` in the directory

If not on linux run the python commands inside the build script.

You must have python, cython and a c++ complier (e.g. g++) installed for this to work.
To install cython check out [anaconda documentation](https://anaconda.org/anaconda/cython) and [cython documentation](https://cython.readthedocs.io/en/latest/src/quickstart/install.html)


The data processing and plotting can be found in the python notebook OAM_Linear_Pulse.ipynb. It is recommend to use jupyter lab and anaconda to view this. There are also the scripts "SFALinear_GenData.py" to generate data on a server and "SFALinear_PlotData.py" to plot this data locally.

## SFALinearPulse.pyx information
A '.pyx' file is a cython file that uses the flexibility of python but with parts that have the speed of c++.

The SFALinearPulse.pyx file codes a python class called SFALinearPulse. As the name suggests it implements the SFA for a linear polarized laser field with a sin^2 pulse.

## Publicaton Specific information
### Ultrafast imaging of molecular chirality with photoelectron vortices
Available at https://arxiv.org/abs/2202.07289
The generation and plotting scripts: PublicationScripts/Gen_PEVD_SFA_Data.py and PublicationScripts/load_and_plot_SFA.py

# Information on the Code
## Class member functions
The only accessible member functions from outside the class are those with cpdef  or def in front of them. Those with only cdef are internal and can not be called outside the class. The most important functions are:

- Af(t) vector potential in z direction
- Ef(t) electric field in z direction
- S(p, theta, phi, t) semi-classical action
- TimesGen(p, theta, phi) find all times of ionization via the saddle point approximation for a specific momentum coordinate
- M(p, theta, phi, tf) Transition amplitude for a final momentum point computed using the saddle point approximaton
- Spectra(s, E, phi,  t , err, limit) Compute the Spectra, integrating over theta, results should be invarient wrt phi, hence its inclusion
- Ml(p, theta, Nphi) output vector of size Nphi of OAM dependent transition amplitude using FFT of M(p, theta, phi, tf)
- d0(p, theta, phi, t) function for the matrix element incorporating the effect of the bound state

### Class member variables
There are the following member variables relevant to the input parameters of the problem and as ways to control the class.

The following are class member that are set when you make an instance and remain fixed after that point (readonly)

- Ip: the ionization potential
- Up: Poneromotive energy
- omega: carrier frequency
- N: number of lasser cycles
- CEP: carrier envelope phase
- Target: determines,which prefactor is chosen
- AlignmentAngle: can rotate some initial states in prefactor d0


Other member varibles include
constant complex iminary unit I, constant Pi and the square root of Pi
