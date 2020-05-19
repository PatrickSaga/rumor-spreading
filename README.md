# Rumor-spreading
Fortran & Python codes to simulate the [spreading of rumors](https://en.wikipedia.org/wiki/Rumor_spread_in_social_network).
This code was developed for the [Stochastic simulation methods](https://ifisc.uib-csic.es/master/programme-syllabus/stochastic-simulation-methods/) subject imparted in the [Complex Systems master course](https://ifisc.uib-csic.es/master/) from the [IFISC](https://ifisc.uib-csic.es/es/).

The simulation is comprised by a Gillespie algorithm for a SIHR model (susceptible-infectious-hibernator-recovered), as introduced in [Zhao et al (2012)](https://www.sciencedirect.com/science/article/pii/S0378437111009058). 

## Code files
- **dranxor.f90**. Random number generator in Fortran90 from "Generation of Gaussian distributed random numbers by using a numerical inversion method", R. Toral, A. Chakrabarti. Computer Physics Communications, 74 (1993) 327-334. Please give credit to the authors of the code R. Toral, A. Chakrabarti, and properly cite them when using or disseminating it. This repository includes it under the free use/dissemination statement of the authors. 
- **rumor.f**. Fortran code for the SIHR dynamics with a Gillespie algorithm
- **rumor_spreading.ipynb**. Python code (IPython) for the SIHR dynamics with a Gillespie algorithm
- **gillespie_plots**. Python script (IPython) to generate the output plots

## Contributors
Patrick Sánchez Galea
