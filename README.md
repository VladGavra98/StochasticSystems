# Stochastic Systems - Spectral Anlysiss & Monte Carlo Simulation


Contains the source code used to analyse the Cessna Citation I response to stochastic inputs.
Parts 1,2,3 and completed in Python:
 - check the windex for choosing the turbulence condition 
  - change the plotting flag for switching between time and frequncy plots

Part 4 (variances) was implemented in Python and MTALAB (only the impulse response). The impulse response method showed bugs in Python due to way lsim works. 
Both frameworks showed simialr resutls should the observation time T be hhigh enough. The variances were calaculted multiple times to ensure consistency.


Main files: simulation.py, H_mehtod.m
Aircraft model class: Cessna_model.py
Testing & verification: testing.py (NOT TO BE CONSIDERED)

 Other linear systems are also considered and/or tested. 
