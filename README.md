# QuantumSimulator
Homeworks/ Project for the course quantum informations and computing 2023/2024 written in C++.
It can be built in two ways:
 - `make` that will generate the executable `quantumSimulator` (which at the moment simulates a time dependent quantum harmonic oscillator, i.e. 5th homework)
 - `make test` that will generate the executable `test` which runs BOOST unit tests.

To build the following dependencies are required
  - libopenblas-dev (for doing matrix operations multithread)
  - liblapack-dev (for LAPACK functions)
  - libboost-all-dev (only for unit tests)
  - fftw configured with multithreading flags, see https://www.fftw.org/fftw3_doc/Multi_002dthreaded-FFTW.html (to perform fast fourier transform)

