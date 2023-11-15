# vitisfft
Template for implementing array-interfaced FFTs with Vitis HLS.

The repository contains the C++ source files needed to use the AMD FFT IP for both forward and inverse complex and real FFTs.

The testbench in [./vitisfft_t.cpp](vitisfft_t.cpp) processes some example signals and, optionally, writes the results thereof to three files that the associated Jupyter notebook can read and analyze.
