# vitisfft

Template for implementing array-interfaced FFTs with Vitis HLS. The repository contains the C++ source files needed to use the AMD FFT IP for both forward and inverse complex FFTs. The cores use streams internally as the pure array-based solution does not implement and/or requires twice the number of memories as the stream-based solution (tested with Vitis HLS 2023.2).

The testbench in [./vitisfft_t.cpp](vitisfft_t.cpp) processes some example signals and, optionally, writes the results thereof to three files that the associated Jupyter notebook can read and analyze.
