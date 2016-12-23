# fourier-series
----

A simple MATLAB script that computes the Fourier series for some kinds of signals. Made originally in 2014.

It shows some basic functions of MATLAB, so may be interesting if you are learning how to use it, especially if you are an Electrical Engineering student.

Also a first test with GitHub.

This script calculates the Fourier Series aproximation of a chosen
signal:
- Rectangular pulse, with option to choose the ON/OFF ratio
- Exponetial function, with option to choose the value for the constant
- Sawtooth function
- Triangular function

After choosing the function type, you choose how many repetitions

Then you choose the form of Fourier Series to be used: 
- Trigonometric form
- Compact Trigonometric form
- Exponential form

The script then prints the first 20 coefficients. 
And plots the function, spectrum and phase (considering the 20 coefficients).

Finally, for the reconstruction, you choose how many coefficients (may be
more than 20).

Then you get the final plot comparing the original and the reconstructed
signal.