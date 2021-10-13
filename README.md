[![](https://img.shields.io/badge/arXiv-2003.07355%20-red.svg)](https://arxiv.org/abs/2003.07355)

# CLASS_ROT: CLASS for ROTATED CMB POWER SPECTRA

A modified version of the publicly available Einstein-Boltzmann code [CLASS](https://github.com/lesgourg/class_public) to implement Rotated CMB Power Spectra.

## CLASS edited by
- Hongbo Cai
- Yilun Guan

## Files


## Installation

After cloning or downloading the repository, compile CLASS_ROT with make

$ make class

## Examples

### Python
Jupyter notebooks with worked out examples in Python can be found [here](https://github.com/mwt5345/class_ede/tree/master/class/notebooks-ede).

### C

CLASS_ROT can be run in C just as normal CLASS. See explanatory-ROT.ini for ROT implementation details.

$ ./class explanatory-ROT.ini

## Modifications to CLASS

Modifications are explained below.

(1) rotation.c

This module computes the rotated CMB anisotropy power spectra. Two parameters involve here
* alpha

The isotropic cosmic rotation angle {\bar{alpha}}

* A_cb
The amplitude of the cosmic rotation power spectrum.
