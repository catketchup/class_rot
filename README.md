[![](https://img.shields.io/badge/arXiv-2111.14199%20-red.svg)](https://arxiv.org/abs/2111.14199)

# CLASS_ROT: CLASS for ROTATED CMB POWER SPECTRA

A modified version of the publicly available Einstein-Boltzmann code [CLASS](https://github.com/lesgourg/class_public) to implement Rotated CMB Power Spectra.

See [Cai et al.](https://arxiv.org/abs/2111.14199) where CLASS_ROT is introduced and tested.

![](https://github.com/catketchup/class_rot/blob/main/figures_ROT/ps_sims.png) <!-- .element height="6%" width="13.5%" -->

## CLASS edited by
- Hongbo Cai; hoc34@pitt.edu
- Yilun Guan; yilun.guan@dunlap.utoronto.ca

## Files


## Installation

After cloning or downloading the repository, compile CLASS_ROT with make

$ make class

## Examples

### Python
Jupyter notebooks with worked out examples in Python can be found [here](https://github.com/catketchup/class_rot/tree/main/notebooks_rot).

### C

CLASS_ROT can be run in C just as normal CLASS. See explanatory-ROT.ini for ROT implementation details.

$ ./class explanatory_ROT.ini

## Modifications to CLASS
Modifications are explained below.

(1) rotation.c:

This module computes the rotated CMB anisotropy power spectra. Two parameters involve here

(2) alpha:

The isotropic cosmic rotation angle

(3) A_cb:
The amplitude of the cosmic rotation power spectrum.

## Cite us
If you use our code in a published work, please cite our paper:

```
@ARTICLE{classrot,
       author = {{Cai}, Hongbo and {Guan}, Yilun},
        title = "{Computing Microwave Background Polarization Power Spectra from Cosmic Birefringence}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics},
         year = 2021,
        month = nov,
          eid = {arXiv:2111.14199},
        pages = {arXiv:2111.14199},
archivePrefix = {arXiv},
       eprint = {2111.14199},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021arXiv211114199C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
