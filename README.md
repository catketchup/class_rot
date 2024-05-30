[![](https://img.shields.io/badge/arXiv-2111.14199%20-red.svg)](https://arxiv.org/abs/2111.14199)

# CLASS_ROT: CLASS for ROTATED CMB POWER SPECTRA

A modified version of the publicly available Einstein-Boltzmann code [CLASS](https://github.com/lesgourg/class_public) to implement Rotated CMB Power Spectra.

See [Cai et al.](https://arxiv.org/abs/2111.14199) where CLASS_ROT is introduced and tested.

![](https://github.com/catketchup/class_rot/blob/main/figures_ROT/ps_sims.png) <!-- .element height="6%" width="13.5%" -->

![](https://github.com/catketchup/class_rot/blob/main/figures_ROT/cl_lens_rot.png) <!-- .element height="6%" width="13.5%" -->

## CLASS edited by
- Hongbo Cai; hoc34@pitt.edu
- Yilun Guan; yilun.guan@dunlap.utoronto.ca

## Files


## Installation

After cloning or downloading the repository, compile CLASS_ROT with make

$ make

## Important Parameters
(1) alpha:

The isotropic cosmic rotation angle

(2) A_cb:

The amplitude of the scale-invariant Gaussian random cosmic rotation power spectrum

(3) input_claa:

The input generic rotation power spectrum

'A_cb' and 'input_claa' should not be given together.

## Examples

### Python
Jupyter notebooks with worked out examples in Python can be found [here](https://github.com/catketchup/class_rot/tree/main/notebooks_rot).

[cl_rot.ipynb](https://github.com/catketchup/class_rot/tree/main/notebooks_rot/cl_rot.ipynb) is for using a scale-invariant rotation power spectrum.

[use_input_claa.ipynb](https://github.com/catketchup/class_rot/tree/main/notebooks_rot/use_input_claa.ipynb) is for using an input generic rotation power spectrum.

[cl_rot_lens.ipynb](https://github.com/catketchup/class_rot/blob/main/notebooks_rot/cl_rot_lens.ipynb) is for generating rotated-lensed Cl,  and comparing it with lensed-raw Cl and rotated-raw Cl.

### C

CLASS_ROT can be run in C just as normal CLASS. See explanatory-ROT.ini for ROT implementation details.

$ ./class explanatory_ROT.ini

## New package added to CLASS
Modifications are explained below.

(1) [rotation.c](https://github.com/catketchup/class_rot/blob/main/source/rotation.c):

This module computes the rotated CMB anisotropy power spectra.


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
