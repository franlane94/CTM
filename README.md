# Cosmological Trajectories Method (CTM) module

A Python module to compute the power spectrum and 2-point correlation function in real space using second-order CTM.

Please see the [Documentation](https://ctm-module.readthedocs.io) for more details

**Any queries please email**: <ctmcodehelp@gmail.com>

## Installation

### Required Modules

[mcfit](https://github.com/eelregit/mcfit) which can be installed using

```
pip install mcfit
```

[classylss](https://classylss.readthedocs.io/en/stable/) which can be installed using

```
conda install -c bccp classylss
```

or

```
pip install classylss
```
***

### Installation of the CTM module

The CTM module can be installed using

```
pip install ctm-module
```

if this method is used then `mcfit` and `classylss` will be automatically installed if they have not previously been installed. To check the code has been installed properly run in the terminal

```
python
from ctm import CTM
```

***

## Parameters
### Cosmological parameters

The following parameters can be specified by the user using

```
CTM(omega0_cdm=0.25, omega0_b=0.05, h=0.7, n_s=0.96, sigma_8=0.8)
```

| Parameter     |  Automatic Value |
| ------------- |:-------------:|
| <img src="https://latex.codecogs.com/gif.latex?\Omega_{cdm}h^2" />  | 0.11933|
| <img src="https://latex.codecogs.com/gif.latex?\Omega_bh^2" />  | 0.02233     |
| <img src="https://latex.codecogs.com/gif.latex?h" />      | 0.6737      |
| <img src="https://latex.codecogs.com/gif.latex?n_s" />      | 0.9665    |
| <img src="https://latex.codecogs.com/gif.latex?\sigma_8" />      | 0.8102  |

### CTM parameters

The following CTM parameters can be specified by the user using

```
CTM().ctm_power(epsilon=0.01, z_init=150.0, k_c=5.0)
```

| Parameter     |  Automatic Value | Description |
| ------------- |:-------------:|:--------------|
| <img src="https://latex.codecogs.com/gif.latex?n_k" />  | 3000| Number of k values used in calculation of the power spectrum |
| <img src="https://latex.codecogs.com/gif.latex?\epsilon_\mathrm{BZ}" />  | 1.0   | The value of the expansion parameter|
| <img src="https://latex.codecogs.com/gif.latex?z_\mathrm{init}" />      | 100.0      | The value of the initial redshift from which the time dependence is integrated from |
| <img src="https://latex.codecogs.com/gif.latex?k_c" />      | 0.0   | The value of the cutoff k value for using an initial Gaussian damped power spectrum |
| <img src="https://latex.codecogs.com/gif.latex?n_\mathrm{val}" />      | 32   | The number of spherical Bessel functions used in the calculation |

The range of values for the CTM parameters are:

<img src="https://latex.codecogs.com/gif.latex?0\leq\epsilon_\mathrm{BZ}\leq1" />
<img src="https://latex.codecogs.com/gif.latex?0\leq{k_c}\leq50" />

## Documentation

### Example Jupyter Notebook

An example notebook can be found in the folder Examples

**Real-space examples**: [CTM_real_space_examples.ipynb](https://github.com/franlane94/CTM/blob/master/Examples/CTM_real_space_examples.ipynb)
