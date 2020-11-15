# General Cosmological Trajectories Method (GCTM) code

A Python code to compute the power spectrum and 2-point correlation function in real and redshift space using the GCTM.

Please see the [Wiki](https://github.com/franlane94/GCTM/wiki) for more details

**Any queries please email**: <flane@roe.ac.uk>

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

### Installation of GCTM code

Once these dependencies have been installed use

```
git clone https://github.com/franlane94/GCTM
```

to install the GCTM code. Then navigate to your .bash_profile (or .profile on Linux machines) script and add the following line

```
export PYTHONPATH="$PYTHONPATH:/DIRECTORY_WHERE_GCTM_IS_INSTALLED"
```
To check the code has been installed properly run

```
python
from GCTM.gctm import GCTM
```
***

## Parameters
### Cosmological parameters

The following parameters can be specified by the user using


| Parameter     |  Automatic Value |
| ------------- |:-------------:|
| <img src="https://latex.codecogs.com/gif.latex?\Omega_{cdm}h^2" />  | 0.11977|
| <img src="https://latex.codecogs.com/gif.latex?\Omega_bh^2" />  | 0.02233     |
| <img src="https://latex.codecogs.com/gif.latex?h" />      | 0.6737      |
| <img src="https://latex.codecogs.com/gif.latex?n_s" />      | 0.9652    |
| <img src="https://latex.codecogs.com/gif.latex?\sigma_8" />      | 0.8101  |

```
GCTM(omega0_cdm=0.25, omega0_b=0.05, h=0.7, n_s=0.96, sigma_8=0.8)
```
### GCTM parameters

The following GCTM parameters can be specified by the user using

| Parameter     |  Automatic Value | Description |
| ------------- |:-------------:|:--------------|
| <img src="https://latex.codecogs.com/gif.latex?n_k" />  | 3000| Number of k values used in calculation of the power spectrum |
| <img src="https://latex.codecogs.com/gif.latex?\epsilon" />  | 1  .0   | The value of the expansion parameter|
| <img src="https://latex.codecogs.com/gif.latex?z_{init}" />      | 99.0      | The value of the initial redshift from which the time dependence is integrated from |
| <img src="https://latex.codecogs.com/gif.latex?k_c" />      | 0.0   | The value of the cutoff k value for using an initial Gaussian damped power spectrum |
| <img src="https://latex.codecogs.com/gif.latex?\mu_k" />      | 0.0 | The value of the line-of-sight parameter for calculating the redshift-space power spectrum |

```
GCTM().gctm_power(epsilon=0.01, z_init=150.0, k_c=5.0)
```

The range of values for the GCTM parameters are:

<img src="https://latex.codecogs.com/gif.latex?0\leq\epsilon\leq1" />
<img src="https://latex.codecogs.com/gif.latex?0\leq{k_c}\leq50" />
<img src="https://latex.codecogs.com/gif.latex?-1\leq\mu_k\leq1" />

## Documentation

### Example Jupyter Notebooks

Example notebooks can be found in the folder Examples

**Real-space examples**: [GCTM_real_space_examples.ipynb](https://github.com/franlane94/GCTM/blob/master/Examples/GCTM_real_space_examples.ipynb)

**Redshift-space examples**: [GCTM_redshift_space_examples.ipynb](https://github.com/franlane94/GCTM/blob/master/Examples/GCTM_redshift_space_examples.ipynb)

### Example Python scripts

Example Python scripts can be found in the folder Examples

**Real-space examples**: [GCTM_real_space_examples.py](https://github.com/franlane94/GCTM/blob/master/Examples/GCTM_real_space_examples.py)

**Redshift-space examples**: [GCTM_redshift_space_examples.py](https://github.com/franlane94/GCTM/blob/master/Examples/GCTM_redshift_space_examples.py)
