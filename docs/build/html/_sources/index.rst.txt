.. GCTM documentation master file, created by
   sphinx-quickstart on Thu Oct 29 11:29:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the GCTM module's documentation!
================================

The General Cosmological Trajectories Method (GCTM) is a perturbative technique to compute the power spectrum and two-point correlation function in both real and redshift space in the mildly to non-linear regime for a range of cosmologies.

Source and installation
-----------------------

The code is available at `Github <https://github.com/franlane94/GCTM>`__ but can also be installed using pip (all of the dependencies will be automatically installed)

.. code-block:: bash

   pip install gctm


To install the code from source, first run the following

.. code-block:: bash

  pip install mcfit
  pip install classylss
  git clone https://github.com/franlane94/GCTM

Then navigate to your .bash_profile (or .profile on Linux machines) and add the following line

.. code-block:: bash

  export PYTHONPATH="$PYTHONPATH:/DIRECTORY_WHERE_GCTM_IS_INSTALLED"


To check the module has been installed correctly run

.. code-block:: bash

  python
  from gctm import GCTM

Getting started
---------------

For a quickstart see the `Jupyter notebooks <https://github.com/franlane94/GCTM/tree/master/Examples>`__. For more detailed examples and information about the code please see the documentation pages.

.. note::

  The units of commonly used variables in this code are:

    * :math:`H_0\ [\mathrm{km}\ \mathrm{s}^{-1}\ \mathrm{Mpc}^{-1}]`
    * :math:`k\ [\mathrm{h}\ \mathrm{Mpc}^{-1}]`
    * :math:`\mathrm{P}\left(k\right)\ [\mathrm{Mpc}^3\ \mathrm{h}^{-3}]`

The default cosmology is `Planck18 <https://arxiv.org/abs/1807.06209>`__ and the values of key parameters shown in the Table below.

+-------------------------+--------+
| Parameter               | Value  |
+=========================+========+
| :math:`H_0`             | 67.66  |
+-------------------------+--------+
| :math:`\Omega_{cdm}h^2` | 0.1193 |
+-------------------------+--------+
| :math:`\Omega_{b}h^2`   | 0.0224 |
+-------------------------+--------+
| :math:`n_s`             | 0.9665 |
+-------------------------+--------+
| :math:`\sigma_8`        | 0.8102 |
+-------------------------+--------+


Citation
--------

The GCTM code is described in Lane et al. (2020). Please cite this paper if you use this module.

Help and issues
---------------

If you encounter any issues please either use the `Github issues <https://github.com/franlane94/GCTM/issues>`__  page, or email gctmcode@gmail.com.

.. toctree::
   :maxdepth: 1
   :hidden:

   basic_cosmology.rst
   gctm.rst
   gctm_redshift_space.rst
