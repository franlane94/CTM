.. CTM documentation master file, created by
   sphinx-quickstart on Thu Oct 29 11:29:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the CTM module's documentation!
================================

The Cosmological Trajectories Method (GCTM) is a perturbative technique to compute the power spectrum and two-point correlation function in real space in the mildly to non-linear regime for a range of cosmologies.

Source and installation
-----------------------

The code is available at `Github <https://github.com/franlane94/CTM>`__ but can also be installed using pip (all of the dependencies will be automatically installed)

.. code-block:: bash

   pip install ctm


To install the code from source, first run the following

.. code-block:: bash

  pip install mcfit
  pip install classylss
  git clone https://github.com/franlane94/CTM

Then navigate to your .bash_profile (or .profile on Linux machines) and add the following line

.. code-block:: bash

  export PYTHONPATH="$PYTHONPATH:/DIRECTORY_WHERE_CTM_IS_INSTALLED"


To check the module has been installed correctly run

.. code-block:: bash

  python
  from ctm import CTM

Getting started
---------------

For a quickstart see the `Jupyter notebook <https://github.com/franlane94/CTM/tree/master/Examples>`__. For more detailed examples and information about the code please see the documentation pages.

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

The CTM module is described in Lane et al. 2021, if you make use of the CTM module, please include a citation to this work in any publications.

Acknowledgements
----------------

Please also cite the following papers whose work has featured in the development of both the CTM and this Python module:

*: `Github <https://github.com/franlane94/CTM>`__


Help and issues
---------------

If you encounter any issues please either use the `Github issues <https://github.com/franlane94/CTM/issues>`__  page, or email ctmcodehelp@gmail.com.

.. toctree::
   :maxdepth: 1
   :hidden:

   basic_cosmology.rst
   ctm.rst
