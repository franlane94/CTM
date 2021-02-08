from setuptools import setup

setup(name='ctmmodule',
      version='1.12',
      description='Module to compute the power spectrum and two-point correlation function using the CTM',
      url='https://ctm-module.readthedocs.io',
      author='Fran Lane',
      author_email='ctmcodehelp@gmail.com',
      license='MIT',
      packages=['ctm', 'ctm.cosmology', 'ctm.misc', 'ctm.power_spectrum'],
      install_requires=['mcfit', 'classylss', 'scipy'])
