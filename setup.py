from setuptools import setup

setup(name='ctm',
      version='1.1',
      description='Module to compute the power spectrum and two-point correlation function using the CTM',
      url='https://ctm.readthedocs.io',
      author='Fran Lane',
      author_email='ctmcode@gmail.com',
      license='MIT',
      packages=['ctm', 'ctm.cosmology', 'ctm.misc', 'ctm.power_spectrum', 'ctm.redshift_space'],
      install_requires=['mcfit', 'classylss', 'scipy'])
