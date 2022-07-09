# model_age_MC
This package gives estimates on model ages and their uncertainties using Monte Carlo bootstrapping. [https://github.com/zidianjun/BOA_calculation](url)

Manual:

  Click the green button "clone or download" to get the whole repository into a local device.
  
Make sure that Python 2 has been already installed, as well as numpy and matplotlib.

Copy and Paste your primary melt composition and cpx-melt P-T data into ./data/primary_melt_composition.txt and ./data/cpx_melt.txt, respectively.

Change the parameters in config.py. The changeable parameters include a melt Fe3+ proportion (Fe3+/FeT), errors of the thermobarometers, intended H2O content range and increment, and the size of unit.

Run "python demo.py" to get the outcome.

The outcome is a figure showing a scatter of H2O content (.wt%) and BOA as well as printing them on the screen.
