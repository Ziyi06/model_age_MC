# model_age_MC
Estimates on zircon Hf model ages and their uncertainties using Monte Carlo bootstrapping, used in Zhu et al. in submission.

This module provides a model age and a εHf(t) along with their errors, by bootstrapping (i) uncertainty in the Hf isotopic composition of the mantle reservoir, (ii) uncertainty in the 176Lu/177Hf for the crustal source region, and (iii) uncertainty in measurements. One can decide which uncertainty(s) to be considered. The minimum requirement of data is a zircon U-Pb age, 176Hf/177Hf, and 176Lu/177Hf.

T_boot() is the core function. The detailed parameter description can be found in the comments. Run "python bootstrap_plot.py" as a demo.

If you already have an array containing a distribution 176Hf/177Hf of arc / depleted mantle, load it and set it to the parameter "mantle_hh" in T_boot(). If not, you can also leave it defaulted as Iizuka's result (AM). Note that the parameter plot can only be True when the inputs are all float. If the inputs are array-like, remember to switch it to False.
