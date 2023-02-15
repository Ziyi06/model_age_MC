# model_age_MC
This module implements the Monte Carlo bootstrapping method to calculate zircon Hf model ages and their uncertainties, as proposed by Zhu et al. (2023) [https://doi.org/10.1016/j.gca.2023.02.005].

This module calculates zircon Hf model ages and their uncertainties by taking into account (i) uncertainty in the Hf isotopic composition of the mantle reservoir, (ii) uncertainty in the 176Lu/177Hf for the crustal source region, and (iii) uncertainty in measurements. The user can choose which uncertainty(s) to include. The minimum data required is a zircon U-Pb age, a 176Hf/177Hf, and a 176Lu/177Hf.

Try running "python bootstrap_plot.py" for a demo. The main function is T_boot(). The parameter descriptions can be found in the code comments.

If you have an array containing a 176Hf/177Hf distribution of arc/depleted mantle, load it and set it to the parameter "mantle_hh" in T_boot(). If not, you can use the recommended data compilation of modern arc mantle or leave it as the default constant value. 

Note: The "plot" parameter can only be set to True when the inputs are floats. If the inputs are arrays, make sure to set "plot" to False.
