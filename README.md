# Scripts to produce CPE templates for the CMS pixel detector
1. Process input Lorentz angle (LA) trees in the `LA_trees_processing` module and obtain track momenta
2. Using track lists from 1. and a TCAD E-field model, simulate clusters with PixelAv in the `pixelav_simulation` module
3. Generate Lorentz trees from simulated clusters, again in the `LA_trees_processing` and plot the simulated E fields to compared them with those from data
4. Generate new TCAD models and repeat steps 2 and 3 until satisfied with agreement
