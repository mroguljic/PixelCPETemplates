# Scripts to produce CPE templates for the CMS pixel detector
1. Process input Lorentz angle (LA) trees in the `LA_trees_processing` module and obtain track momenta and electric field (and some other useful data)
2. Using track lists from 1. and a TCAD E-field model, simulate clusters with PixelAv in the `pixelav_simulation` module
3. Generate Lorentz trees from simulated clusters, again in the `LA_trees_processing` and plot the simulated E fields to compare them with those from data
4. Generate new TCAD models and repeat steps 2 and 3 until satisfied with agreement
   - Processing of TCAD simulation is described in `tcad_file_processing`
   - We can run a "test_diode" run which is quick nad offers general insight into the E field. The field can be plotted directly using the output of TCAD, without doing pixelav simulation
   - Once we have a good agreement with the quick simulation, we can run TCAD + pixelav for full simulation
5. Simulate clusters with a wide range of angles for producing templates
   - We run the "Pixelav simulation for template production" part in  the `pixelav_simulation` module
   - Have to make sure the pixel2.init (e field) is the one we obtained for the desired model in the `tcad_file_processing` module
6. Produce templates
   - Copy the simulated clusters (from the `finished` directory) to your local machine and run `gunzip` on them
   - Follow https://github.com/OzAmram/PixelTemplateProduction 
