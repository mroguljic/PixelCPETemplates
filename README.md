# Scripts to produce CPE templates for the CMS pixel detector
1. Process input Lorentz angle (LA) trees in the `LA_trees_processing` module and obtain track momenta and electric field (and some other useful data)
2. Using track lists from 1. and a TCAD E-field model, simulate clusters with PixelAv in the `pixelav_simulation` module
3. Generate Lorentz trees from simulated clusters, again in the `LA_trees_processing` and plot the simulated E fields to compare them with those from data
4. Generate new TCAD models and repeat steps 2 and 3 until satisfied with agreement
   - Processing of TCAD simulation is described in `tcad_file_processing`
   - We can run a "test_diode" run which is quick nad offers general insight into the E field. The field can be plotted directly using the output of TCAD, without doing pixelav simulation
   - Once we have a good agreement with the quick simulation, we can run TCAD + pixelav for full simulation.
5. Compare the simulated clusters based on the selected model with the measured clusters.
   - Described in "Comparing simulated and measured clusters" in the `LA_trees_processing` module
7. Simulate clusters with a wide range of angles for producing templates
   - We run the "Pixelav simulation for template production" part in  the `pixelav_simulation` module
   - Have to make sure the pixel2.init (e field) is the one we obtained for the desired model in the `tcad_file_processing` module
8. Produce templates using the `template_production` module
   - Copy the simulated clusters (from the `finished` directory) to your local machine and run `gunzip` on them
   - Create `pix_2t.proc` file in local directory where the clusters are stored. Example of the file is stored in `template_production`
   - From the same directory, run gen_xy_template followed by gen_zp_template (the compiled binaries should be in `template_production/bin`
9. Test the new templates
   - Simulate tracks at uniform spread of angles (instructions in "Pixelav simulation after template production" part of `pixelav_simulation` module)
   - Follow "Testing the templates in `LA_trees_processing` module
