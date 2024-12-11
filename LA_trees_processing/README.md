## Calculating E-field and producing track angles for PixelAV simulation

```
git clone git@github.com:mroguljic/PixelCPETemplates.git
cd PixelCPETemplates/LA_tress_processing
./linkrootc++ calculateEfieldFromNtpl_template_2024b_prob
./calculateEfieldFromNtpl_template_2024b_prob #For the three prompts input 2, 2 and 0, respectively.
```
* You may need to adjust the boost libraries include path in the `linkrootc++`. Replace the `/opt/homebrew/Cellar/boost/1.85.0/include` path with the path to your local boost installation
* Only two Lorentz trees file are included in this example. The Run number of the Lorentz trees, the number of files to analyze, etc. can be adjusted in `analyze_LT2.txt` file. Only the first line matters
* Folder containing the templates used to reconstruct the hits should be specified in `calculateEfieldFromNtpl_template_2024b_prob.cxx` by changing the local variable: `std::string templates_dir = "templates_dir/"`
* IDs of the templates to be used for reconstruction are stored in local variables named bpl[1-4]m/p
* Output tracks are saved in files such as `c_R380947_cotangles_L1B.txt`. These files are contain track information (momenta) that are used in pixelav_simulation (next step)

## Generating and processing Lorentz trees from PixelAV simulation

* Copy simulated clusters to the `simulated_clusters` directory. File should be named such as `template_events_d28721.out`
* Edit the first line in `q_dist_2t.txt` to point to the input file with clusters. The first number if the CMS run number while the second one is the pixelav number. The script will try to find a cluster file `simulated_clusters/template_events_d$pixel_av_number.out`

```
./linkrootc++ make_ntuple_v1b
./make_ntuple_v1b #1109, -2, 0
```

Edit the file name at the beginning of `calculateLorentzAngleFromNtpl_calibrate_2017temp4.C` and run it with `root calculateLorentzAngleFromNtpl_calibrate_2017temp4.C`. This produces `LA_from_Ntpl_output/c_lorentzFit.txt` file. 

The contents of `c_lorenzFit.txt` files should be copied to a file for plotting, in this case `ef_comp2024_L1_IOV2a`. Below each copied content, insert a line such as `R380947_L1U green solid` or `dj0192g_f red dashed` to set the plotting style. To plot the comparisons of electric fields, run: `python tdr.py ef_comp2024_L1_IOV2a`
* Files with E fields based on simulated clusters are stored in `LA_from_Ntpl_output` while those calculated directly from data (Lorentz trees) are stored in `Efield_output`
* The first line in simulated_clusters .out file contains information about layer and E field model used in simulation.
* Model naming and tunes are kept in the `pixelav_lor_calibrations.txt` file

## Comparing simulated and measured clusters
* Edit the first line in `q_dist_2t.txt` to point to the input file with clusters. If this was done in the step of "generating and processing Lorentz trees from PixelAV simulation", it might already be correct
* Compile the two scripts for creating histograms:
    *  ./linkrootc++ test_q_dist_pt_bpix21_2t
    *  ./linkrootc++ test_standard_pt_bpix21_2t
*  `./test_standard_pt_bpix21_2t` asks for input. We can do: `0. 0. 0. 1109 1 2 0 0` for example
*  `./test_q_dist_pt_bpix21_2t` also asks for input. We can do: `0. 4 0. 1109 1 2 0 0` for example
* We look at the two pdfs and see if the agreement is acceptable. We proceed with simulating the clusters required for template production or work on tweaking the TCAD model. We focus on cluster sizes and charges.
* The agreement can be improved by tuning the parameters in `q_dist_2t.txt` 


## Testing the templates
After creating a new set of templates, we have a folder such as `PixelCPETemplates/templates_32600`. Besides the templates, it should have clusters simulated at uniform spread of angles (instructions in "Pixelav simulation after template production" part of pixelav_simulation module).
In `template_test_code.cxx`, edit the templates_dir variable to point to the folder with the templates and insert a new line at the beginning of `q_dist_2t.txt`. Then do:
 ```
./linkrootc++ template_test_code
template_test_code #Input it 1, 1116 1, 0
```
In the output file, we look at residuals for temp (template) and cmssw (generic). The offsets in x for the two should be close. If it is not, we can tweak "generror_summary" file. "54.4700" and "27.230000" are the charge width and mean Lorentz drift correction (usually set to half the charge width, but does not necessarily need to be). We want to keep the charge width = 2*drift correction while tweaking



