## Calculating E-field and producing track angles for PixelAV simulation

```
git clone git@github.com:mroguljic/PixelCPETemplates.git
cd PixelCPETemplates/LA_tress_processing
./linkrootc++ calculateEfieldFromNtpl_template_2024b_prob
./calculateEfieldFromNtpl_template_2024b_prob #For the three prompts input 2, 2 and 0, respectively.
```
* You may need to adjust the boost libraries include path in the `linkrootc++`. Replace the `/opt/homebrew/Cellar/boost/1.85.0/include` path with the path to your local boost installation
* Only two Lorentz trees file are included in this example. The Run number of the Lorentz trees, the number of files to analyze, etc. can be adjusted in `analyze_LT2.txt` file. Only the first line matters
* File containing the templates used to reconstruct the hits  should be specified in `calculateEfieldFromNtpl_template_2024b_prob.cxx` by changing the local variable: `std::string templates_dir = "templates_dir/"`
* Output tracks are saved in files such as `c_R380947_cotangles_L1B.txt`. These files are needed for further steps of template production

## Generating and processing Lorentz trees from PixelAV simulation

* Copy simulated clusters to the `simulated_clusters` directory. File should be named such as `template_events_d28721.out`
* Edit the first line in `q_dist_2t.txt` to point to the input file with clusters. The first number if the CMS run number while the second one is the pixelav number. The script will try to find a cluster file `simulated_clusters/template_events_d$pixel_av_number.out`

```
./linkrootc++ make_ntuple_v1b
./make_ntuple_v1b #1109, -2, 0
```

Edit the file name at the beginning of `calculateLorentzAngleFromNtpl_calibrate_2017temp4.C` and run it with `root calculateLorentzAngleFromNtpl_calibrate_2017temp4.C`
