## Calculating E-field and producing track angles for PixelAV simulation

```
git clone git@github.com:mroguljic/PixelCPETemplates.git
cd PixelCPETemplates/LA_tress_processing
./linkrootc++ calculateEfieldFromNtpl_template_2024b_prob
./calculateEfieldFromNtpl_template_2024b_prob #For the three prompts input 2, 2 and 0, respectively.
```
* You may need to adjust the boost libraries include path in the `linkrootc++`. Replace the `/opt/homebrew/Cellar/boost/1.85.0/include` path with the path to your local boost installation.
* Only two Lorentz trees file are included in this example. The Run number of the Lorentz trees, the number of files to analyze, etc. can be adjusted in `analyze_LT2.txt` file. Only the first line matters
* File containing the templates used to reconstruct the hits  should be specified in `calculateEfieldFromNtpl_template_2024b_prob.cxx` by changing the local variable: `std::string templates_dir = "templates_dir/"`

