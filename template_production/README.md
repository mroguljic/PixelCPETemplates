# Pixel Template Production
Code for production of pixel templates for CMS, adapted from: https://github.com/OzAmram/PixelTemplateProduction

## Compiling

Everything is run independent from the cmssw environment but makes use of some of the code in it. Current versions of cmssw no longer support the standalone pixel template code so we store the compatible version in `cmssw_code`


The script **fetch\_cmssw\_code.sh** grabs the latest version of all the needed pixel code from the [cmssw github](https://github.com/cms-sw/cmssw).
If you want to grab from a branch other than the cmssw master (eg to test some
changes), you can change the `branch` variable in the script to point to a different branch. **Current cmssw version cannot be run standalone so if you do this, the code will not compile**

The code shared with CMSSW uses the boost libraries. Change the include path in `src/Makefile` to match your `boostlib` include folder

## Making Templates

Enter the `PixelCPETemplates/templates_XXXXX` folder. It should have been created when transferring the simulated clusters from remote to local machine ("Pixelav simulation for template production" of the pixelav_simulation module). This folder should have all 205 `template_events_*out` files. One also needs to manually create `pix_2t.proc` file. Then, from that folder, execute:
```
../template_production/bin/gen_xy_template
../template_production/bin/gen_zp_template
```

**Below might be outdated!**

There are two simple bash scripts that run the necessary executables to make templates: **make\_1d\_templates.sh** and **make\_2d\_templates.sh**. 

They should be run inside a directory containing pixelav events and a config file named `pix_2t.proc`. They take 1 argument, which is the location of the bin/ directory. 

Many files are produced when making templates but for the 1d production the real output is one template file named `template_summary_zpXXXX.out`and one gen errors file named `generror_summary_zpXXXX.out`. For 2d production it is one 2d template file named `template_summary2D_zpXXXX.out`. (The XXXX will be the starting file number in your config). 

The format `pix_2t.proc` is as follows:

> start\_file nfiles noise thresh1 thresh2 thresh1\_noise\_frac common\_noise\_frac gain\_noise\_frac readout\_noise frontend\_type

> use\_l1\_offset write\_header xtalk\_frac xtalk\_spread do\_cluster\_healing


> id NTy NTyx NTxx DType Bfield VBias temp fluenc qscale 

Note that NTy is not used by the 2D templates so its value doesn't matter, but to keep the format consistent something must be there. 
Using 0 for xtalk\_frac will turn off cross talk. 
Extra parameters on any of the lines will be ignored. 

An example config for 1D barrel templates is: 

> 58401 205 250. 1600. 1600. 0.073 0.080 0.080 350. 0

> 0 1 0.0 0.0 0

> 900 60 5 29 0 3.8 125. 263. 0. 1.


Additional lines can be added which turn on the creation of the
PixelResolutionHistograms and define their binning. 

For **gen_zp_template** (which stores resolutions in 2x 1D bins of fixed width) this additional line should be structured as:

> name CotBetaBinWidth CotBetaStart nCotBetaBins CotAlphaBinWidth CotAlphaStart nCotAlphaBins


For **gen_zp_template2d** (which stores resolutions in 2D bins) these additional lines should be structured as:

> name
> CotBetaBinEdge1  CotBetaBinEdge2  CotBetaBinEdge3  CotBetaBinEdge4  ...
> CotAlphaBinEdge1  CotAlphaBinEdge2  CotAlphaBinEdge3  CotAlphaBinEdge4  ...

Where the 2nd and 3rd rows can be as long as desired and list the bin edges for
the CotBeta and CotAlpha binning.
Examples can be found in the config\_db directory.


## Description of Executables
**gen\_xy\_template** : Makes 1D projections of average charge distributions from pixelav events. Also does charge variance fits. Makes files like `ptemp_XXXX.txt` and `ztemp_XXXX.txt`.

**gen\_xy\_template2d** : Makes 2d projections of average charge distributions from pixelav events.  Makes files like `zptemp_XXXX.txt`.

**gen\_zp\_template**: Uses the 1D projections and pixelav events. Runs the generic and template reco (using CMSSW code) on pixelav events to get resolution different algorithms and compute corrections to be saved in templates. Outputs one template file named `template_summary_zpXXXX.out`and one gen errors file named `generror_summary_zpXXXX.out`.

**gen\_zp\_template2d**: Uses the 1D and 2D projections and pixelav events. Runs the 2D template reco (using CMSSW code) on pixelav events to get resolution. Outputs one 2d template file named `template_summary2D_zpXXXX.out`.

**compare\_templates**: Takes in the file names of two templates and checks that all numerical values are the same within some threshold (default is 10^-5). It lists any discrepancies with the line number for investigation. Useful for testing changes. 

**test_template**: Uses pre-made 1D templates to run local version of CMSSW 1D template reco and makes various plots. Useful for testing a new set of 1D templates. 
Should be run a directory with template\_events, generror and template\_summary files. Also takes a config called `test_params.txt`.
The first line of the config is the same as the `pix_2t.proc` but without the
nfiles parameter (because it will only use one file). 
The second line has six  parameters:

> nFile use_l1_offset xtalk_frac xtalk_noise do_cluster_healing do_IBCs do_2d_templates

nFile is the file number of the template (the XXXXX) and the
second is the `use_l1_offset` parameter, xtalk_frac is the fractional charge sharing for the cross talk, the xtalk_noise is the gaussian spread on the central value of the charge sharing, do_cluster_healing turns on the cluster healing before the CPEs, do_IBCs turns on the irradiation bias corrections for the generic algorithm, do_2d_templates tests the 2d reco as well and produces some additional plots (2d templates must be provided).



An example `test_params.txt` config is:
> 58606 150. 1600. 1600. 0.073 0.080 0.08 350. 0

> 58401 0 0.0 0.0 1 0 0


