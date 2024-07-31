## Pixelav simulation of clusters based on LA trees

The goal of this module is to take tracks, obtained from Lorentz angle trees, and simulate corresponding clusters using PixelAv. Besides the tracks, a required input is the electric field of the sensor that is obtained using TCAD simulations.

Jobs are submitted from the `batch` folder using .pbs files. For example:
```
cd batch
sbatch slurm_submit.pbs
squeue -u $USER
```

The files that need to be provided are listed in the .pbs file. We need a file with the electric field (needs to be called `pixel2.init`) and a file with the weighting potential 8 (`wgt_pot.init`), both obtained from TCAD. The third file is a list of tracks that is obtained by processing LA trees.

The number of jobs sent is specified in the `#SBATCH --array=1-14` lines. This line would submit 14 jobs. We typically process 30 thousand tracks per job so one can count the number of tracks in the `track_list.txt` file to gauge the number of required jobs. The last number in the second line of the `pixel2.init` file sets the initial (PixelAv) run number for the jobs. This is used for bookeeping as the output of the jobs will have `template_events_d$RUN$JOB_ID.out.gz` format. Typically, we separate jobs for different simulations by 100 runs.

The output is stored in `finished` directory. Files from a single run should be concatenated
```
gunzip template*
# For easier bookkeeping, add layer + CMS run number to the end of the first line of first file. For example: L1U 380947 
cat template_events_d285??.out > template_events_d28721.out 
```

## Pixelav simulation for template productio
Here we simulate tracks at predefined angles (expressed as a grid of track sizes when projected to two local axes of the sensor). The procedure is the same as above except that we use an additional file, `runlist.init`, containing the grid of angles. The corresponding submission file is `smake_template_2f_bpix2.pbs`
