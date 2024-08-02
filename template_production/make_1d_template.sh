#!/bin/sh
# A simple script to do both steps of making 1D templates
# Should be run in a directory with pixelav events (named template_events_dXXXXX.out) 
# Requires a config file named pix_2t.proc in the same directory
# Takes one argument, which is the location of the directory with executables in it 

exec_dir=$1

./${exec_dir}gen_xy_template
./${exec_dir}gen_zp_template
