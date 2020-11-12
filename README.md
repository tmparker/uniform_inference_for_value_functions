## Replication files for "Uniform inference for value functions"

First commit 2020-07-03.

This repository contains files to replicate all the figures in the paper.  This
includes both the simulation experiment results and the empirical illustration.

### The experiment directories confband and sd

These directories contain the code that was used to run simulation experiments.
They depend on the file `bound.cpp`, which is one level up.  The files that did
the simulations are `sim_cb.R` and `sim_sd.R`.  These files were run with 100 
random seeds and 10 repetitions each using the `rlecluyer`  R package (the SLURM
shell script used to run the simulations is also included in each directory, as
is a file called `glue.R` which was used to paste all the results together at
the end).  The resulting pvalue data is packaged in `.rda` files labeled by 
sample size.  The files `cband_powercurves.R` and `sd_powercurves.R` produce pdf 
files of the same name.

### The experiment directory par_choose

The directory `par_choose` is a little different &mdash; it is like the two
experiments but it contains its own `bound_localpower.cpp` file that was used in
choosing the parameter "K" that appeared in the formula for $a_n$ in the main
text.

### The empirical directory

This directory contains the code used to produce all the example figures.  The
file `lalonde.RData` contains the original observations.  The file
`ex_lalonde.R` has all the R commands used to produce the rest of the figure
files that are in the directory (this file depends on `bound.cpp` in the
higher-level directory).  It also produces the side product 
`lalonde_cdf_data.rda`, which contains the data produced by `ex_lalonde.R` up to
that point (before the figures get plotted).

### Common files

Everything depends on the file `bound.cpp`, which contains C++ code implementing
the methods for the examples.
