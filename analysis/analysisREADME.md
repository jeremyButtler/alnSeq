# Use

This folder contains the scripts and genomes needed to
  replicated my quick benchmarking of alnSeq.

Note: My benchmarking was done a Linux Debian machine, so
  I do not know if all of the programs will install on
  OpenBSD.

# Installing programs for benchmarking

You will have to install emboss separately, but this
  is easily installed by a local package manager.
  `apt-get install emboss` for Debian. The remaining
  packages can be downloaded and compiled by using the
  make file in this directory. The Makefile downloads
  bio-alignment, Complete-Striped-Smith-Waterman-Library,
  and seq-align. You can remove these folders and the 
  compiled binaries with `make clean`.

# Running the analysis

Once the programs are installed you can run this analysis
  with the benchAln.sh bash script and the
  alnSeq-benchmarking.r R script. You will needed R, with
  the ggplot2, ggpubr, data.table, viridis, and svglite
  libraries to build the graphs.

```
# Commands to run this analysis

# install emboss
make
make run
```
