# MidrangeFilter
Implementation of the Midrange filter and code used in the publication "Aggregating SNPs Improves Filtering for False Positive Associations Post-Imputation", currently under Review

The folder MidrangeFilter contains an implementation of the new method and an example data set with instructions how to run it.

The folder code_publication_midrange_filter contains functions and files used to produce the simulations presented in the publication. To run this code, file_paths.R and configuration files need to be adjusted to the specific system the simulation is supposed to run on. The file example_simulation.R contains the main function calls to run the simulation for one example simulation setting. 

Hapgen2: https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html

Beagle5.2: https://faculty.washington.edu/browning/beagle/b5_2.html

ImputAccur: https://gitlab.gwdg.de/kolja.thormann1/imputationquality.git

Bcftools: https://samtools.github.io/bcftools/
