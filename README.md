Reversible Jump PDMP Samplers
=============================

# Data

<!--
Provide a short (< 100 words), high-level description of the data
-->

Examples do not use external data. R scripts to simulate data and simulated .Rdata files are available in the `data` directory. See the Readme in the `data` directory for details.

# Code

R code implementing the experiments is located in the `code` directory -- see the readme located there for details. Subdirectories contain code for reproducing the tables and figures, the simulation study and generating the data. All simulations were run on the Lyra high computing cluster at QUT using 3780 x 64bit Intel Xeon Cores with the SUSE Linux operating system (<https://cms.qut.edu.au/__data/assets/pdf_file/0012/388785/high-performance-computing.pdf>). The package rjpdmp should be installed as described in the `code/README.md`.

<!--
Provide a short (< 100 words), high-level description of the code. If necessary, more details can be provided in files that accompany the code. If no code is provided, please state this and say why (e.g., if the paper contains no computational work).
-->
#### Version of primary software used

Due to version conflict with the R package rstan on the computing cluster we used two versions of R for the experiments. R version 3.4.1 was used for the logistic regression examples and R version 3.6.2 was used for the robust regression examples (where stan was a competitor).
Details on the R version number of hours, memory and number of parallel jobs are given in the `code/Job Submission Scripts` directory.

#### Libraries and dependencies used by the code

The R package `rjpdmp` created to implement the methods is given in the `code` directory. Packages rstan_2.19.2 and nimble_0.10.1 were used in simulation studies. Additional packages and dependencies include:

1. Rcpp_1.0.6 (for R 3.4.1) and Rcpp_1.0.6 (for R 3.6.2)
2. data.table_1.12.8
3. ggplot2_3.3.5
4. gridExtra_2.3   
5. Rmisc_1.5.1     
6. rstan_2.19.2
7. nimble_0.10.1
8. devtools_2.4.2
9. RcppArmadillo_0.10.2.1.0

<!--
Include version numbers (e.g., version numbers for any R or Python packages used)
-->

## Workflow


<!--
Indicate where the materials (generally including the code, unless in a separate location and indicated in the previous section) are available. We strongly encourage authors to place their materials (but not large datasets) in a Git repository hosted on a site such as GitHub, GitLab, or BitBucket. If the repository is private during the review process, please indicate the location where it will be available publicly upon publication, and also include the materials as a zip file (e.g., obtained directly from the Git hosting site) as supplementary materials.
-->

All code, output and data are recorded in the GitHub. Simulated data is
located in the `data` directory and may also be reproduced using the
scripts from `code/data`. The package rjpdmp should be installed as
described in the `code/README.md`. Scripts which run the methods on the
different datasets are located in the `code/Generate Output` directory
and their outputs are located in the `output` directory. Scripts to
generate figures 1-5 and the set of 4 tables from the paper can be found
in directory `code/Figures and Tables`. The figures and tables may be reproduced by simply running the code located in `code/Figures and Tables` they do not require re-simulating all results in the output directory.

