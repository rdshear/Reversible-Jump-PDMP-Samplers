## Code for simulation studies

Code for the simulation studies was implemented in R and Rcpp. The package developed for implementing the novel methods may be found in this directory and can be installed via the commands `devtools::load_all("code/rjpdmp")` followed by `devtools::install("code/rjpdmp")` (later versions of the package will be available on CRAN). Simulations were implemented using a high performance computing cluster -- the job scripts used are given in directory `code\Job Submission Scripts`. Each method was run for a fixed computational budget and the number of iterations and final estimates were recorded (see output directory). Scripts for running the methods to generate the output may be found in the directory `code\Generate Output`. Scripts for processing the output and generating tables and figures may be found in directory `code\Figures and Tables`. 


