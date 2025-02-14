Untitled
================
2022-06-22

<!--HOW TO COMPLETE THIS FORM:-->
<!--
1. Checkboxes in this document appear as follows: 

- [ ] This is a checkbox 

To check a checkbox, replace [ ] by [x], as follows: 

- [x] This is a checked checkbox 

Note that older versions of RStudio (versions lower than 1.3) may not create a formatted checkbox but will leave the original characters, i.e., literally "[ ]" or "[x]". It's fine to submit a PDF in this form.
 
2. For text answers, simply type the relevant text in the areas indicated. A blank line starts a new paragraph. 
 
3. Comments (like these instructions) provide additional instructions throughout the form. There is no need to remove them; they will not appear in the compiled document. 

4. If you are comfortable with Markdown syntax, you may choose to include any Markdown-compliant formatting in the form. For example, you may wish to include R code chunks and compile this document in R Markdown.
-->

This form documents the artifacts associated with the article (i.e., the
data and code supporting the computational findings) and describes how
to reproduce the findings.

# Part 1: Data

-   [x] This paper does not involve analysis of external data (i.e., no
    data are used or the only data are generated by the authors via
    simulation in their code).

<!--
If box above is checked and if no simulated/synthetic data files are provided by the authors, please skip directly to the Code section. Otherwise, continue.
-->

-   [ ] I certify that the author(s) of the manuscript have legitimate
    access to and permission to use the data used in this manuscript.

<!-- If data are simulated using random number generation, please be sure to set the random number seed in the code you provide -->

## Abstract

<!--
Provide a short (< 100 words), high-level description of the data
-->

The data used in this paper may be reproduced via the data simulation
scripts in the code folder or may be read directly from the R data files
under the data directory. Details about the scripts corresponding to the
generated data are given in the Readme section of the GitHub in the
`data` directory.

## Availability

-   [x] Data **are** publicly available.
-   [ ] Data **cannot be made** publicly available.

If the data are publicly available, see the *Publicly available data*
section. Otherwise, see the *Non-publicly available data* section,
below.

### Publicly available data

-   [ ] Data are available online at:

-   [x] Data are available as part of the paper’s supplementary
    material.

-   [ ] Data are publicly available by request, following the process
    described here:

-   [ ] Data are or will be made available through some other mechanism,
    described here:

<!-- If data are available by request to the authors or some other data owner, please make sure to explain the process of requesting access to the data. -->

### Non-publicly available data

<!--
The Journal of the American Statistical Association requires authors to make data accompanying their papers available to the scientific community except in cases where: 1) public sharing of data would be impossible, 2) suitable synthetic data are provided which allow the main analyses to be replicated (recognizing that results may differ from the "real" data analyses), and 3) the scientific value of the results and methods outweigh the lack of reproducibility.

Please discuss the lack of publicly available data. For example:
-   why data sharing is not possible,
-   what synthetic data are provided, and 
-   why the value of the paper's scientific contribution outweighs the lack of reproducibility.
-->

## Description

### File format(s)

<!--
Check all that apply
-->

-   [x] CSV or other plain text.
-   [x] Software-specific binary format (.Rda, Python pickle, etc.):
    pkcle
-   [ ] Standardized binary format (e.g., netCDF, HDF5, etc.):
-   [ ] Other (please specify):

### Data dictionary

<!--
A data dictionary provides information that allows users to understand the meaning, format, and use of the data.
-->

-   [x] Provided by authors in the following file(s): `data/README.md`
    from the GitHub and Section 5 of the manuscript
-   [ ] Data file(s) is(are) self-describing (e.g., netCDF files)
-   [ ] Available at the following URL:

### Additional Information (optional)

<!-- 
OPTIONAL: Provide any additional details that would be helpful in understanding the data. If relevant, please provide unique identifier/DOI/version information and/or license/terms of use.
-->

# Part 2: Code

## Abstract

The package developed for implementing the novel methods may be found in
the `code` directory and can be installed via the command
`R CMD INSTALL rjpdmp_2.0.0.tar.gz` (later versions of the package will
be available on CRAN). Simulations were implemented using a high
performance computing cluster see README in `code` directory of GitHub
for details.

<!--
Provide a short (< 100 words), high-level description of the code. If necessary, more details can be provided in files that accompany the code. If no code is provided, please state this and say why (e.g., if the paper contains no computational work).
-->

## Description

### Code format(s)

<!--
Check all that apply
-->

-   [x] Script files
    -   [x] R
    -   [ ] Python
    -   [ ] Matlab
    -   [ ] Other:
-   [x] Package
    -   [x] R
    -   [ ] Python
    -   [ ] MATLAB toolbox
    -   [ ] Other:
-   [ ] Reproducible report
    -   [ ] R Markdown
    -   [ ] Jupyter notebook
    -   [ ] Other:
-   [x] Shell script
-   [ ] Other (please specify):

### Supporting software requirements

#### Version of primary software used

Due to version conflict with the R package rstan on the computing
cluster we used two versions of R for the experiments. R version 3.4.1
was used for the logistic regression examples and R version 3.6.2 was
used for the robust regression examples (where stan was a competitor).
Details on the R version number of hours, memory and number of parallel
jobs are given in the `code\Job Submission Scripts` directory.

#### Libraries and dependencies used by the code

The R package `rjpdmp` created to implement the methods is given in the
`code` directory. Additional packages and dependencies include:

1.  Rcpp_1.0.6 (for R 3.4.1) and Rcpp_1.0.6 (for R 3.6.2)
2.  data.table_1.12.8
3.  ggplot2_3.3.5
4.  gridExtra_2.3  
5.  Rmisc_1.5.1  
6.  rstan_2.19.2
7.  nimble_0.10.1
8.  devtools_2.4.2
9.  RcppArmadillo_0.10.2.1.0

<!--
Include version numbers (e.g., version numbers for any R or Python packages used)
-->

### Supporting system/hardware requirements (optional)

Details about the hardware used in the computing cluser is given here:
<https://cms.qut.edu.au/__data/assets/pdf_file/0012/388785/high-performance-computing.pdf>.
Operating system is SUSE Linux 3780 x 64bit Intel Xeon Cores. <!--
OPTIONAL: System/hardware requirements including operating system with version number, access to cluster, GPUs, etc.
-->

### Parallelization used

-   [ ] No parallel code used
-   [ ] Multi-core parallelization on a single machine/node
    -   Number of cores used:
-   [x] Multi-machine/multi-node parallelization
    -   Number of nodes and cores used: 60 cores were used. See
        `code\Job Submission Scripts` for cluser submission files

### License

-   [x] MIT License (default)
-   [ ] BSD
-   [ ] GPL v3.0
-   [ ] Creative Commons
-   [ ] Other: (please specify)

### Additional information (optional)

<!--
OPTIONAL: By default, submitted code will be published on the JASA GitHub repository (http://github.com/JASA-ACS) as well as in the supplementary material. Authors are encouraged to also make their code available in a public code repository, such as on GitHub, GitLab, or BitBucket. If relevant, please provide unique identifier/DOI/version information (e.g., a Git commit ID, branch, release, or tag). If the code and workflow are provided together, this section may be omitted, with information provided in the "Location" section below.
-->

# Part 3: Reproducibility workflow

<!--
The materials provided should provide a straightforward way for reviewers and readers to reproduce analyses with as few steps as possible. 
-->

## Scope

The provided workflow reproduces:

-   [x] Any numbers provided in text in the paper
-   [x] The computational method(s) presented in the paper (i.e., code
    is provided that implements the method(s))
-   [x] All tables and figures in the paper
-   [ ] Selected tables and figures in the paper, as explained and
    justified below:

## Workflow

### Location

The workflow is available:

<!--
Check all that apply, and in the case of a Git repository include unique identifier, such as specific commit ID, branch, release, or tag.
-->

-   [ ] As part of the paper’s supplementary material.
-   [x] In this Git repository:
-   [ ] Other (please specify):

All code, output and data are recorded in the GitHub. Simulated data is
located in the `data` directory and may also be reproduced using the
scripts from `code/data`. The package rjpdmp should be installed as
described in the `code/README.md`. Scripts which the run methods on the
different datasets are located in the `code/Generate Output` directory
and their outputs are located in the `output` directory. Scripts to
generate figures 1-5 and the set of 4 tables from the paper can be found
in directory `code/Figures and Tables`.

<!--
Indicate where the materials (generally including the code, unless in a separate location and indicated in the previous section) are available. We strongly encourage authors to place their materials (but not large datasets) in a Git repository hosted on a site such as GitHub, GitLab, or BitBucket. If the repository is private during the review process, please indicate the location where it will be available publicly upon publication, and also include the materials as a zip file (e.g., obtained directly from the Git hosting site) as supplementary materials.
-->

### Format(s)

<!--
Check all that apply
-->

-   [ ] Single master code file
-   [ ] Wrapper (shell) script(s)
-   [ ] Self-contained R Markdown file, Jupyter notebook, or other
    literate programming approach
-   [x] Text file (e.g., a readme-style file) that documents workflow
-   [ ] Makefile
-   [ ] Other (more detail in *Instructions* below)

### Instructions

<!--
Describe how to use the materials provided to reproduce analyses in the manuscript. Additional details can be provided in file(s) accompanying the reproducibility materials. If no workflow is provided, please state this and say why (e.g., if the paper contains no computational work).
-->

### Expected run-time

Approximate time needed to reproduce the analyses on a standard desktop
machine:

-   [ ] \< 1 minute
-   [ ] 1-10 minutes
-   [ ] 10-60 minutes
-   [ ] 1-8 hours
-   [ ] \> 8 hours
-   [x] Not feasible to run on a desktop machine, as described here:

Full Monte Carlo simulation studies involved running at least 4 to 6
methods for 2 minutes 100 times on over 60 datasets. The total
computation time was well over 1,200 hours and implementation was
carried out using a computing cluster. A single run of the methods on a
single dataset may be reproduced using the recorded number of iterations
from the cluster run and the same seed for a much lower computational
cost (See the *Additional information* section below for an example).
The output from the cluster has been included and the figures and tables
may be reproduced using the relevant files from the
`code\Figures and Tables` directory.

### Additional information (optional)

As the methods were implemented using a computing cluster full
reproduction may be computationally challenging. The code below
demonstrates an example of replicating the results for a single dataset
with the same seed and computational load used on the cluster.

``` r
library(rjpdmp)
## Prior used in the application
rj_val <- 0.6
prior_on <- 10/p
prior_sigma2 <- 10

maxT <- 200
maxI <- 10^6

## Data set to reproduce
n = 100; p = 100; sim =1
beta0 <- rep(0,p)
gamma_0 <- rep(0,p)

name <- paste0("output/Logistic/replication_logit_n%d_p%d_sim_%d.RData")
name <- sprintf(name,n,p,sim)
load(name)
datname <- paste0("data/data_Logit_n%d_p%d_sim_%d.RData")
datname <- sprintf(datname,n,p,sim)
load(datname)

## Example replication from the Zig-Zag sampler for a given seed
seed <- 3
set.seed(seed);system.time(zigzag_fit <- zigzag_logit(nmax = nIters[3,seed], dataX = data$dataX, datay = data$dataY,
                                                     prior_sigma2 = prior_sigma2,theta0 = gamma_0,
                                                     x0 = beta0, rj_val = rj_val,ppi = prior_on,
                                                     maxTime = maxT))

## Check reproducibility of marginal probability of inculsion
mp <- model_probabilities(zigzag_fit$times,zigzag_fit$positions, marginals = 1:p, 
                    burnin = ceiling(0.1*nIters[3,seed]))$marginal_prob
mp - results[,1,3,seed]

## Check reproducibility of marginal means
mm <- marginal_mean(zigzag_fit$times,zigzag_fit$positions,zigzag_fit$theta, marginals = 1:p,
                    burnin = ceiling(0.1*nIters[3,seed]))
mm - results[,2,3,seed] 

## Example replication of Gibbs sampler for a given seed
set.seed(seed);system.time(gibbs_fit <- gibbs_logit(datay = data$dataY, dataX = data$dataX, beta = beta0, 
                                                   gamma = gamma_0,ppi = prior_on, prior_sigma2 = prior_sigma2,
                                                   nsamples = nIters[1,seed], maxTime = maxT))
## Check reproducible marginal probability of inculsion
mp <- rowMeans(gibbs_fit$gamma[,-c(1:ceiling(0.1*nIters[1,seed]))])
mp - results[,1,1,seed]

## Check reproducibility of marginal means
mm <- rowMeans(gibbs_fit$beta[,-c(1:ceiling(0.1*nIters[1,seed]))]*gibbs_fit$gamma[,-c(1:ceiling(0.1*nIters[1,seed]))])
mm - results[,2,1,seed] 
```

<!--
OPTIONAL: Additional documentation provided (e.g., R package vignettes, demos or other examples) that show how to use the provided code/software in other settings.
-->

# Notes (optional)

<!--
OPTIONAL: Any other relevant information not covered on this form. If reproducibility materials are not publicly available at the time of submission, please provide information here on how the reviewers can view the materials.
-->
