#include <RcppArmadillo.h>
using namespace Rcpp;
#include "modelprior.h"
#include "likelyhood.h"
#include "pdmp_ss.h"
#include "zigzag_logit.h"
#include "zigzag_rr.h"
#include "bps_logit.h"
#include "bps_rr.h"
#include "hmc_logit.h"
