# source("http://egret.psychol.cam.ac.uk/rlib/INSTALL_MANY_MODULES.R")

# Also, in Ubuntu: need modules listed in R_EXTRAS variable at http://egret.psychol.cam.ac.uk/techniques/scripts/rnc_ubuntu_setup

PREFERRED_CRAN_REPOSITORY = c(CRAN="http://cran.ma.imperial.ac.uk")

install_if_absent <- function(
  libname,
  repos=PREFERRED_CRAN_REPOSITORY,
  type=getOption("pkgType")
) {
  if (require(libname, character.only=TRUE)) {
    cat("Loaded library:", libname, "\n")
  } else {
    cat("Attempting to install package: ", libname, " (type=", type, ")\n", sep="")
    install.packages(libname, repos=repos, dependencies=TRUE, type=type)
    if (require(libname, character.only = TRUE)) {
      cat("Loaded newly installed library:", libname, "\n")
    } else {
      cat("FAILED TO INSTALL PACKAGE:", libname, "\n")
    }
  }
}

DESIRED_LIBRARIES = c(
  "arm",
  "afex",
  "arrayhelpers",
  # "base64",
  # "BMS",
  "bridgesampling",
  "broom",
  "Cairo",
  "car",
  "coda",
  "data.table",
  "diagram",
  "DiagrammeR",
  "doParallel",
  "doSNOW",
  # WHEN RCPP VERSION CATCHES UP # "dplyr",  # Oct 2014
  # "extrafont", # use Cairo instead
  "ez",
  "foreign",
  "gdata",
  "ggplot2",
  "gplots",
  "grid",
  "gridExtra",
  "gtools",
  "Hmisc",
  "inline", # for RStan
  "Kmisc",  # including cat.cb, to write to the clipboard
  "languageR", # removed from CRAN in 2013; then back 2013-12-12 (needs R 3.0.2) see http://cran.r-project.org/web/packages/languageR/index.html
  "lattice",
  "lme4",
  "LMERConvenienceFunctions",
  "lmerTest",
  "lsmeans",
  "ltm",
  "lubridate",
  "maptools",
  "matrixStats",
  "MCMCglmm",
  "nlme",
  "nortest",
  "parallel",
  "plyr",
  "plotrix",
  "popbio",
  "proto",
  # "R2jags",
  # "Rcmdr", # use RStudio instead,
  "Rcpp", # for RStan
  "raster",
  "readODS",
  "readr",
  "readxl",
  "reshape",
  "reshape2",
  "rgeos",
  "rgdal",
  # "rjags", # also needs some command-line prerequisites
  "rstan", # if fails: install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)
  # ... see https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
  "RMySQL",
  "RODBC",
  "semver",
  "shinystan",
  "snow",
  "sp",  # spatial, for maps
  "sqldf",
  "stringr",
  "TTR",
  "visreg",
  "XLConnect",
  "xlsx",
  "xtermStyle",
  "zoo",
  "MASS",
  "magrittr",
  "dplyr",
  "plyr",
  "boot",
  "sjlabelled",
  # folloings are the almost all packages related to bayes in r, maybe repeat with previous
  "abc",
  "abn",
  "AdMit",
  "arm",
  "AtelieR",
  "BaBooN",
  "BACCO",
  "BaM",
  "bamlss",
  "BART",
  "BAS",
  "BayesDA",
  "BayesFactor",
  "bayesGARCH",
  "bayesImageS",
  "bayesm",
  "bayesmeta",
  "bayesmix",
  "bayesQR",
  "BayesSummaryStatLM",
  "bayesSurv",
  "BayesTree",
  "BayesValidate",
  "BayesVarSel",
  "BayesX",
  "BayHaz",
  "BAYSTAR",
  "bbemkr",
  "BCBCSF",
  "BCE",
  "bclust",
  "bcp",
  "BDgraph",
  "BLR",
  "BMA",
  "Bmix",
  "bmixture",
  "BMS",
  "bnlearn",
  "boa",
  "Bolstad",
  "Boom",
  "BoomSpikeSlab",
  "bqtl",
  "bridgesampling",
  "brms",
  "bsamGP",
  "bspec",
  "bspmma",
  "BSquare",
  "bsts",
  "BVS",
  "catnet",
  "coalescentMCMC",
  "coda",
  "dclone",
  "deBInfer",
  "dlm",
  "DPpackage",
  "EbayesThresh",
  "ebdbNet",
  "eco",
  "eigenmodel",
  "ensembleBMA",
  "evdbayes",
  "exactLoglinTest",
  "factorQR",
  "FME",
  "geoR",
  "geoRglm",
  "ggmcmc",
  "gRain",
  "growcurves",
  "hbsae",
  "HI",
  "Hmisc",
  "iterLap",
  "LaplacesDemon",
  "LearnBayes",
  "lmm",
  "MasterBayes",
  "matchingMarkets",
  "mcmc",
  "MCMCglmm",
  "MCMCpack",
  "MCMCvis",
  "mgcv",
  "mlogitBMA",
  "MNP",
  "mombf",
  "monomvn",
  "NetworkChange",
  "nimble",
  "openEBGM",
  "pacbpred",
  "PAWL",
  "predmixcor",
  "PReMiuM",
  "prevalence",
  "profdpm",
  "pscl",
  "R2BayesX",
  "R2jags",
  "R2WinBUGS",
  "ramps",
  "revdbayes",
  "RJaCGH",
  "rjags",
  "RSGHB",
  "RSGHB",
  "rstan",
  "rstiefel",
  "runjags",
  "Runuran",
  "RxCEcolInf",
  "SamplerCompare",
  "SampleSizeMeans",
  "SampleSizeProportions",
  "sbgcop",
  "SimpleTable",
  "sna",
  "spBayes",
  "spikeslab",
  "spikeSlabGAM",
  "spTimer",
  "ssgraph",
  "stochvol",
  "tgp",
  "zic",
  "rstanarm",
  "bayesplot",
  "shinystan",
  "rstan",
  "shiny",
  "shinystan",
  "parallel",
  "bridgesampling",
  "bayesplot",
  "tinytex",
  "survminer",
  "precrec",
  "PRROC",
  "moments",
  "devtools",
  "forecast",
  "MARSS",
  "vars",
  "metafor"
)

# First, update packages:
update.packages(checkBuilt=TRUE, ask=FALSE)

# Now install:
for (libname in DESIRED_LIBRARIES) {
  install_if_absent(libname)
}

# rstan used to require special procedures, but doesn't any more:
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

#cat("
#
#===============================================================================
#DONE MAJOR PACKAGES. Now trying RStan...
#===============================================================================
#")
#
## RStan is a package that's a bit different:
#install_if_absent(
#    "rstan",
#    repos="https://github.com/stan-dev/rstan",
#    type="source")
## ... formerly (1) http://wiki.stan.googlecode.com/git/R
## ... formerly (2) http://wiki.rstan-repo.googlecode.com/git/

library(devtools)
install_github("kassambara/factoextra")
