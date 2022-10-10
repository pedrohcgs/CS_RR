#---------------------------------------------------
#      Effect of Minimum Wage on Teen Employments
#      Codes based on Callaway and Sant'Anna (2021)
#      Main goal: show how we can do sensitivy analysis
#      for potential violations of parallel trends.
#      Combine Callaway and Sant'Anna (2021) with 
#      Rambachan and Roth (2022)
#---------------------------------------------------
# Install some packages
library(devtools)

install_github("bcallaway11/BMisc", dependencies = TRUE)
install_github("bcallaway11/did", dependencies = TRUE)
install_github("asheshrambachan/HonestDiD", dependencies = TRUE)
#-----------------------------------------------------------------------------
# Load packages
#-----------------------------------------------------------------------------
# Libraries
# Load libraries
library(ggplot2)
library(here)
library(foreign)
library(tidyverse)
library(dplyr)
library(did)
library(HonestDiD)
#---------------------------------------------------------------------------------------
# Aux Functions to integrate did and HonesDiD packages

## -----------------------------------------------------------------------------

#' @title honest_did
#'
#' @description a function to compute a sensitivity analysis
#'  using the approach of Rambachan and Roth (2022)
#' @param es an event study
honest_did <- function(es, ...) {
  UseMethod("honest_did", es)
}


#' @title honest_did.AGGTEobj
#'
#' @description a function to compute a sensitivity analysis
#'  using the approach of Rambachan and Roth (2022) when
#'  the event study is estimating using the `did` package
#'
#' @param e event time to compute the sensitivity analysis for.
#'  The default value is `e=0` corresponding to the "on impact"
#'  effect of participating in the treatment.
#' @param type Options are "smoothness" (which conducts a
#'  sensitivity analysis allowing for violations of linear trends
#'  in pre-treatment periods) or "relative_magnitude" (which
#'  conducts a sensitivity analysis based on the relative magnitudes
#'  of deviations from parallel trends in pre-treatment periods).
#' @inheritParams HonestDiD::createSensitivityResults
#' @inheritParams HonestDid::createSensitivityResults_relativeMagnitudes

honest_did.AGGTEobj <- function(es,
                                e=0,
                                type=c("smoothness", "relative_magnitude"),
                                method=NULL,
                                bound="deviation from parallel trends",
                                Mvec=NULL,
                                Mbarvec=NULL,
                                monotonicityDirection=NULL,
                                biasDirection=NULL,
                                alpha=0.05,
                                parallel=FALSE,
                                gridPoints=10^3,
                                grid.ub=NA,
                                grid.lb=NA,
                                ...) {
  
  type <- type[1]
  
  # make sure that user is passing in an event study
  if (es$type != "dynamic") {
    stop("need to pass in an event study")
  }
  
  # check if used universal base period and warn otherwise
  if (es$DIDparams$base_period != "universal") {
    warning("it is recommended to use a universal base period for honest_did")
  }
  
  # recover influence function for event study estimates
  es_inf_func <- es$inf.function$dynamic.inf.func.e
  
  # recover variance-covariance matrix
  n <- nrow(es_inf_func)
  V <- t(es_inf_func) %*% es_inf_func / (n*n) 
  
  
  nperiods <- nrow(V)
  npre <- sum(1*(es$egt < 0))
  npost <- nperiods - npre
  
  baseVec1 <- basisVector(index=(e+1),size=npost)
  
  orig_ci <- constructOriginalCS(betahat = es$att.egt,
                                 sigma = V, numPrePeriods = npre,
                                 numPostPeriods = npost,
                                 l_vec = baseVec1)
  type <- type[1]
  
  # make sure that user is passing in an event study
  if (es$type != "dynamic") {
    stop("need to pass in an event study")
  }
  
  # check if used universal base period and warn otherwise
  if (es$DIDparams$base_period != "universal") {
    stop("Use a universal base period for honest_did")
  }
  
  # recover influence function for event study estimates
  es_inf_func <- es$inf.function$dynamic.inf.func.e
  
  # recover variance-covariance matrix
  n <- nrow(es_inf_func)
  V <- t(es_inf_func) %*% es_inf_func / (n*n) 
  
  #Remove the coefficient normalized to zero
  referencePeriodIndex <- which(es$egt == -1)
  V <- V[-referencePeriodIndex,-referencePeriodIndex]
  beta <- es$att.egt[-referencePeriodIndex]
  
  nperiods <- nrow(V) 
  npre <- sum(1*(es$egt < -1))
  npost <- nperiods - npre
  
  baseVec1 <- basisVector(index=(e+1),size=npost)
  
  orig_ci <- constructOriginalCS(betahat = beta,
                                 sigma = V, 
                                 numPrePeriods = npre,
                                 numPostPeriods = npost,
                                 l_vec = baseVec1)
  
  if (type=="relative_magnitude") {
    if (is.null(method)) method <- "C-LF"
    robust_ci <- createSensitivityResults_relativeMagnitudes(betahat = beta, 
                                                             sigma = V, 
                                                             numPrePeriods = npre, 
                                                             numPostPeriods = npost,
                                                             bound=bound,
                                                             method=method,
                                                             l_vec = baseVec1,
                                                             Mbarvec = Mbarvec,
                                                             monotonicityDirection=monotonicityDirection,
                                                             biasDirection=biasDirection,
                                                             alpha=alpha,
                                                             gridPoints=100,
                                                             parallel=parallel)
    
  } else if (type=="smoothness") {
    robust_ci <- createSensitivityResults(betahat = beta,
                                          sigma = V, 
                                          numPrePeriods = npre, 
                                          numPostPeriods = npost,
                                          method=method,
                                          l_vec = baseVec1,
                                          monotonicityDirection=monotonicityDirection,
                                          biasDirection=biasDirection,
                                          alpha=alpha,
                                          parallel=parallel)
  }
  
  return(list(robust_ci=robust_ci, orig_ci=orig_ci, type=type))
}



##-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Load data used in Callaway and Sant'Anna (2021) application
min_wage <- readRDS((here("data",'min_wage_CS.rds')))

#--------------------------------------------------------------------------------------------
# Callaway and Sant'Anna (2021) procedure
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Formula for covarites
xformla <- ~ region + (medinc + pop ) + I(pop^2) + I(medinc^2)  + white + hs  + pov
#--------------------------------------------------------------------------------------------
# Using covariates and DR DiD with never-treated as comparison group
# Fix the reference time periods
CS_never_cond <- did::att_gt(yname="lemp",
                             tname="year",
                             idname="countyreal",
                             gname="first.treat",
                             #xformla=~1,
                             xformla = xformla,
                             control_group="nevertreated",
                             data = min_wage,
                             panel = TRUE,
                             base_period="universal",
                             bstrap = TRUE,
                             cband = TRUE)
# Now, compute event study
CS_es_never_cond <- aggte(CS_never_cond, type = "dynamic",
                          min_e = -5, max_e = 5)
summary(CS_es_never_cond)
# Plot event study
fig_CS <- ggdid(CS_es_never_cond,
      title = "Event-study aggregation \n DiD based on conditional PTA and using never-treated as comparison group ")
# Save event study plot
ggsave(here("plots","CS_es_nv_X_fixed_baseline.png"),
       fig_CS,  
       dpi = 500,
       width = 14, 
       height = 7)
#--------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------
# code for running honest_did
hd_cs_smooth_never <- honest_did(CS_es_never_cond,
                           type="smoothness")
hd_cs_smooth_never


hd_cs_rm_never <- honest_did(CS_es_never_cond, type="relative_magnitude")
hd_cs_rm_never
# Drop 0 as that is not really allowed.
hd_cs_rm_never$robust_ci <- hd_cs_rm_never$robust_ci[-1,]

## -----------------------------------------------------------------------------
# make sensitivity analysis plots
cs_HDiD_smooth <- createSensitivityPlot(hd_cs_smooth_never$robust_ci,
                      hd_cs_smooth_never$orig_ci)

cs_HDiD_relmag <- createSensitivityPlot_relativeMagnitudes(hd_cs_rm_never$robust_ci,
                                         hd_cs_rm_never$orig_ci)
## -----------------------------------------------------------------------------
# Save plots
ggsave(here("plots","cs_HDiD_smooth.png"),
       cs_HDiD_smooth,  
       dpi = 500,
       width = 14, 
       height = 7)


# Save plots
ggsave(here("plots","cs_HDiD_relmag.png"),
       cs_HDiD_relmag,  
       dpi = 500,
       width = 14, 
       height = 7)
