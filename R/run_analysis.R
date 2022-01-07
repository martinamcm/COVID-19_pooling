########################################################################################
# Run analysis in McMenamin et al. (2022) COVID-19 testing strategies in Hong Kong
#
# Author: Martina McMenamin
# Date: 2022-01-05
#######################################################################################

# Load libraries


pacman::p_load(
  tidyverse, dplyr, viridis, ggExtra, ggpubr, scales, ggplot2, here, boot, 
  cubature, distrEx, tvgeom, deSolve, EpiModel, gtable, grid
  )

# Load scripts

source(here('R', 'sens_pool.R'))
source(here('R', 'expected_tests.R'))
source(here('R', 'estimate_Rt.R'))
source(here('R', 'outbreak_sim_paper.R'))


