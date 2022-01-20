########################################################################################
# Run analysis in McMenamin et al. (2022) COVID-19 testing strategies in Hong Kong
#
# Author: Martina McMenamin
# Date: 2022-01-05
#######################################################################################

# Load libraries


pacman::p_load(
  tidyverse, dplyr, viridis, ggExtra, ggpubr, scales, ggplot2, here, boot, 
  cubature, distrEx, tvgeom, deSolve, EpiModel, gtable, grid, tidyr,
  reshape, graphics, hrbrthemes, incidence, remotes, zoo, lubridate, conflicted
  )

conflict_prefer("mutate", "dplyr")

# Load scripts

source(file.path("R/sens_pool.R"))
source(file.path('R/expected_tests.R'))
source(file.path('R/estimate_Rt.R'))
# source(file.path('R/outbreak_sim_paper.R'))
# source(file.path('R/outbreak_stoch.R'))


