########################################################################################
# Run analysis in McMenamin et al. (2022) COVID-19 testing strategies in Hong Kong
#
# Author: Martina McMenamin
# Date: 2022-01-05
#######################################################################################

# Load libraries

library('pacman')
p_load(tidyverse, dplyr, viridis, ggExtra, ggpubr, scales, ggplot2, here, boot, cubature, distrEx)

# Load scripts

source('sens_pool.R')
source('expected_tests.R')
source('estimate_Rt.R')
source('outbreak_sim_paper.R')


