# Author: Kevin See
# Purpose: prep and run DABOM
# Created: 4/1/20
# Last Modified: 11/7/2022
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(DABOM)
library(tidyverse)
library(rjags)
library(magrittr)
library(lubridate)
library(here)

#-----------------------------------------------------------------
# load configuration and site_df data
load(here('analysis/data/derived_data/site_config.rda'))

#-----------------------------------------------------------------
# Load required DABOM data
#-----------------------------------------------------------------
# set year
yr = 2022

# for(yr in 2011:2019) {
#   cat(paste("Working on", yr, "\n\n"))

# # load and filter biological data
# bio_df = read_rds(here('analysis/data/derived_data',
#                        'Bio_Data_2011_2021.rds')) %>%
#   filter(year == yr)

# load processed detection histories
load(here('analysis/data/derived_data/PITcleanr',
          paste0('UC_Steelhead_', yr, '.rda')))

# filter to keep only the observations you want to keep
filter_obs = prepped_ch %>%
  mutate(user_keep_obs = if_else(is.na(user_keep_obs),
                                 auto_keep_obs,
                                 user_keep_obs)) %>%
  filter(user_keep_obs)


# determine origin of each fish
fish_origin = bio_df %>%
  filter(tag_code %in% unique(filter_obs$tag_code)) %>%
  select(tag_code, origin) %>%
  distinct()

# file path to the default and initial model
basic_modNm = here('analysis/model_files', "PRA_DABOM.txt")

writeDABOM(file_name = basic_modNm,
           parent_child = parent_child,
           configuration = configuration,
           time_varying = F)

#------------------------------------------------------------------------------
# Alter default model code for species and year of
# interest; sets prior for some detection node efficiencies at 0 or 100%
# based on actual tag detection data; 0% if no tags were seen
#------------------------------------------------------------------------------

# filepath for specific JAGS model code for species and year
mod_path = here('analysis/model_files',
                paste0('PRA_Steelhead_', yr, '.txt'))

# writes species and year specific jags code
fixNoFishNodes(init_file = basic_modNm,
               file_name = mod_path,
               filter_ch = filter_obs,
               parent_child = parent_child,
               configuration = configuration,
               fish_origin = fish_origin)


setInitialValues = function(filter_ch,
                            parent_child,
                            configuration) {
  
  stopifnot(exprs = {
    !is.null(filter_ch)
    !is.null(parent_child)
    !is.null(configuration)
  })
  
  # how many child sites does each parent site have?
  parent_info = parent_child %>%
    group_by(parent, parent_rkm) %>%
    mutate(n_child = n_distinct(child))
  
  # determine parent site, and what branch number the tag must have taken
  branch_df = getNodeInfo(parent_child,
                          configuration) %>%
    select(site_code, parent_site, child_num) %>%
    distinct()
  
  # build node order
  no = parent_child %>%
    PITcleanr::addParentChildNodes(configuration = configuration) %>%
    PITcleanr::buildNodeOrder()
  
  # how many branches at each branching node?
  n_branch_list = setBranchNums(parent_child) %>%
    # add a black box
    map(.f = function(x) x + 1)
  
  # look at estimated spawn location, and the sites tag must have crossed to get there
  spawn_node = estimateFinalLoc(filter_ch) %>%
    select(tag_code, final_node) %>%
    distinct() %>%
    mutate(spawn_site = if_else(grepl("B0$", final_node) &
                                  nchar(final_node) >= 5,
                                str_remove(final_node, "B0"),
                                final_node),
           spawn_site = if_else(grepl("A0$", spawn_site) &
                                  nchar(spawn_site) >= 5,
                                str_remove(spawn_site, "A0"),
                                spawn_site)) %>%
    left_join(no %>%
                select(final_node = node,
                       spawn_path = path),
              by = "final_node") %>%
    separate_rows(spawn_path) %>%
    rename(node = spawn_path) %>%
    left_join(no %>%
                select(node,
                       node_order),
              by = "node") %>%
    mutate(site_code = if_else(grepl("B0$", node) &
                                 nchar(node) >= 5,
                               str_remove(node, "B0"),
                               node),
           site_code = if_else(grepl("A0$", site_code) &
                                 nchar(site_code) >= 5,
                               str_remove(site_code, "A0"),
                               site_code)) %>%
    arrange(tag_code, node_order)
  
  
  # each tag passed each of these sites
  tag_sites = spawn_node %>%
    select(tag_code,
           spawn_site,
           site_code) %>%
    distinct() %>%
    group_by(tag_code) %>%
    mutate(lead_site = lead(site_code)) %>%
    ungroup() %>%
    left_join(branch_df %>%
                rename(lead_site = site_code,
                       site_code = parent_site),
              by = c("site_code", "lead_site")) %>%
    left_join(n_branch_list %>%
                unlist() %>%
                enframe(name = "site_code",
                        value = "max_branch"),
              by = "site_code") %>%
    mutate(child_num = if_else(is.na(child_num) & !is.na(max_branch),
                               as.integer(max_branch),
                               child_num))
  
  # construct some initial values lists
  a_list = tag_sites %>%
    filter(site_code %in% names(n_branch_list)) %>%
    split(list(.$site_code)) %>%
    map(.f = function(x) {
      not_there = max(x$max_branch) + 1
      
      x %>%
        complete(tag_code = unique(filter_ch$tag_code),
                 fill = list(child_num = not_there)) %>%
        arrange(tag_code) %>%
        select(tag_code, site_code, child_num) %>%
        pull(child_num)
    }) %>%
    rlang::set_names(nm = function(x) paste0("a_", x))
  
  # any eta initial values needed?
  n_eta <- parent_child %>%
    dplyr::count(parent,
                 name = "n_child") %>%
    filter(n_child == 1) %>%
    nrow()
  
  if(n_eta > 0) {
    eta_list = parent_child %>%
      dplyr::count(parent,
                   name = "n_child") %>%
      filter(n_child == 1) %>%
      select(parent) %>%
      split(list(.$parent)) %>%
      map(.f = function(x) {
        tag_sites %>%
          filter(site_code == x$parent) %>%
          mutate(seen = if_else(!is.na(lead_site),
                                1, 0)) %>%
          complete(tag_code = unique(filter_ch$tag_code),
                   fill = list(seen = 0)) %>%
          arrange(tag_code) %>%
          pull(seen)
      }) %>%
      rlang::set_names(nm = function(x) paste0("eta_", x))
  } else {
    eta_list = NULL
  }
  
  jags_inits <- function() {
    y = c(a_list,
          eta_list)
    return(y)
  }
  
  return(jags_inits)
}
#------------------------------------------------------------------------------
# Creates a function to spit out initial values for MCMC chains
init_fnc = setInitialValues(filter_obs,
                            parent_child,
                            configuration)

createJAGSinputs = function(filter_ch = NULL,
                            parent_child = NULL,
                            configuration = NULL,
                            fish_origin = NULL) {
  
  stopifnot(exprs = {
    !is.null(filter_ch)
    !is.null(parent_child)
    !is.null(configuration)
  })
  
  if(is.null(fish_origin)) {
    fish_origin = filter_ch %>%
      select(tag_code) %>%
      distinct() %>%
      mutate(origin = "W")
  }
  
  
  # determine starting point (root_site)
  root_site = PITcleanr::buildNodeOrder(parent_child) %>%
    filter(node_order == 1) %>%
    pull(node)
  
  # get the column names of the capture history matrix
  col_nms = defineDabomColNms(root_site = root_site,
                              parent_child = parent_child,
                              configuration = configuration) %>%
    unlist() %>%
    as.vector()
  
  # create capture history
  cap_hist = createDABOMcapHist(filter_ch = filter_ch,
                                parent_child = parent_child,
                                configuration = configuration,
                                split_matrices = F)
  
  # what kind of fish (wild or hatchery)
  fish_type = cap_hist %>%
    select(tag_code) %>%
    left_join(fish_origin %>%
                mutate(origin = recode(origin,
                                       "W" = 1,
                                       "H" = 2)),
              by = "tag_code") %>%
    pull(origin)
  
  # how many branches at each branching node?
  n_branch_list = setBranchNums(parent_child) %>%
    rlang::set_names(nm = function(x) paste0("n_branch_", x)) %>%
    # add a black box
    map(.f = function(x) x + 1)
  
  # set Dirichlet vectors
  init_val_func = setInitialValues(filter_ch = filter_ch,
                                   parent_child = parent_child,
                                   configuration = configuration)
  init_mats = init_val_func()
  
  # init_mats[stringr::str_remove(names(init_mats), '^a_') %in% stringr::str_remove(names(n_branch_list), 'n_branch_')]
  
  dirich_df = n_branch_list %>%
    unlist() %>%
    tibble::enframe(name = "site",
                    value = "n_brnch") %>%
    mutate(site = stringr::str_remove(site, 'n_branch_')) %>%
    left_join(init_mats %>%
                enframe(name = 'site',
                        value = 'inits') %>%
                mutate(site = stringr::str_remove(site, '^a_')),
              by = 'site') %>%
    mutate(dirch_vec = purrr::map2(n_brnch,
                                   inits,
                                   .f = function(x, y) {
                                     c(createDirichletVector(x,
                                                             table(y[fish_type == 1]),
                                                             initial_one = F,
                                                             final_one = T),
                                       createDirichletVector(x,
                                                             table(y[fish_type == 2]),
                                                             initial_one = F,
                                                             final_one = T)) %>%
                                       matrix(nrow = 2,
                                              byrow = T)
                                   }))
  
  dirich_vecs = dirich_df$dirch_vec %>%
    rlang::set_names(paste0(dirich_df$site, '_dirch_vec'))
  
  
  jags_list = c(list(n_fish = nrow(cap_hist),
                     # vector of zeros, large enough to match any element of dabom_list
                     zero_vec = rep(0, max(unlist(n_branch_list)) + 1),
                     cap_hist = cap_hist %>%
                       select(any_of(col_nms)),
                     fish_type = fish_type),
                n_branch_list,
                dirich_vecs)
  
  return(jags_list)
  
  
}
# Create all the input data for the JAGS model
jags_data = createJAGSinputs(filter_ch = filter_obs,
                             parent_child = parent_child,
                             configuration = configuration,
                             fish_origin = fish_origin)

# Tell JAGS which parameters in the model that it should save.
jags_params = setSavedParams(model_file = mod_path,
                             time_varying = F)


# run the model
jags = jags.model(mod_path,
                  data = jags_data,
                  inits = init_fnc,
                  # n.chains = 1,
                  # n.adapt = 5)
                  n.chains = 4,
                  n.adapt = 5000)


#--------------------------------------
# test the MCMC outcome and summary functions
dabom_mod = coda.samples(jags,
                         jags_params,
                         # n.iter = 10)
                         n.iter = 5000,
                         thin = 10)


save(dabom_mod, jags_data, filter_obs, bio_df,
     file = here("analysis/data/derived_data/model_fits",
                 paste0('PRA_DABOM_Steelhead_', yr,'.rda')))

# rm(dabom_mod, jags_data, filter_obs)
# }

#------------------------------------------------------------------------------
# diagnostics
#------------------------------------------------------------------------------
# load model run
load(here("analysis/data/derived_data/model_fits",
          paste0('PRA_DABOM_Steelhead_', yr,'.rda')))

# using mcmcr package
library(mcmcr)

# pull out mcmc.list object
my_mod = dabom_mod

#---------------------------------------
# using mcmcr
anyNA(my_mod)
my_mcmcr = as.mcmcr(my_mod)

# get Rhat statistics for all parameters
rhat_df = rhat(my_mcmcr,
               by = 'parameter',
               as_df = T) %>%
  full_join(esr(my_mcmcr,
                by = 'parameter',
                as_df = T)) %>%
  mutate(type = if_else(grepl('_p$', parameter),
                        'Detection',
                        if_else(grepl('^psi', parameter) |
                                  grepl('^phi', parameter),
                                'Movement',
                                'Other')))

# plot histogram of Rhat statistics
rhat_df %>%
  ggplot(aes(x = rhat)) +
  geom_histogram(fill = 'blue',
                 # bins = 40) +
                 binwidth = 0.001) +
  facet_wrap(~ type,
             scales = 'free')

# which parameters have converged and which haven't?
convg_df = converged(my_mcmcr,
                     by = 'parameter',
                     as_df = T)

janitor::tabyl(convg_df,
               converged)

# look at parameters that have not converged
convg_df %>%
  # filter(!converged) %>%
  left_join(rhat_df) %>%
  arrange(esr)

#---------------------------------------
# using postpack
library(postpack)

# what parameters were tracked?
get_params(my_mod,
           type = 'base_only')

# some summary statistics
post_summ(my_mod,
          '_p$') %>%
  t() %>%
  as_tibble(rownames = 'param')

post_summ(my_mod,
          '^phi') %>%
  t() %>%
  as_tibble(rownames = 'param') %>%
  filter(mean > 0)

post_summ(my_mod,
          '^psi') %>%
  t() %>%
  as_tibble(rownames = 'param') %>%
  filter(mean > 0)




param_chk = c('psi_RRF')
param_chk = convg_df %>%
  filter(!converged) %>%
  pull(parameter)

post_summ(my_mod,
          param_chk) %>%
  t() %>%
  as_tibble(rownames = 'param') %>%
  mutate(cv = sd / mean) %>%
  arrange(desc(cv))


diag_plots(post = my_mod,
           p = param_chk,
           save = F,
           file = here('outgoing/PRA_diagnostics.pdf'))

# calculate Brooks-Gelman-Rubin Potential Scale Reduction Factor (Rhat)
# if ratio is close to 1, the chains have converged to the same distribution
# <1.10 is generally considered converged
post_summ(my_mod,
          # '_p$',
          get_params(my_mod,
                     type = 'base_only'),
          neff = T, # effective sample size
          Rhat = T)[c("Rhat", "neff"),] %>%
  t() %>%
  as_tibble(rownames = 'param') %>%
  filter(!is.na(Rhat)) %>%
  arrange(neff)

# find and remove params where Rhat == "NaN"
all_params = get_params(my_mod,
                        type = 'base_only')


post_summ_nas = post_summ(my_mod,
                          # '_p$',
                          all_params,
                          neff = T, # effective sample size
                          Rhat = T)[c("Rhat", "neff"),] %>%
  t() %>%
  as.data.frame() %>%
  as_tibble(rownames = 'param') %>%
  filter(Rhat == "NaN") %>%
  pull(param)

param_chk = get_params(my_mod,
                       type = 'base_only')[grep('_p$', get_params(my_mod, type = 'base_only'))]
param_chk = param_chk[!param_chk %in% post_summ_nas]

# diagnostic plots for remaining params
diag_plots(post = my_mod,
           p = param_chk,
           ext_device = T)

# save plots
diag_plots(post = my_mod,
           p = param_chk,
           save = T,
           file = here('../output/DABOM_trace_plots.pdf'))

