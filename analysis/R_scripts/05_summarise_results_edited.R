for(yr in 2011:2022) {
  
  cat(paste("Working on", yr, "\n\n"))
  
  #   if(yr %in% c(2011:2015, 2018)) {
  #     dam_cnt_name == "PriestRapids"
  #   } else {
  #     dam_cnt_name == "RockIsland"
  #   }
  
  
  # load compressed detections and biological data
  load(here('analysis/data/derived_data/PITcleanr',
            paste0('UC_Steelhead_', yr, '.rda')))
  
  # load JAGS MCMC results
  load(here("analysis/data/derived_data/model_fits",
            paste0('PRA_DABOM_Steelhead_', yr,'.rda')))
  
  
  # estimate final spawning location
  tag_summ = summarizeTagData(filter_obs,
                              bio_df %>%
                                filter(tag_code %in% unique(filter_obs$tag_code)) %>%
                                group_by(tag_code) %>%
                                slice(1) %>%
                                ungroup())
  
  # look at which branch each tag was assigned to for spawning
  brnch_df = buildNodeOrder(addParentChildNodes(parent_child, configuration)) %>%
    separate(col = path,
             into = paste("step", 1:max(.$node_order), sep = "_"),
             remove = F) %>%
    mutate(group = if_else(node == "PRA",
                           "Start",
                           if_else(grepl('LWE', path) | node %in% c("CLK"),
                                   "Wenatchee",
                                   if_else(grepl("ENL", path),
                                           "Entiat",
                                           if_else(grepl("LMR", path),
                                                   "Methow",
                                                   if_else(grepl("OKL", path) | node %in% c("FST"),
                                                           "Okanogan",
                                                           if_else(step_2 != "RIA" & !is.na(step_2),
                                                                   "BelowPriest",
                                                                   if_else(node == "WEA",
                                                                           "WellsPool",
                                                                           "Other")))))))) %>%
    select(-starts_with("step")) %>%
    mutate(group = factor(group,
                          levels = c("Wenatchee",
                                     "Entiat",
                                     "Methow",
                                     "Okanogan",
                                     "BelowPriest",
                                     "WellsPool",
                                     "Start",
                                     "Other")))
  
  tag_summ %<>%
    left_join(brnch_df,
              by = c("final_node" = "node"))
  
  
  
  # summarize detection probabilities
  detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                     filter_ch = filter_obs)
  
  # which sites had detection probabilities fixed at 0% or 100%
  detect_summ %>%
    filter(sd == 0)
  
  # look at all the other sites
  detect_summ %>%
    filter(sd > 0) %>%
    # arrange(desc(n_tags))
    arrange(desc(sd))
  
  # compile all movement probabilities, and multiply them appropriately
  trans_df = compileTransProbs_PRA(dabom_mod,
                                   parent_child) %>%
    mutate(origin = recode(origin,
                           "2" = "H",
                           "1" = "W"))
  
  # summarize transition probabilities
  trans_summ = trans_df %>%
    group_by(origin, param) %>%
    summarise(mean = mean(value),
              median = median(value),
              mode = estMode(value),
              sd = sd(value),
              skew = moments::skewness(value),
              kurtosis = moments::kurtosis(value),
              lowerCI = coda::HPDinterval(coda::as.mcmc(value))[,1],
              upperCI = coda::HPDinterval(coda::as.mcmc(value))[,2],
              .groups = "drop") %>%
    mutate(across(c(mean, median, mode, sd, matches('CI$')),
                  ~ if_else(. < 0, 0, .)))
  
  #-----------------------------------------------------------------
  # total escapement past Priest, by origin
  start_date = paste0(yr-1, '0601')
  end_date = paste0(yr, '0531')
  
  
  # # start with PIT-tag based reascension data
  org_escape = queryPITtagData(damPIT = 'PRA',
                               spp = "Steelhead",
                               start_date = start_date,
                               end_date = end_date) %>%
    filter(!str_detect(TagId, "000.0")) %>%
    mutate(SpawnYear = yr) %>%
    mutate(across(TagIdAscentCount,
                  tidyr::replace_na,
                  0)) %>%
    mutate(ReAscent = ifelse(TagIdAscentCount > 1, T, F)) %>%
    group_by(Species, SpawnYear, Date) %>%
    summarise(tot_tags = n_distinct(TagId),
              reascent_tags = n_distinct(TagId[ReAscent]),
              .groups = "drop") %>%
    group_by(Species, SpawnYear) %>%
    summarise(across(matches('tags'),
                     sum,
                     na.rm = T),
              .groups = "drop") %>%
    mutate(reasc_rate = reascent_tags / tot_tags,
           reasc_rate_se = sqrt(reasc_rate * (1 - reasc_rate) / tot_tags)) %>%
    # add window counts
    bind_cols(getWindowCounts(dam = 'PRD',
                              spp = "Steelhead",
                              start_date = start_date,
                              end_date = end_date) %>%
                summarise_at(vars(win_cnt),
                             list(sum),
                             na.rm = T) %>%
                select(tot_win_cnt = win_cnt)) %>%
    mutate(adj_win_cnt = tot_win_cnt * (1 - reasc_rate),
           adj_win_cnt_se = tot_win_cnt * reasc_rate_se) %>%
    bind_cols(bio_df %>%
                group_by(origin) %>%
                summarise(n_tags = n_distinct(tag_code)) %>%
                mutate(prop = n_tags / sum(n_tags),
                       prop_se = sqrt((prop * (1 - prop)) / sum(n_tags)))) %>%
    rowwise() %>%
    mutate(tot_escp = adj_win_cnt * prop,
           tot_escp_se = msm::deltamethod(~ x1 * x2,
                                          mean = c(adj_win_cnt, prop),
                                          cov = diag(c(adj_win_cnt_se, prop_se)^2))) %>%
    select(Species, SpawnYear, origin, reasc_rate, matches('escp'))
  print(org_escape)
  
  #----------------------------------------------------------------
  # better way to use upstream dam counts
  # add origin to prepped capture histories
  pit_obs = prepped_ch %>%
    left_join(bio_df %>%
                select(tag_code,
                       origin)) %>%
    select(tag_code, origin,
           everything())
  
  # get the total dam counts from various dams
  dam_escp_df = tibble(year = yr,
                       dam = as_factor(c("PriestRapids",
                                         "RockIsland")),
                       # "RockyReach"),
                       dam_code = c("PRD",
                                    "RIS"),
                       # "RRH"),
                       pit_code = c("PRA",
                                    "RIA")) %>%
    # "RRF")) %>%
    rowwise() %>%
    mutate(win_cnt = map_dbl(dam_code,
                             .f = function(x) {
                               suppressMessages(STADEM::getWindowCounts(dam = x,
                                                                        spp = "Steelhead",
                                                                        start_date = start_date,
                                                                        end_date = end_date)) %>%
                                 summarise_at(vars(win_cnt),
                                              list(sum),
                                              na.rm = T) %>%
                                 pull(win_cnt)
                             })) %>%
    # add re-ascension data
    mutate(reasc_df = map(pit_code,
                          .f = function(x) {
                            dart_df = try(suppressMessages(queryPITtagData(damPIT = x,
                                                                           spp = "Steelhead",
                                                                           start_date = start_date,
                                                                           end_date = end_date)))
                            if(class(dart_df)[1] == "try-error") {
                              # return(NA)
                              tibble(origin = c("H", "W"),
                                     tot_tags = NA,
                                     reascent_tags = NA) %>%
                                return()
                            } else {
                              dart_df %>%
                                filter(!str_detect(TagId, "000.0")) %>%
                                mutate(SpawnYear = yr) %>%
                                mutate(across(TagIdAscentCount,
                                              as.numeric)) %>%
                                mutate(across(TagIdAscentCount,
                                              tidyr::replace_na,
                                              0)) %>%
                                mutate(ReAscent = ifelse(TagIdAscentCount > 1, T, F)) %>%
                                filter(RearType %in% c("H", "W")) %>%
                                group_by(Species,
                                         origin = RearType,
                                         SpawnYear, Date) %>%
                                summarise(tot_tags = n_distinct(TagId),
                                          reascent_tags = n_distinct(TagId[ReAscent]),
                                          .groups = "drop") %>%
                                group_by(Species,
                                         origin,
                                         SpawnYear) %>%
                                summarise(across(matches('tags'),
                                                 sum,
                                                 na.rm = T),
                                          .groups = "drop") %>%
                                mutate(reasc_rate = reascent_tags / tot_tags,
                                       reasc_rate_se = sqrt(reasc_rate * (1 - reasc_rate) / tot_tags)) %>%
                                select(-Species,
                                       -SpawnYear) %>%
                                return()
                            }
                          })) %>%
    unnest(reasc_df) %>%
    # if no re-ascension data, use Priest Rapids
    arrange(year,
            origin,
            dam) %>%
    fill(reasc_rate,
         reasc_rate_se,
         .direction = "down") %>%
    arrange(year,
            dam,
            origin) %>%
    # adjust for re-ascension
    mutate(adj_win_cnt = win_cnt * (1 - reasc_rate),
           adj_win_cnt_se = win_cnt * reasc_rate_se) %>%
    # crossing(origin = c("W", "H")) %>%
    # break down hatchery and wild based on upstream tags
    mutate(n_obs = map2_int(pit_code,
                            origin,
                            .f = function(x, y) {
                              pit_obs %>%
                                filter(origin == y) %>%
                                summarize(n_obs = n_distinct(tag_code[node == x])) %>%
                                pull(n_obs)
                            })) %>%
    mutate(n_upstrm = map2_int(pit_code,
                               origin,
                               .f = function(x, y) {
                                 pit_obs %>%
                                   filter(origin == y) %>%
                                   summarize(n_path = n_distinct(tag_code[str_detect(path, paste0(" ", x))])) %>%
                                   pull(n_path)
                               }),
           n_upstrm = if_else(dam == "PriestRapids",
                              n_obs,
                              n_upstrm)) %>%
    # generate proportion hatchery/wild based on upstream tags
    group_by(dam) %>%
    mutate(prop_org = n_upstrm / sum(n_upstrm),
           prop_org_se = sqrt((prop_org * (1 - prop_org)) / sum(n_upstrm))) %>%
    # generate proportion hatchery/wild based on tags detected at dam
    # mutate(prop_org = n_obs / sum(n_obs),
    #        prop_org_se = sqrt((prop_org * (1 - prop_org)) / sum(n_obs))) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(tot_escp = win_cnt * (1 - reasc_rate) * prop_org,
           tot_escp_se = msm::deltamethod(~ x1 * (1 - x2) * x3,
                                          mean = c(win_cnt,
                                                   reasc_rate,
                                                   prop_org),
                                          cov = diag(c(0,
                                                       reasc_rate_se,
                                                       prop_org_se)^2))) %>%
    ungroup() %>%
    # add estimated transistion parameters
    left_join(trans_summ %>%
                filter(param %in% c("RIA", "RRF")) %>%
                select(origin,
                       pit_code = param,
                       trans_est = mean,
                       trans_se = sd) %>%
                bind_rows(crossing(origin = c("H", "W"),
                                   pit_code = "PRA",
                                   trans_est = 1,
                                   trans_se = 0))) %>%
    rowwise() %>%
    # calculate equivalent adjusted count at Priest
    mutate(priest_cnt = tot_escp / trans_est,
           priest_cnt_se = msm::deltamethod(~ x1 / x2,
                                            mean = c(tot_escp,
                                                     trans_est),
                                            cov = diag(c(tot_escp_se,
                                                         trans_se)^2))) %>%
    ungroup()
  
  
  for(dam_cnt_name in c("PriestRapids",
                        "RockIsland")) {
    
    cat(paste("\t Using", dam_cnt_name, " dam \n\n"))
    
    org_escape <- dam_escp_df %>%
      filter(dam == dam_cnt_name) %>%
      mutate(Species = "Steelhead") %>%
      select(Species,
             SpawnYear = year,
             dam,
             origin,
             reasc_rate,
             tot_escp = priest_cnt,
             tot_escp_se = priest_cnt_se)
    
    
    # translate movement estimates to escapement
    escape_post = trans_df %>%
      left_join(org_escape %>%
                  group_by(origin) %>%
                  summarise(tot_esc_samp = map2(tot_escp,
                                                tot_escp_se,
                                                .f = function(x, y) {
                                                  tibble(tot_escp = rnorm(max(trans_df$iter),
                                                                          mean = x,
                                                                          sd = y)) %>%
                                                    mutate(iter = 1:n())
                                                })) %>%
                  unnest(cols = tot_esc_samp)) %>%
      mutate(escp = value * tot_escp)
    
    escape_summ = escape_post %>%
      group_by(origin, location = param) %>%
      summarise(mean = mean(escp),
                median = median(escp),
                mode = estMode(escp),
                sd = sd(escp),
                skew = moments::skewness(escp),
                kurtosis = moments::kurtosis(escp),
                lowerCI = coda::HPDinterval(coda::as.mcmc(escp))[,1],
                upperCI = coda::HPDinterval(coda::as.mcmc(escp))[,2],
                .groups = 'drop') %>%
      mutate(across(c(mean, median, mode, sd, matches('CI$')),
                    ~ if_else(. < 0, 0, .))) %>%
      mutate(across(c(mean, median, mode, sd, skew, kurtosis, matches('CI$')),
                    round,
                    digits = 2)) %>%
      arrange(desc(origin), location) %>%
      tibble::add_column(species = "Steelhead",
                         spawn_year = yr,
                         .before = 0)
    
    #-----------------------------------------------------------------
    # save some of these objects
    save(tag_summ,
         bio_df,
         trans_df,
         trans_summ,
         dam_escp_df,
         org_escape,
         escape_post,
         escape_summ,
         detect_summ,
         configuration,
         flowlines,
         parent_child,
         sites_sf,
         file = here("analysis/data/derived_data/estimates",
                     dam_cnt_name,
                     paste0("UC_Sthd_DABOM_", yr, ".rda")))
  }
}
comp_df_trial <- left_join(comp_df, org_escape, by="origin")
comp_df_trial %<>% mutate(adj_win_cnt = win_cnt * (1 - reasc_rate))
comp_df_trial %>%
  ggplot(aes(x = dam.x,
             y = mean,
             color = "DABOM")) +
  geom_errorbar(aes(ymin = lowerCI,
                    ymax = upperCI),
                width = 0) +
  geom_point(size = 3) +
  geom_point(aes(y = win_cnt,
                 color = "Dam Count"),
             size = 3,
             position = position_dodge(width = 1)) +
  geom_point(aes(y = adj_win_cnt,
                 color = "Adj. Dam Count"),
             size = 3,
             position = position_dodge(width = 1)) +
  scale_color_manual(values = c("DABOM" = "gray20",
                                "Dam Count" = "red",
                                "Adj. Dam Count" = "blue"),
                     name = "Source") +
  facet_wrap(~ origin,
             scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom") +
  labs(x = "Dam",
       y = "Estimate",
       title = paste("Steelhead", yr),
       subtitle = "Using Rock Island Counts")
