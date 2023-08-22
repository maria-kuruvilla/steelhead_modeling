# Author: Kevin See
# Purpose: format results to be saved
# Created: 11/30/22
# Last Modified: 4/12/2023
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
# library(UCSthdReddObsErr)
library(PITcleanr)
library(janitor)
library(readxl)
library(writexl)
library(magrittr)
library(msm)
library(here)
library(DescTools)
library(sroem)
devtools::install_github("KevinSee/sroem")


#####
#editing sroem function

correlate_rchs <- function(redd_df = NULL,
                           date_nm = "SurveyDate",
                           cor_redd_nm = "NewRedds",
                           reach_nm = "Reach",
                           use = "pairwise.complete.obs",
                           ...) {
  if (is.null(redd_df)) {
    stop("redd data must be supplied")
  }
  
  n_rchs <- redd_df |>
    pull({{ reach_nm }}) %>%
    n_distinct()
  
  if (n_rchs == 1) {
    rch_nm <- redd_df |>
      pull({{ reach_nm }}) %>%
      unique()
    cor_mat <- matrix(1, 1, 1,
                      dimnames = list(
                        rch_nm,
                        rch_nm
                      )
    )
  } else {
    # check the class of survey date column
    date_class <- redd_df |>
      pull({{ date_nm }}) |>
      class()
    redd_df <- redd_df |>
        mutate(across(
          {{ date_nm }},
          as.POSIXct
        ))
    
    
    cor_mat <- redd_df |>
      rename(
        redds = {{ cor_redd_nm }},
        reach = {{ reach_nm }},
        surv_date = {{ date_nm }}
      ) |>
      mutate(surv_period = lubridate::week(surv_date)) |>
      group_by(reach, surv_period) %>%
      summarize(across(redds,
                       mean,
                       na.rm = T
      )) %>%
      pivot_wider(
        names_from = reach,
        values_from = redds,
        names_sort = T
      ) |>
      select(-surv_period) |>
      # corrr::correlate(...)
      stats::cor(
        use = use,
        ...
      )
  }
  
  return(cor_mat)
}

summarize_redds <- function(redd_df = NULL,
                            species = c("Steelhead", "Spring Chinook"),
                            group_vars = c(river, reach, index, survey_type),
                            summ_vars = c("river", "location"),
                            new_redd_nm = "new_redds",
                            vis_redd_nm = "visible_redds",
                            net_err_nm = "net_error",
                            net_se_nm = "net_error_se",
                            min_non0_wks = 3,
                            min_redds = 2,
                            gauc = NULL,
                            add_zeros = F,
                            use_cor = F,
                            date_nm = "survey_date",
                            cor_redd_nm = "new_redds",
                            reach_nm = "reach",
                            use = "pairwise.complete.obs",
                            method = c(
                              "pearson",
                              "kendall",
                              "spearman"
                            ),
                            ...) {
  if (is.null(redd_df)) {
    stop("redd data must be supplied")
  }
  
  # estimate redds for each reach using net error estimate and GAUC
  rch_est <- estimate_redds(
    redd_df,
    species = species,
    group_vars = unique(c(group_vars, summ_vars)),
    new_redd_nm = new_redd_nm,
    vis_redd_nm = vis_redd_nm,
    net_err_nm = net_err_nm,
    net_se_nm = net_se_nm,
    min_non0_wks = min_non0_wks,
    min_redds = min_redds,
    gauc = gauc,
    add_zeros = add_zeros
  )
  
  # summarize results across summ_vars
  all_rchs <- redd_df %>%
    pull({{ reach_nm }}) %>%
    unique()
  
  summ_est <- rch_est %>%
    dplyr::filter(.data[[reach_nm]] %in% all_rchs) %>%
    dplyr::group_by(dplyr::across({{ summ_vars }})) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      n_rchs = purrr::map_dbl(data,
                              .f = function(x) {
                                x %>%
                                  rename(rch = {{ reach_nm }}) %>%
                                  summarize(n_rchs = n_distinct(rch)) %>%
                                  pull(n_rchs) %>%
                                  return()
                              }
      ),
      strm_obs = purrr::map_dbl(data,
                                .f = function(x) {
                                  sum(x$tot_feat, na.rm = T)
                                }
      ),
      strm_est = purrr::map_dbl(data,
                                .f = function(x) {
                                  sum(x$redd_est, na.rm = T)
                                }
      ),
      strm_se = purrr::map_dbl(data,
                               .f = function(x) {
                                 sqrt(sum(x$redd_se, na.rm = T)^2)
                               }
      )
    )
  
  
  return_list <- list(
    rch_est = rch_est,
    summ_est = summ_est
  )
  
  if (use_cor) {
    # generate correlation matrix
    cor_df <- redd_df %>%
      dplyr::group_by(across({{ summ_vars }})) %>%
      tidyr::nest() %>%
      dplyr::summarize(
        cor_mat = purrr::map(data,
                             .f = function(x) {
                               correlate_rchs(x,
                                              date_nm = date_nm,
                                              cor_redd_nm = cor_redd_nm,
                                              reach_nm = reach_nm,
                                              use = use,
                                              method = method,
                                              ...
                               )
                             }
        ),
        .groups = "drop"
      )
    
    return_list$summ_est <- summ_est %>%
      rename(strm_se_naive = strm_se) %>%
      dplyr::left_join(cor_df,
                       by = {{ summ_vars }}
      ) %>%
      # use correlations between reaches to get appropriate standard error
      dplyr::mutate(strm_se = purrr::map2_dbl(data,
                                              cor_mat,
                                              .f = function(x, y) {
                                                se <- try(
                                                  deltamethod(as.formula(paste("~", paste0("x", 1:nrow(x), collapse = "+"))),
                                                              mean = x$redd_est,
                                                              cov = diag(
                                                                x = x$redd_se,
                                                                nrow = nrow(y)
                                                              ) %*% y %*%
                                                                diag(
                                                                  x = x$redd_se,
                                                                  nrow = nrow(y)
                                                                )
                                                  ),
                                                  silent = T
                                                )
                                                
                                                if (class(se)[1] == "try-error" | is.na(se)) {
                                                  se <- sqrt(sum(x$redd_se^2, na.rm = T))
                                                }
                                                
                                                
                                                return(se)
                                              }
      ))
    
    return_list <- c(
      return_list,
      list(cor_df = cor_df)
    )
  }
  
  return(return_list)
}


query_redd_data <- function(
    redd_file_path = "data",
    redd_file_name = "STHD_Wenatchee_Redd_Surveys.xlsx",
    experience_path = "data",
    experience_file_name = "STHD_Surveyor_Experience.xlsx",
    query_year = yrs) {
  
  data_file = paste(redd_file_path,
                    redd_file_name,
                    sep = "/")
  
  if(!file.exists(data_file)) {
    stop("File not found.")
  }
  
  data_list <- readxl::excel_sheets(data_file) |>
    as.list() |>
    rlang::set_names() |>
    purrr::map(.f = purrr::quietly(function(x) {
      readxl::read_excel(data_file,
                         sheet = x) |>
        janitor::clean_names()
    })) |>
    purrr::map("result")
  
  
  redd_surv_df <- data_list$`Redd Surveys` |>
    dplyr::mutate(
      across(
        c(survey_type,
          river),
        stringr::str_to_title
      )
    )
  
  # to-date, only grab experience covariates from Wenatchee data
  if(stringr::str_detect(redd_file_name, "Wenatchee")) {
    
    # get experience data
    exp_df <- suppressMessages(readxl::read_excel(paste(experience_path,
                                                        experience_file_name,
                                                        sep = "/"),
                                                  sheet = "Experience",
                                                  skip = 1)) |>
      dplyr::rename(basin = `...1`,
                    surveyor_initials = `...2`) |>
      tidyr::pivot_longer(-c(basin:surveyor_initials),
                          names_to = "spawn_year",
                          values_to = "experience") |>
      dplyr::mutate(
        dplyr::across(
          spawn_year,
          as.numeric
        )
      )
    
    redd_surv_df <- redd_surv_df |>
      dplyr::mutate(
        dplyr::across(
          c(surveyor1,
            surveyor2),
          ~ stringr::str_remove(.,
                                "\\ \\([:alpha:]+\\)"))) |>
      dplyr::left_join(exp_df |>
                         dplyr::select(spawn_year,
                                       surveyor1 = surveyor_initials,
                                       exp1 = experience),
                       by = c("spawn_year", "surveyor1")) |>
      dplyr::left_join(exp_df |>
                         dplyr::select(spawn_year,
                                       surveyor2 = surveyor_initials,
                                       exp2 = experience),
                       by = c("spawn_year", "surveyor2")) |>
      dplyr::rowwise() |>
      dplyr::mutate(exp_sp_total = mean(c(exp1, exp2), na.rm = T)) |>
      dplyr::ungroup()
  }
  
  redd_df <- redd_surv_df |>
    dplyr::left_join(data_list$`Reach Length` |>
                       dplyr::group_by(river,
                                       reach,
                                       type, index) |>
                       dplyr::summarize(
                         dplyr::across(length_km,
                                       ~ sum(.))),
                     by = c("river", "reach", "index")) |>
    dplyr::left_join(data_list$`Thalweg CV` |>
                       dplyr::select(river,
                                     reach,
                                     mean_thalweg_cv),
                     by = c("river", "reach")) |>
    dplyr::left_join(data_list$`Discharge Gages` |>
                       dplyr::select(reach,
                                     usgs_site_code = site_code),
                     by = "reach") |>
    # left_join(data_list$Discharge,
    #           by = c("spawn_year", "river", "reach", "index", "survey_type", "survey_date")) |>
    dplyr::mutate(
      dplyr::across(
        c(reach,
          river),
        as.factor),
      dplyr::across(
        reach,
        forcats::fct_relevel,
        "W10",
        "C1", "N1", "P1",
        after = Inf),
      dplyr::across(
        reach,
        forcats::fct_relevel,
        "MH1", "T1", "WN1",
        after = Inf)) |>
    dplyr::filter(spawn_year %in% query_year) |>
    # calculate redd density and log of experience
    dplyr::mutate(naive_density_km = visible_redds / length_km,
                  exp_sp_total_log = log(exp_sp_total))
  
  if(nrow(redd_df) == 0) {
    message(paste("No redd data found for",
                  paste(query_year, collapse = ", "),
                  ".\n"))
    return(NULL)
  }
  
  message("\t Querying USGS for discharge data\n")
  
  # query USGS for mean daily discharge data
  discharge_df <- redd_df |>
    dplyr::filter(!is.na(usgs_site_code)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      dataRetrieval::readNWISdv(usgs_site_code,
                                parameterCd = "00060", # discharge
                                startDate = as.character(lubridate::ymd(survey_date)),
                                endDate = as.character(lubridate::ymd(survey_date)),
                                statCd = "00003" # mean
      )) |>
    dplyr::ungroup() |>
    dplyr::rename(mean_discharge = X_00060_00003) |>
    dplyr::select(-c(agency_cd:Date, X_00060_00003_cd))
  
  # adjust discharge for W10: Plain - Chiwawa discharge
  if("W10" %in% unique(redd_df$reach)) {
    plain_code <- redd_df |>
      dplyr::filter(reach == "W9") |>
      dplyr::pull(usgs_site_code) |>
      unique()
    
    w10_discharge <- discharge_df |>
      dplyr::filter(reach == "W10") |>
      dplyr::rowwise() |>
      dplyr::mutate(
        dataRetrieval::readNWISdv(plain_code,
                                  parameterCd = "00060", # discharge
                                  startDate = as.character(lubridate::ymd(survey_date)),
                                  endDate = as.character(lubridate::ymd(survey_date)),
                                  statCd = "00003" # mean
        )) |>
      dplyr::ungroup() |>
      dplyr::rename(plain_discharge = X_00060_00003) |>
      dplyr::select(-c(agency_cd:Date, X_00060_00003_cd)) |>
      dplyr::mutate(mean_discharge = plain_discharge - mean_discharge) |>
      dplyr::select(-plain_discharge)
    
    # put it all back together
    discharge_df |>
      dplyr::filter(reach != "W10") |>
      dplyr::bind_rows(w10_discharge) -> discharge_df
    
  }
  
  redd_data <- discharge_df |>
    dplyr::bind_rows(redd_df |>
                       dplyr::filter(is.na(usgs_site_code))) |>
    dplyr::arrange(spawn_year,
                   river,
                   reach,
                   index,
                   survey_date)
  
  # identical(nrow(redd_df),
  #           nrow(redd_data))
  
  
  return(redd_data)
  
}

prep_wen_sthd_data <- function(
    redd_file_path = "T:/DFW-Team FP Upper Columbia Escapement - General/UC_Sthd/inputs/Redd Data",
    redd_file_name = "Wenatchee_Redd_Surveys.xlsx",
    experience_path = redd_file_path,
    experience_file_name = redd_file_name,
    dabom_file_path = "T:/DFW-Team FP Upper Columbia Escapement - General/UC_Sthd/estimates",
    dabom_file_name = "UC_STHD_Model_Output.xlsx",
    brood_file_path = "T:/DFW-Team FP Upper Columbia Escapement - General/UC_Sthd/inputs/Bio Data/Sex and Origin PRD-Brood Comparison Data",
    brood_file_name = "STHD UC Brood Collections_2011 to current.xlsx",
    removal_file_path = "T:/DFW-Team FP Upper Columbia Escapement - General/UC_Sthd/inputs/Fish Removals",
    removal_file_name = "Master_STHD_Removals_2.18.23.MH.xlsx",
    n_observers = "two",
    query_year = lubridate::year(lubridate::today()) - 1,
    save_rda = F,
    save_by_year = T,
    save_file_path = here::here("analysis/data/derived_data"),
    save_file_name = NULL
) {
  
  message("\t Gathering redd data.\n")
  
  # load data for selected years
  redd_df_all <- query_redd_data(redd_file_path,
                                 redd_file_name,
                                 experience_path,
                                 experience_file_name,
                                 query_year = query_year)
  
  if(!is.null(redd_df_all)) {
    # divide reaches into various location categories
    redd_df_all <- redd_df_all |>
      dplyr::mutate(location = dplyr::if_else(reach %in% paste0("W", 8:10),
                                              "Above Tumwater",
                                              dplyr::if_else(reach %in% paste0("W", 1:7),
                                                             "Below Tumwater",
                                                             "Tributaries")))
    
    # predict net error
    redd_df_all <- redd_df_all |>
      sroem::predict_neterr(species = "Steelhead",
                            num_obs = n_observers)
  }
  
  #-----------------------------------------------------------------
  # load data on error calls for sex at Priest Rapids when fish were tagged
  
  message("\t Pulling PIT tag data.\n\n")
  
  # get info on tags detected somewhere in the Wenatchee
  wen_tags_all <- readxl::read_excel(paste(dabom_file_path,
                                           dabom_file_name,
                                           sep = "/"),
                                     sheet = "Tag Summary") |>
    janitor::clean_names() |>
    dplyr::filter(str_detect(path, "LWE"),
                  spawn_year %in% query_year) |>
    dplyr::mutate(location = dplyr::if_else(spawn_node %in% c('TUM', 'UWE'),
                                            'Above Tumwater',
                                            dplyr::if_else(str_detect(spawn_node, "^LWE"),
                                                           "Below Tumwater",
                                                           dplyr::if_else(str_detect(path, "CHL"),
                                                                          "Chiwawa",
                                                                          dplyr::if_else(str_detect(path, "NAL"),
                                                                                         "Nason",
                                                                                         dplyr::if_else(str_detect(path, "PES"),
                                                                                                        "Peshastin",
                                                                                                        "Other Tributaries"))))),
                  dplyr::across(location,
                                factor,
                                levels = c("Below Tumwater",
                                           'Above Tumwater',
                                           "Peshastin",
                                           "Nason",
                                           "Chiwawa",
                                           'Other Tributaries'))) |>
    dplyr::select(spawn_year,
                  tag_code,
                  location,
                  origin,
                  sex)
  
  #-------------------------------------------------------
  # generate fish / redd and pHOS for different areas
  fpr_all = wen_tags_all |>
    dplyr::group_by(spawn_year,
                    location) |>
    dplyr::summarize(n_male = n_distinct(tag_code[sex == "M"]),
                     n_female = n_distinct(tag_code[sex == "F"]),
                     n_sexed = n_male + n_female,
                     n_wild = n_distinct(tag_code[origin == "W"]),
                     n_hatch = n_distinct(tag_code[origin == "H"]),
                     n_origin = n_wild + n_hatch,
                     .groups = "drop") |>
    dplyr::mutate(prop_m = n_male / n_sexed,
                  prop_se = sqrt((prop_m * (1 - prop_m)) / (n_sexed)),
                  fpr = (prop_m) / (1 - prop_m) + 1) |>
    dplyr::rowwise() |>
    dplyr::mutate(fpr_se = msm::deltamethod(~ x1 / (1 - x1) + 1,
                                            mean = prop_m,
                                            cov = prop_se^2)) |>
    dplyr::ungroup() |>
    dplyr::mutate(phos = n_hatch / n_origin,
                  phos_se = sqrt((phos * (1 - phos)) / (n_origin)))
  
  message("\t Adjusting fish/redd.\n")
  
  # adjust fish / redd for errors in Priest sex calls
  # the excel file contains rounded numbers, so re-calculate
  # various statistics for use in analyses
  # estimate error rate for each sex
  sex_err_rate <- readxl::read_excel(paste(dabom_file_path,
                                           dabom_file_name,
                                           sep = "/"),
                                     sheet = "Tag Summary") |>
    janitor::clean_names() |>
    select(spawn_year,
           tag_code,
           sex_field = sex) |>
    inner_join(read_excel(paste(brood_file_path,
                                brood_file_name,
                                sep = "/"),
                          sheet = "Brood Collected_PIT Tagged Only") |>
                 clean_names() |>
                 rename(tag_code = recaptured_pit) |>
                 select(spawn_year,
                        tag_code,
                        sex_final) |>
                 distinct(),
               by = c("spawn_year",
                      "tag_code")) |>
    filter(!is.na(sex_final),
           !is.na(sex_field)) |>
    mutate(agree = if_else(sex_field == sex_final,
                           T, F)) |>
    group_by(spawn_year,
             sex = sex_field) |>
    summarize(n_tags = n_distinct(tag_code),
              n_true = sum(agree),
              n_false = sum(!agree),
              .groups = "drop") |>
    mutate(binom_ci = map2(n_false,
                           n_tags,
                           .f = function(x, y) {
                             DescTools::BinomCI(x, y) |>
                               as_tibble()
                           })) |>
    unnest(binom_ci) |>
    clean_names() |>
    rename(perc_false = est,
           lowerci = lwr_ci,
           upperci = upr_ci) |>
    mutate(perc_se = sqrt((perc_false * (1 - perc_false)) / n_tags)) |>
    relocate(perc_se,
             .after = "perc_false")
  
  #
  #
  # sex_err_rate <- readxl::read_excel(paste(dabom_file_path,
  #                                          dabom_file_name,
  #                                          sep = "/"),
  #                                    sheet = "Sex Error Rates") |>
  #   janitor::clean_names() |>
  #   dplyr::select(spawn_year:n_false) |>
  #   dplyr::filter(spawn_year %in% query_year) |>
  #   dplyr::rowwise() |>
  #   dplyr::mutate(binom_ci = map2(n_false,
  #                                 n_tags,
  #                                 .f = function(x, y) {
  #                                   DescTools::BinomCI(x, y) |>
  #                                     as_tibble()
  #                                 })) |>
  #   tidyr::unnest(binom_ci) |>
  #   janitor::clean_names() |>
  #   dplyr::rename(perc_false = est) |>
  #   dplyr::mutate(perc_se = sqrt((perc_false * (1 - perc_false)) / n_tags))
  
  adj_fpr <- fpr_all |>
    dplyr::select(spawn_year,
                  location,
                  n_male,
                  n_female) |>
    tidyr::pivot_longer(cols = c(n_male,
                                 n_female),
                        names_to = "sex",
                        values_to = "n_fish") |>
    dplyr::mutate(
      dplyr::across(sex,
                    str_remove,
                    "^n_"),
      dplyr::across(sex,
                    str_to_title)) |>
    dplyr::mutate(
      dplyr::across(sex,
                    recode,
                    "Male" = "M",
                    "Female" = "F")) |>
    dplyr::left_join(sex_err_rate |>
                       dplyr::select(spawn_year,
                                     sex,
                                     dplyr::starts_with("perc_")),
                     by = c("spawn_year", "sex")) |>
    tidyr::pivot_wider(names_from = sex,
                       values_from = c(n_fish,
                                       perc_false,
                                       perc_se)) |>
    dplyr::mutate(true_male = n_fish_M - (n_fish_M * perc_false_M) + (n_fish_F * perc_false_F),
                  true_female = n_fish_F - (n_fish_F * perc_false_F) + (n_fish_M * perc_false_M),
                  dplyr::across(starts_with("true"),
                                janitor::round_half_up)) |>
    dplyr::rowwise() |>
    dplyr::mutate(true_m_se = msm::deltamethod(~ x1 - (x1 * x2) + (x3 * x4),
                                               mean = c(n_fish_M,
                                                        perc_false_M,
                                                        n_fish_F,
                                                        perc_false_F),
                                               cov = diag(c(0,
                                                            perc_se_M,
                                                            0,
                                                            perc_se_F)^2)),
                  true_f_se = msm::deltamethod(~ x1 - (x1 * x2) + (x3 * x4),
                                               mean = c(n_fish_F,
                                                        perc_false_F,
                                                        n_fish_M,
                                                        perc_false_M),
                                               cov = diag(c(0,
                                                            perc_se_F,
                                                            0,
                                                            perc_se_M)^2))) |>
    dplyr::mutate(n_sexed = true_male + true_female,
                  prop_m = true_male / (true_male + true_female),
                  prop_se = msm::deltamethod(~ x1 / (x1 + x2),
                                             mean = c(true_male,
                                                      true_female),
                                             cov = diag(c(true_m_se,
                                                          true_f_se)^2)),
                  fpr = (prop_m) / (1 - prop_m) + 1,
                  fpr_se = msm::deltamethod(~ x1 / (1 - x1) + 1,
                                            mean = prop_m,
                                            cov = prop_se^2)) |>
    dplyr::ungroup() |>
    dplyr::rename(n_male = true_male,
                  n_female = true_female) |>
    dplyr::left_join(fpr_all |>
                       dplyr::select(spawn_year,
                                     location,
                                     n_wild,
                                     n_hatch,
                                     n_origin,
                                     starts_with("phos")),
                     by = c("spawn_year", "location")) |>
    dplyr::select(dplyr::any_of(names(fpr_all)))
  
  # # look at changes to fish/redd
  # fpr_all |>
  #   select(spawn_year,
  #          location,
  #          old_fpr = fpr) |>
  #   left_join(adj_fpr |>
  #               select(spawn_year,
  #                      location,
  #                      adj_fpr = fpr))
  
  # if any fpr values are Inf, use the older ones
  if(sum(adj_fpr$fpr == Inf) > 0) {
    adj_fpr <- adj_fpr |>
      dplyr::left_join(fpr_all |>
                         dplyr::select(location,
                                       old_fpr = fpr,
                                       old_se = fpr_se)) |>
      dplyr::mutate(fpr = dplyr::if_else(is.na(fpr) | fpr == Inf,
                                         old_fpr,
                                         fpr),
                    fpr_se = dplyr::if_else(is.na(fpr_se) | fpr_se == Inf,
                                            old_se,
                                            fpr_se)) |>
      dplyr::select(-dplyr::starts_with("old"))
  }
  
  fpr_all <- adj_fpr
  
  rm(adj_fpr)
  
  #-----------------------------------------------------------------
  # read in data about known removals of fish prior to spawning
  if(file.exists(paste(removal_file_path,
                       removal_file_name,
                       sep = "/"))) {
    if(str_detect(removal_file_name, "csv$")) {
      removal_df <- readr::read_csv(paste(removal_file_path,
                                          removal_file_name,
                                          sep = "/")) |>
        janitor::clean_names() |>
        dplyr::filter(subbasin == "Wenatchee",
                      spawn_year %in% query_year)
    }
    if(str_detect(removal_file_name, "xls$") |
       str_detect(removal_file_name, "xlsx$")) {
      
      removal_df <- readxl::read_excel(paste(removal_file_path,
                                             removal_file_name,
                                             sep = "/"),
                                       skip = 3,
                                       col_names = c("run cycle",
                                                     "spawn year",
                                                     "population",
                                                     "removal location",
                                                     "agency",
                                                     "adult trapping surplus H",
                                                     "brood collections H",
                                                     "harvest H",
                                                     "adult trapping surplus W",
                                                     "brood collections W",
                                                     "harvest W",
                                                     "adult trapping surplus T",
                                                     "brood collections T",
                                                     "harvest T")) |>
        janitor::clean_names() |>
        dplyr::mutate(
          across(
            c(ends_with("_h"),
              ends_with("_w"),
              ends_with("_t")),
            as.numeric
          )) |>
        dplyr::filter(!is.na(spawn_year)) |>
        dplyr::select(spawn_year:harvest_w) |>
        tidyr::pivot_longer(cols = c(ends_with("_h"),
                                     ends_with("_w")),
                            names_to = "source",
                            values_to = "removed") |>
        dplyr::mutate(origin = str_sub(source, -1)) |>
        dplyr::relocate(origin,
                        .before = "removed") |>
        dplyr::filter(origin %in% c("h", "w")) |>
        dplyr::mutate(
          across(
            source,
            str_remove,
            "_h$"),
          across(
            source,
            str_remove,
            "_w$"),
          across(
            source,
            ~ str_to_title(str_replace_all(., "_", " "))
          ),
          across(
            origin,
            recode,
            "h" = "Hatchery",
            "w" = "Natural"
          )
        ) |>
        dplyr::filter(spawn_year %in% query_year,
                      population == "Wenatchee")
    }
    
  } else {
    message("Removal data not found.\n")
    removal_df <- NULL
  }
  
  
  #-----------------------------------------------------------------
  # pull in some estimates from DABOM
  
  message("\t Gathering PIT escapement estimates.\n")
  
  all_escp = readxl::read_excel(paste(dabom_file_path,
                                      dabom_file_name,
                                      sep = "/"),
                                sheet = "Run Escp All Locations") |>
    janitor::clean_names() |>
    dplyr::filter(spawn_year %in% query_year,
                  location %in% c('ICL',
                                  'PES',
                                  'MCL',
                                  'CHM',
                                  'CHW',
                                  'CHL',
                                  'NAL',
                                  'LWN',
                                  'WTL',
                                  'LWE',
                                  'LWE_bb',
                                  'TUM_bb',
                                  'UWE_bb'))
  
  # pull out estimates of tributary spawners from DABOM
  trib_spawners_all = all_escp |>
    dplyr::filter(location %in% c('ICL',
                                  'PES',
                                  'MCL',
                                  'CHM',
                                  'CHW',
                                  'CHL',
                                  'NAL',
                                  'LWN',
                                  'WTL')) |>
    dplyr::select(spawn_year,
                  origin,
                  location,
                  spawners = estimate,
                  spawners_se = se) |>
    dplyr::mutate(
      dplyr::across(origin,
                    recode,
                    "W" = "Natural",
                    "H" = "Hatchery"),
      dplyr::across(location,
                    recode,
                    'CHL' = 'Chiwawa',
                    'CHM' = 'Chumstick',
                    'CHW' = 'Chiwaukum',
                    'ICL' = 'Icicle',
                    'LWN' = 'Little Wenatchee',
                    'MCL' = 'Mission',
                    'NAL' = 'Nason',
                    'PES' = 'Peshastin',
                    'WTL' = 'White River')) |>
    dplyr::arrange(location, origin)
  
  # pull out mainstem escapement estimates
  escp_wen_all = all_escp |>
    dplyr::filter(location %in% c('LWE',
                                  'LWE_bb',
                                  'TUM_bb',
                                  'UWE_bb')) |>
    dplyr::mutate(
      dplyr::across(location,
                    recode,
                    'LWE' = 'Wen_all',
                    'LWE_bb' = 'Below Tumwater',
                    'TUM_bb' = "Above Tumwater",
                    'UWE_bb' = 'Above Tumwater')) |>
    dplyr::mutate(
      dplyr::across(origin,
                    recode,
                    "W" = "Natural",
                    "H" = "Hatchery")) |>
    dplyr::group_by(spawn_year,
                    location,
                    origin) |>
    dplyr::summarise(
      dplyr::across(estimate,
                    sum),
      dplyr::across(se,
                    ~ sqrt(sum(.^2))),
      .groups = "drop")
  
  #-----------------------------------------------------------------
  # save
  if(save_rda & save_by_year) {
    for(yr in query_year) {
      message(paste("Saving data from spawn year",
                    yr,
                    ".\n\n"))
      
      if(!is.null(redd_df_all)) {
        redd_df <- redd_df_all |>
          dplyr::filter(spawn_year == yr)
      } else {
        redd_df <- NULL
      }
      
      wen_tags <- wen_tags_all |>
        dplyr::filter(spawn_year == yr)
      
      sex_err <- sex_err_rate |>
        dplyr::filter(spawn_year == yr)
      
      fpr_df <- fpr_all |>
        dplyr::filter(spawn_year == yr)
      
      trib_spawners <- trib_spawners_all |>
        dplyr::filter(spawn_year == yr)
      
      escp_wen <- escp_wen_all |>
        dplyr::filter(spawn_year == yr)
      
      rem_df <- removal_df |>
        dplyr::filter(spawn_year == yr)
      
      if(is.null(save_file_name)) {
        file_nm = paste0('wen_', yr, '.rda')
      } else {
        file_nm = save_file_name
      }
      
      save(redd_df,
           wen_tags,
           sex_err,
           fpr_df,
           trib_spawners,
           escp_wen,
           rem_df,
           file = paste(save_file_path,
                        file_nm,
                        sep = "/"))
      rm(file_nm)
      
    }
  }
  else {
    if(!is.null(redd_df_all)) {
      redd_df <- redd_df_all
    } else {
      redd_df <- NULL
    }
    
    wen_tags <- wen_tags_all
    
    sex_err <- sex_err_rate
    
    fpr_df <- fpr_all
    
    trib_spawners <- trib_spawners_all
    
    escp_wen <- escp_wen_all
    
    rem_df <- removal_df
    
    
    if(save_rda & !save_by_year) {
      if(is.null(save_file_name)) {
        if(length(query_year) > 1) {
          save_file_name <- paste0('wen_',
                                   paste(min(query_year),
                                         max(query_year),
                                         sep = "-"),
                                   '.rda')
        } else {
          save_file_name <- paste0('wen_',
                                   query_year,
                                   '.rda')
        }
      }
      save(redd_df,
           wen_tags,
           sex_err,
           fpr_df,
           trib_spawners,
           escp_wen,
           rem_df,
           file = paste(save_file_path,
                        save_file_name,
                        sep = "/"))
    } else {
      
      tmp_file <- tempfile(fileext = ".rda")
      
      save(redd_df,
           wen_tags,
           sex_err,
           fpr_df,
           trib_spawners,
           escp_wen,
           rem_df,
           file = tmp_file)
      
      load(tmp_file,
           envir = .GlobalEnv)
      
      file.remove(tmp_file)
      rm(tmp_file)
    }
  }
  
}

######


# load age table
data("age_table")

#-----------------------------------------------------------------
# function to fix all table names
makeTableNms = function(df) {
  df |>
    rlang::set_names(nm = function(x) {
      x |>
        str_replace_all("_", " ") |>
        str_to_title() |>
        str_replace("Nos", "NOS") |>
        str_replace("Hos", "HOS") |>
        str_replace("Se$", "SE") |>
        str_replace("Sd$", "SD") |>
        str_replace("Cv$", "CV") |>
        str_replace("Fpr", "FpR") |>
        str_replace("Phos", "pHOS") |>
        str_replace("Lowerci", "LCI") |>
        str_replace("Upperci", "UCI") |>
        str_replace("Cwt$", "CWT")
    }) %>%
    return()
}

#-----------------------------------------------------------------
# set the highest year to be included
max_yr = 2022
# max_yr = lubridate::year(lubridate::today()) - 1

#-----------------------------------------------------------------
# generate tables of redds and spawners

# set some thresholds
# minimum number of total redds observed
min_redds = 2
# minimum number of weeks with at least one new redd observed
min_non0_wks = 3


#-----------------------------------------------------------------
spwn_est <- crossing(population = c("Wenatchee"),
                                    #"Methow"),
                     spawn_year = 2014:max_yr) |>
  mutate(run_year = spawn_year - 1) %>%
  select(run_year, spawn_year, population) %>%
  filter(!(population == "Methow" &
             spawn_year < 2021)) %>%
  mutate(results_list = map2(spawn_year,
                             population,
                             .f = possibly(function(yr, pop) {

                               message(paste("Prepping data from", pop, yr, "\n\n"))

                               rm(redd_df,
                                  wen_tags,
                                  met_tags,
                                  sex_err,
                                  fpr_df,
                                  trib_spawners,
                                  escp_wen,
                                  escp_met,
                                  escp_est,
                                  rem_df)

                               if(pop == "Wenatchee") {
                                 prep_wen_sthd_data(
                                   redd_file_path =
                                     here("data"),
                                   redd_file_name = "STHD_Wenatchee_Redd_Surveys.xlsx",
                                   experience_path = here("data"),
                                   experience_file_name = "STHD_Surveyor_Experience.xlsx",
                                   dabom_file_path =
                                     here("data"),
                                   dabom_file_name = "UC_STHD_Model_Output.xlsx",
                                   brood_file_path =
                                     here("data"),
                                   brood_file_name = "STHD_UC Brood Collections_2011 to current.xlsx",
                                   removal_file_path =
                                     here("data"),
                                   removal_file_name = "STHD_Removals.xlsx",
                                   n_observers = "two",
                                   query_year = yr,
                                   save_rda = F,
                                   save_by_year = T,
                                   save_file_path = here::here("output"),
                                   save_file_name = NULL
                                 )
                                 escp_est = escp_wen
                                 rm(escp_wen)
                               } else if(pop == "Methow") {
                                 prep_met_sthd_data(query_year = yr,
                                                    n_observers = "two",
                                                    save_rda = F)
                                 escp_est = escp_met
                                 rm(escp_met)
                               } else {
                                 return(NULL)
                               }
                               # for 2020, we'll need to do something different
                               if(yr != 2020) {

                                 results_lst <- summarize_redds(redd_df,
                                                                species = "Steelhead",
                                                                group_vars = c("river", "reach", "index", "survey_type"),
                                                                summ_vars = c("river", "location", "index"),
                                                                min_non0_wks = min_non0_wks,
                                                                min_redds = min_redds,
                                                                gauc = T,
                                                                add_zeros = T,
                                                                use_cor = T)


                                 redd_results <- results_lst$rch_est |>
                                   arrange(river,
                                           reach,
                                           index) %>%
                                   relocate(river, reach, location,
                                            .before = 1) |>
                                   mutate(location = if_else(location == "Tributaries",
                                                             as.character(river),
                                                             location))

                                 spwn_rch = redd_results %>%
                                   select(river, location, reach, index, redd_est, redd_se) %>%
                                   left_join(fpr_df %>%
                                               select(location, starts_with('fpr')),
                                             by = "location") %>%
                                   left_join(fpr_df %>%
                                               select(location, starts_with('phos')),
                                             by = "location") %>%
                                   rowwise() %>%
                                   mutate(Tot_Spawners = redd_est * fpr,
                                          Hatchery = redd_est * fpr * phos,
                                          Natural = redd_est * fpr * (1 - phos),
                                          Hatchery_SE = ifelse(Tot_Spawners > 0,
                                                               deltamethod(~ x1 * x2 * x3,
                                                                           mean = c(redd_est, fpr, phos),
                                                                           cov = diag(c(redd_se, fpr_se, phos_se)^2)),
                                                               NA),
                                          Natural_SE = ifelse(Tot_Spawners > 0,
                                                              deltamethod(~ x1 * x2 * (1 - x3),
                                                                          mean = c(redd_est, fpr, phos),
                                                                          cov = diag(c(redd_se, fpr_se, phos_se)^2)),
                                                              NA)) %>%
                                   ungroup() %>%
                                   mutate(reach = factor(reach,
                                                         levels = c(paste0('W', 1:10),
                                                                    "C1",
                                                                    "N1",
                                                                    "P1",
                                                                    paste0("MRW", 1:8)))) %>%
                                   select(river, location, reach,
                                          type = index,
                                          Hatchery:Natural_SE) %>%
                                   mutate(across(type,
                                                 recode,
                                                 'Y' = 'Index',
                                                 'N' = 'Non-Index'),
                                          across(location,
                                                 fct_recode,
                                                 "Below Tumwater (mainstem)" = "Below Tumwater",
                                                 "Above Tumwater (mainstem)" = "Above Tumwater")) %>%
                                   bind_rows(trib_spawners %>%
                                               rename(spawners_SE = spawners_se) %>%
                                               pivot_wider(names_from = origin,
                                                           values_from = c(spawners,
                                                                           spawners_SE),
                                                           names_glue = "{origin}_{.value}") %>%
                                               rlang::set_names(function(x) str_remove(x, "_spawners")) %>%
                                               mutate(river = location) %>%
                                               mutate(type = 'DABOM') %>%
                                               select(-spawn_year)) %>%
                                   mutate(location = factor(location, levels = c('Below Tumwater (mainstem)',
                                                                                 'Above Tumwater (mainstem)',
                                                                                 'Icicle',
                                                                                 'Peshastin',
                                                                                 'Mission',
                                                                                 'Chumstick',
                                                                                 'Chiwaukum',
                                                                                 'Chiwawa',
                                                                                 'Nason',
                                                                                 'Little Wenatchee',
                                                                                 'White River',
                                                                                 "Lower Methow",
                                                                                 "Gold",
                                                                                 "Libby",
                                                                                 "Methow Fish Hatchery",
                                                                                 "Upper Methow",
                                                                                 "Twisp",
                                                                                 "Chewuch",
                                                                                 "Spring Creek",
                                                                                 "Beaver"))) %>%
                                   arrange(location, reach, type) %>%
                                   mutate(h_cv = Hatchery_SE / Hatchery,
                                          n_cv = Natural_SE / Natural)


                                 redd_spwn_yr <- redd_results |>
                                   mutate(reach_type = recode(index,
                                                              'Y' = 'Index',
                                                              'N' = 'Non-Index')) |>
                                   left_join(fpr_df %>%
                                               select(location,
                                                      starts_with("fpr")),
                                             by = "location") |>
                                   select(river,
                                          reach,
                                          reach_type,
                                          observed_redds = tot_feat,
                                          redd_est,
                                          redd_se,
                                          fpr,
                                          fpr_se) %>%
                                   mutate(redd_cv = redd_se / redd_est) %>%
                                   relocate(redd_cv,
                                            .before = fpr) %>%
                                   left_join(spwn_rch %>%
                                               select(river,
                                                      reach,
                                                      reach_type = type,
                                                      NOS_est = Natural,
                                                      NOS_se = Natural_SE,
                                                      NOS_cv = n_cv,
                                                      HOS_est = Hatchery,
                                                      HOS_se = Hatchery_SE,
                                                      HOS_cv = h_cv),
                                             by = c("river", "reach", "reach_type"))

                                 spwn_yr <- spwn_rch |>
                                   group_by(river) |>
                                   summarize(across(c(nos_est = Natural,
                                                      hos_est = Hatchery),
                                                    sum,
                                                    na.rm = T),
                                             across(c(nos_se = Natural_SE,
                                                      hos_se = Hatchery_SE),
                                                    ~ sqrt(sum(.^2, na.rm = T))),
                                             .groups = "drop") %>%
                                   mutate(nos_cv = nos_se / nos_est,
                                          hos_cv = hos_se / hos_est) %>%
                                   mutate(total_spawners = nos_est + hos_est) |>
                                   mutate(across(river,
                                                 recode,
                                                 "Wenatchee" = "Wenatchee (mainstem)",
                                                 "Methow" = "Methow (mainstem)"),
                                          across(river,
                                                 as.factor),
                                          across(river,
                                                 fct_relevel,
                                                 "Wenatchee (mainstem)",
                                                 "Methow (mainstem)",
                                                 after = Inf)) |>
                                   arrange(river) |>
                                   select(river,
                                          total_spawners,
                                          starts_with("nos"),
                                          starts_with("hos"))

                                 spwn_yr %<>%
                                   bind_rows(spwn_yr %>%
                                               summarize(river = "Total",
                                                         across(c(total_spawners,
                                                                  nos_est,
                                                                  hos_est),
                                                                sum),
                                                         across(ends_with("_se"),
                                                                ~ sqrt(sum(.^2, na.rm = T))),
                                                         nos_cv = nos_se / nos_est,
                                                         hos_cv = hos_se / hos_est,
                                                         .groups = "drop")) %>%
                                   rowwise() %>%
                                   mutate(phos = hos_est / (hos_est + nos_est),
                                          phos_se = deltamethod(~ x1 / (x1 + x2),
                                                                mean = c(hos_est,
                                                                         nos_est),
                                                                cov = diag(c(hos_se,
                                                                             nos_se)^2))) %>%
                                   ungroup() %>%
                                   mutate(across(river,
                                                 as.factor),
                                          across(river,
                                                 fct_relevel,
                                                 c("Wenatchee (mainstem)",
                                                   "Methow (mainstem)",
                                                   "Total"),
                                                 after = Inf)) %>%
                                   arrange(river)


                                 list(results_lst = results_lst,
                                      fpr = fpr_df,
                                      redd_results = redd_results,
                                      spwn_rch = spwn_rch,
                                      redd_spwn_yr = redd_spwn_yr,
                                      spwn_yr = spwn_yr) %>%
                                   return()

                               } else {
                                 # data from WDFW (Nate Fuchs' radio telemetry study)
                                 rt_df = tibble(subbasin = "Wenatchee",
                                                year = rep(2015:2016, each = 2),
                                                origin = rep(c("Hatchery", "Natural"), 2),
                                                ow_fish = c(20, 25, 4, 12),
                                                surv_fish = c(16, 24, 3, 9)) |>
                                   bind_rows(tibble(subbasin = "Methow",
                                                    year = rep(2015:2016, each = 2),
                                                    origin = rep(c("Hatchery", "Natural"), 2),
                                                    ow_fish = c(26, 16, 64, 13),
                                                    surv_fish = c(21, 16, 60, 12))) |>
                                   mutate(phi = surv_fish / ow_fish,
                                          phi_se = sqrt((phi * (1 - phi))/ ow_fish))

                                 # add years together
                                 rt_df %<>%
                                   bind_rows(rt_df %>%
                                               group_by(subbasin, origin) %>%
                                               summarize(across(c(ow_fish, surv_fish),
                                                                sum),
                                                         .groups = "drop") %>%
                                               mutate(phi = surv_fish / ow_fish,
                                                      phi_se = sqrt((phi * (1 - phi)) / ow_fish))) %>%
                                   mutate(across(year,
                                                 as.factor)) %>%
                                   mutate(across(year,
                                                 fct_explicit_na,
                                                 na_level = "Total")) %>%
                                   arrange(subbasin,
                                           year,
                                           origin)

                                 spwn_yr <- escp_est %>%
                                   filter(str_detect(location, "_all$")) |>
                                   mutate(subbasin = if_else(str_detect(location,
                                                                        "Met"),
                                                             "Methow",
                                                             if_else(str_detect(location,
                                                                                "Wen"),
                                                                     "Wenatchee",
                                                                     NA_character_))) %>%
                                   relocate(subbasin, .before = 1) %>%
                                   left_join(rem_df %>%
                                               group_by(origin) %>%
                                               summarize(across(removed,
                                                                sum),
                                                         .groups = "drop"),
                                             by = "origin") %>%
                                   mutate(across(removed,
                                                 replace_na,
                                                 0)) |>
                                   rowwise() %>%
                                   mutate(escp = max(0, estimate - removed)) %>%
                                   ungroup() %>%
                                   left_join(rt_df %>%
                                               filter(year == "Total") %>%
                                               select(-year,
                                                      -ends_with("fish")),
                                             by = c("subbasin", "origin")) %>%
                                   rowwise() %>%
                                   mutate(all_spwn = escp * phi) %>%
                                   mutate(all_spwn_se = msm::deltamethod(~ x1 * x2,
                                                                         mean = c(escp, phi),
                                                                         cov = diag(c(se, phi_se)^2))) %>%
                                   ungroup() %>%
                                   left_join(trib_spawners %>%
                                               group_by(origin) %>%
                                               summarize(across(spawners,
                                                                sum),
                                                         across(spawners_se,
                                                                ~ sqrt(sum(.^2))),
                                                         .groups = "drop") %>%
                                               rename(trib_spwn = spawners,
                                                      trib_se = spawners_se),
                                             by = "origin") %>%
                                   mutate(main_spwn = all_spwn - trib_spwn,
                                          main_spwn_se = sqrt(all_spwn_se^2 + trib_se^2),
                                          across(main_spwn,
                                                 ~ if_else(. < 0, 0, .))) %>%
                                   mutate(location = paste(subbasin, "(mainstem)"),
                                          river = location,
                                          type = "RT") %>%
                                   select(river,
                                          origin,
                                          location,
                                          type,
                                          spawners = main_spwn,
                                          spawners_se = main_spwn_se) %>%
                                   bind_rows(trib_spawners %>%
                                               mutate(river = location) |>
                                               mutate(type = 'DABOM')) |>
                                   pivot_wider(names_from = origin,
                                               values_from = c(spawners,
                                                               spawners_se),
                                               names_glue = "{origin}_{.value}") |>
                                   rlang::set_names(function(x) str_remove(x, "_spawners")) |>
                                   mutate(across(
                                     location,
                                     factor,
                                     levels = c('Below_TUM',
                                                'TUM_bb',
                                                'Icicle',
                                                'Peshastin',
                                                'Mission',
                                                'Chumstick',
                                                'Chiwaukum',
                                                'Chiwawa',
                                                'Nason',
                                                'Little Wenatchee',
                                                'White River',
                                                "Wenatchee (mainstem)",
                                                "Lower Methow",
                                                "Gold",
                                                "Libby",
                                                "Methow Fish Hatchery",
                                                "Upper Methow",
                                                "Twisp",
                                                "Chewuch",
                                                "Spring Creek",
                                                "Beaver"))) |>
                                   arrange(location, type) |>
                                   mutate(h_cv = Hatchery_se / Hatchery,
                                          n_cv = Natural_se / Natural,
                                          total_spawners = Natural + Hatchery) |>
                                   select(river,
                                          # reach,
                                          # reach_type = type,
                                          total_spawners,
                                          nos_est = Natural,
                                          nos_se = Natural_se,
                                          nos_cv = n_cv,
                                          hos_est = Hatchery,
                                          hos_se = Hatchery_se,
                                          hos_cv = h_cv)

                                 spwn_yr %<>%
                                   bind_rows(spwn_yr %>%
                                               summarize(river = "Total",
                                                         across(c(total_spawners,
                                                                  nos_est,
                                                                  hos_est),
                                                                sum),
                                                         across(ends_with("_se"),
                                                                ~ sqrt(sum(.^2, na.rm = T))),
                                                         nos_cv = nos_se / nos_est,
                                                         hos_cv = hos_se / hos_est,
                                                         .groups = "drop")) %>%
                                   rowwise() %>%
                                   mutate(phos = hos_est / (hos_est + nos_est),
                                          phos_se = deltamethod(~ x1 / (x1 + x2),
                                                                mean = c(hos_est,
                                                                         nos_est),
                                                                cov = diag(c(hos_se,
                                                                             nos_se)^2))) %>%
                                   ungroup() %>%
                                   mutate(across(river,
                                                 as.factor),
                                          across(river,
                                                 fct_relevel,
                                                 c("Wenatchee (mainstem)",
                                                   "Methow (mainstem)",
                                                   "Total"),
                                                 after = Inf)) %>%
                                   arrange(river)

                                 list(results_lst = NA,
                                      fpr = fpr_df,
                                      redd_results = NA,
                                      spwn_rch = NA,
                                      redd_spwn_yr = NA,
                                      spwn_yr = spwn_yr) %>%
                                   return()
                               }
                             },
                             otherwise = NULL)))

redd_spwn <- spwn_est %>%
  mutate(redd_spwn_yr = map(results_list,
                            function(x) x[["redd_spwn_yr"]])) %>%
  select(-results_list) %>%
  unnest(redd_spwn_yr) %>%
  filter(!is.na(observed_redds)) %>%
  makeTableNms()

spwn_df <- spwn_est %>%
  mutate(spwn_yr = map(results_list,
                       function(x) x[["spwn_yr"]])) %>%
  select(-results_list) %>%
  unnest(spwn_yr) %>%
  makeTableNms()

#-----------------------------------------------------------------
# generate tables based on DABOM output
dabom_est <- crossing(spawn_year = c(2011:max_yr)) %>%
  mutate(run_year = spawn_year - 1) %>%
  select(run_year, spawn_year) %>%
  mutate(dam_cnt_name = if_else(spawn_year %in% c(2011:2015, 2018),
                                "PriestRapids",
                                "RockIsland")) %>%
  mutate(all_results = map2(spawn_year,
                            dam_cnt_name,
                            .f = function(yr, dam_cnt_name) {

                              message(paste("\n\t Estimates from", yr, "\n"))

                              # load data
                              load(here("analysis/data/derived_data/estimates",
                                        dam_cnt_name,
                                        paste0("UC_Sthd_DABOM_", yr, ".rda")))

                              # age proportions that year
                              age_prop_yr <- tag_summ |>
                                filter(!group %in% c("Start", "Other")) |>
                                mutate(across(group,
                                              recode,
                                              "BelowPriest" = "Below Priest",
                                              "WellsPool" = "Wells Pool"),
                                       across(origin,
                                              recode,
                                              "H" = "Hatchery",
                                              "W" = "Natural")) |>
                                mutate(run_year = year - 1) |>
                                select(tag_code,
                                       population = group,
                                       origin,
                                       age) |>
                                left_join(age_table |>
                                            select(age,
                                                   fresh_water,
                                                   salt_water,
                                                   total_age),
                                          by = "age") |>
                                pivot_longer(cols = fresh_water:total_age,
                                             names_to = "age_type",
                                             values_to = "age_label") |>
                                group_by(population,
                                         origin,
                                         age_type) |>
                                count(age_label) |>
                                filter(!is.na(age_label)) |>
                                mutate(sample_size = sum(n),
                                       prop = n / sample_size) |>
                                ungroup() |>
                                # add data for all Priest tags
                                bind_rows(tag_summ |>
                                            mutate(across(origin,
                                                          recode,
                                                          "H" = "Hatchery",
                                                          "W" = "Natural")) |>
                                            mutate(run_year = year - 1) |>
                                            select(tag_code,
                                                   population = group,
                                                   origin,
                                                   age) |>
                                            left_join(age_table |>
                                                        select(age,
                                                               fresh_water,
                                                               salt_water,
                                                               total_age),
                                                      by = "age") |>
                                            pivot_longer(cols = fresh_water:total_age,
                                                         names_to = "age_type",
                                                         values_to = "age_label") |>
                                            group_by(origin,
                                                     age_type) |>
                                            count(age_label) |>
                                            filter(!is.na(age_label)) |>
                                            mutate(sample_size = sum(n),
                                                   prop = n / sample_size) |>
                                            ungroup() |>
                                            add_column(population = "Priest Rapids Tags",
                                                       .before = 0))

                              # population escapement that year
                              pop_escp_yr <- escape_post %>%
                                filter(param %in% c('LWE',
                                                    'ENL',
                                                    'LMR',
                                                    'OKL', 'FST',
                                                    'ICH', 'JD1', 'JDA', 'PRH', 'PRO', 'PRV', 'RSH', 'TMF')) %>%
                                mutate(param = recode(param,
                                                      'LWE' = 'Wenatchee',
                                                      'ENL' = 'Entiat',
                                                      'LMR' = 'Methow',
                                                      'OKL' = 'Okanogan',
                                                      'FST' = 'Okanogan',
                                                      'ICH' = 'Below Priest',
                                                      'JD1' = 'Below Priest',
                                                      'JDA' = 'Below Priest',
                                                      'PRH' = 'Below Priest',
                                                      'PRO' = 'Below Priest',
                                                      'PRV' = 'Below Priest',
                                                      'RSH' = 'Below Priest',
                                                      'TMF' = 'Below Priest')) %>%
                                mutate(param = factor(param,
                                                      levels = c("Wenatchee",
                                                                 "Entiat",
                                                                 "Methow",
                                                                 "Okanogan",
                                                                 "Below Priest"))) %>%
                                group_by(chain, iter, origin, param) %>%
                                summarize(across(escp,
                                                 sum),
                                          .groups = "drop") %>%
                                group_by(origin,
                                         population = param) %>%
                                summarise(mean = mean(escp),
                                          median = median(escp),
                                          # mode = estMode(escp),
                                          sd = sd(escp),
                                          .groups = 'drop') %>%
                                mutate(across(c(mean, median, sd),
                                              ~ if_else(. < 0, 0, .))) %>%
                                select(-mean) %>%
                                rename(est = median,
                                       se = sd) %>%
                                mutate(across(origin,
                                              recode,
                                              "H" = "hos",
                                              "W" = "nos"),
                                       cv = se / est) %>%
                                pivot_wider(names_from = origin,
                                            values_from = c(est, se, cv),
                                            names_glue = "{origin}_{.value}") %>%
                                mutate(total_escapement = nos_est + hos_est,
                                       phos = hos_est / (hos_est + nos_est)) %>%
                                select(population,
                                       total_escapement,
                                       starts_with("nos"),
                                       starts_with("hos"),
                                       phos)

                              # all escapement locations that year
                              escp_yr = escape_summ %>%
                                mutate(across(origin,
                                              recode,
                                              "H" = "Hatchery",
                                              "W" = "Natural")) %>%
                                select(-species, -spawn_year,
                                       -skew, -kurtosis) %>%
                                mutate(across(mean:mode,
                                              janitor::round_half_up),
                                       across(mean:mode,
                                              as.integer)) %>%
                                rename(estimate = median,
                                       se = sd) %>%
                                select(-mean, -mode)

                              # detection probabilities that year
                              detect_yr <- detect_summ %>%
                                select(node,
                                       n_tags,
                                       estimate = median,
                                       se = sd,
                                       ends_with("CI"))

                              # sex proportions that year
                              sexed_tags <- tag_summ %>%
                                filter(group %in% c("Wenatchee",
                                                    "Entiat",
                                                    "Methow",
                                                    "Okanogan",
                                                    "BelowPriest")) %>%
                                mutate(across(group,
                                              fct_recode,
                                              "Below Priest" = "BelowPriest"),
                                       across(group,
                                              fct_drop)) %>%
                                bind_rows(tag_summ %>%
                                            mutate(group = "Priest Rapids Tags")) %>%
                                group_by(population = group,
                                         origin,
                                         sex) %>%
                                summarise(n_tags = n_distinct(tag_code[!is.na(sex)]),
                                          .groups = "drop") %>%
                                filter(!is.na(sex)) %>%
                                full_join(expand(.,
                                                 population,
                                                 origin,
                                                 sex),
                                          by = c("population", "origin", "sex")) %>%
                                mutate(across(n_tags,
                                              replace_na,
                                              0))

                              sex_yr <- sexed_tags %>%
                                bind_rows(sexed_tags %>%
                                            group_by(population, sex) %>%
                                            summarize(across(n_tags,
                                                             sum),
                                                      .groups = "drop") %>%
                                            mutate(origin = "All")) %>%
                                mutate(origin = factor(origin,
                                                       levels = c("W", "H", "All")),
                                       across(origin,
                                              fct_recode,
                                              "Natural" = "W",
                                              "Hatchery" = "H")) %>%
                                mutate(group = paste(population, origin, sep = "_")) %>%
                                arrange(population, origin, sex) %>%
                                group_by(group) %>%
                                mutate(total_sexed = sum(n_tags)) %>%
                                mutate(prop = n_tags / total_sexed,
                                       prop_se = sqrt((prop * (1 - prop)) / total_sexed)) %>%
                                ungroup() %>%
                                left_join(pop_escp_yr %>%
                                            rename(total_est = total_escapement) %>%
                                            rowwise() %>%
                                            mutate(total_se = sqrt(sum(c(nos_se,
                                                                         hos_se)^2))) %>%
                                            ungroup() %>%
                                            select(population,
                                                   ends_with("est"),
                                                   ends_with("se")) %>%
                                            pivot_longer(c(ends_with("est"),
                                                           ends_with("se"))) %>%
                                            mutate(type = if_else(str_detect(name, "_se$"),
                                                                  "se",
                                                                  "est"),
                                                   origin = str_split(name, "_", simplify = T)[,1],
                                                   across(origin,
                                                          recode,
                                                          "total" = "All",
                                                          "nos" = "Natural",
                                                          "hos" = "Hatchery")) %>%
                                            select(-name) %>%
                                            pivot_wider(names_from = type,
                                                        values_from = value),
                                          by = c("population", "origin")) %>%
                                rowwise() %>%
                                mutate(across(est,
                                              ~ . * prop),
                                       se = deltamethod(~ x1 * x2,
                                                        mean = c(est, prop),
                                                        cov = diag(c(se, prop_se)^2))) %>%
                                ungroup() %>%
                                select(population,
                                       origin,
                                       sex,
                                       n_tags,
                                       total_sexed,
                                       starts_with("prop"),
                                       escp = est,
                                       escp_se = se)


                              # convert lengths from mm to cm if needed
                              if(max(tag_summ$fork_length, na.rm = T) > 100) {
                                tag_summ %<>%
                                  mutate(across(fork_length,
                                                measurements::conv_unit,
                                                from = "mm",
                                                to = "cm"))
                              }

                              size_yr <- tag_summ |>
                                filter(!group %in% c("Start", "Other")) |>
                                mutate(across(group,
                                              fct_drop)) |>
                                mutate(across(group,
                                              recode,
                                              "BelowPriest" = "Below Priest",
                                              "WellsPool" = "Wells Pool"),
                                       across(origin,
                                              recode,
                                              "H" = "Hatchery",
                                              "W" = "Natural")) |>
                                select(tag_code,
                                       population = group,
                                       origin,
                                       age,
                                       fork_length) |>
                                left_join(age_table |>
                                            select(age,
                                                   salt_water),
                                          by = "age") |>
                                filter(!is.na(fork_length),
                                       !is.na(salt_water)) |>
                                group_by(population,
                                         origin,
                                         salt_water) |>
                                summarize(across(fork_length,
                                                 list(mean = mean,
                                                      sd = sd,
                                                      n = length),
                                                 .names = "{.fn}"),
                                          .groups = "drop") %>%
                                rename(salt_water_age = salt_water) |>
                                # add data for all Priest Rapids tags
                                bind_rows(tag_summ |>
                                            mutate(across(origin,
                                                          recode,
                                                          "H" = "Hatchery",
                                                          "W" = "Natural")) |>
                                            select(tag_code,
                                                   origin,
                                                   age,
                                                   fork_length) |>
                                            left_join(age_table |>
                                                        select(age,
                                                               salt_water),
                                                      by = "age") |>
                                            filter(!is.na(fork_length),
                                                   !is.na(salt_water)) |>
                                            group_by(origin,
                                                     salt_water) |>
                                            summarize(across(fork_length,
                                                             list(mean = mean,
                                                                  sd = sd,
                                                                  n = length),
                                                             .names = "{.fn}"),
                                                      .groups = "drop") %>%
                                            rename(salt_water_age = salt_water) |>
                                            add_column(population = "Priest Rapids Tags",
                                                       .before = 0))

                              # proportions / escapement of mark groups, by population
                              mark_tag_summ = tag_summ |>
                                mutate(cwt = if_else(is.na(cwt), F,
                                                     if_else(cwt %in% c("SN", "BD"),
                                                             T, NA)),
                                       ad_clip = if_else(is.na(ad_clip) | origin == "W",
                                                         F,
                                                         if_else(ad_clip == "AD", T, NA))) |>
                                mutate(ad_clip_chr = recode(as.character(ad_clip),
                                                            "TRUE" = "AD",
                                                            "FALSE" = "AI"),
                                       cwt_chr = recode(as.character(cwt),
                                                        "TRUE" = "CWT",
                                                        "FALSE" = "noCWT")) |>
                                tidyr::unite("mark_grp", ad_clip_chr, cwt_chr, remove = T) |>
                                mutate(mark_grp = if_else(origin == "W",
                                                          "Wild",
                                                          mark_grp))

                              # proportions of each type of fish (H/W, Ad-clip, CWT combinations) past each site
                              # determine which group/population each site is in
                              site_pop_df = buildNodeOrder(parent_child) |>
                                mutate(group = if_else(node == "PRA",
                                                       "Start",
                                                       if_else(str_detect(path, 'LWE') | node %in% c("CLK"),
                                                               "Wenatchee",
                                                               if_else(str_detect(path, "ENL"),
                                                                       "Entiat",
                                                                       if_else(str_detect(path, "LMR"),
                                                                               "Methow",
                                                                               if_else(str_detect(path, "OKL") | node %in% c("FST"),
                                                                                       "Okanogan",
                                                                                       if_else(str_detect(path, "RIA", negate = T) & str_length(path > 3),
                                                                                               "BelowPriest",
                                                                                               if_else(node == "WEA",
                                                                                                       "WellsPool",
                                                                                                       "Other")))))))) |>
                                mutate(group = factor(group,
                                                      levels = c("Wenatchee",
                                                                 "Entiat",
                                                                 "Methow",
                                                                 "Okanogan",
                                                                 "BelowPriest",
                                                                 "WellsPool",
                                                                 "Start",
                                                                 "Other"))) |>
                                rename(site_code = node,
                                       site_order = node_order,
                                       population = group)


                              mark_grp_prop = mark_tag_summ |>
                                mutate(spawn_site = str_remove(final_node, "B0$"),
                                       spawn_site = str_remove(spawn_site, "A0$"),
                                       spawn_site = recode(spawn_site,
                                                           "S" = "SA0")) |>
                                left_join(site_pop_df |>
                                            select(end_loc = site_code,
                                                   population,
                                                   path) |>
                                            separate(path,
                                                     into = paste("site", 1:8, sep = "_")) |>
                                            pivot_longer(starts_with('site_'),
                                                         names_to = "node_order",
                                                         values_to = "site_code") |>
                                            filter(!is.na(site_code)) |>
                                            select(-node_order),
                                          by = c("spawn_site" = "end_loc"),
                                          multiple = "all") |>
                                group_by(origin,
                                         population = group,
                                         site_code, ad_clip, cwt, mark_grp) |>
                                summarise(n_tags = n_distinct(tag_code),
                                          .groups = "drop") |>
                                right_join(crossing(site_code = union(parent_child$parent, parent_child$child),
                                                    expand(mark_tag_summ, nesting(origin, ad_clip, cwt, mark_grp))) |>
                                             left_join(site_pop_df |>
                                                         select(site_code,
                                                                population),
                                                       by = "site_code") |>
                                             relocate(population, .after = 1),
                                           by = c("origin", "population", "site_code", "ad_clip", "cwt", "mark_grp")) |>
                                arrange(site_code, mark_grp) |>
                                mutate(across(n_tags,
                                              replace_na,
                                              0)) |>
                                group_by(origin,
                                         site_code) |>
                                mutate(tot_tags = sum(n_tags),
                                       prop = n_tags / tot_tags,
                                       # using normal approximation
                                       prop_se = sqrt((prop * (1 - prop))/tot_tags)) |>
                                ungroup() |>
                                arrange(population,
                                        site_code,
                                        origin,
                                        ad_clip,
                                        cwt)

                              # generate posterior samples of mark proportions
                              escape_post <- escape_post |>
                                mutate(n_iter = max(iter),
                                       iter = (chain - 1) * n_iter + iter)
                              set.seed(6)
                              n_iter = max(escape_post$iter)
                              prop_samps = mark_grp_prop |>
                                filter(tot_tags > 0) |>
                                group_by(origin,
                                         population,
                                         site_code,
                                         tot_tags) |>
                                nest() |>
                                mutate(samp = map(data,
                                                  .f = function(x) {
                                                    rmultinom(n_iter, sum(x$n_tags), x$prop) |>
                                                      set_colnames(1:n_iter) |>
                                                      set_rownames(x$mark_grp) |>
                                                      as_tibble(rownames = 'mark_grp') |>
                                                      pivot_longer(-1,
                                                                   names_to = "iter",
                                                                   values_to = "n_tags") |>
                                                      mutate(across(iter,
                                                                    as.integer)) |>
                                                      group_by(iter) |>
                                                      mutate(prop = n_tags / sum(n_tags)) |>
                                                      arrange(iter, mark_grp) |>
                                                      ungroup()
                                                  })) |>
                                ungroup() |>
                                select(-data) |>
                                unnest(samp) |>
                                left_join(mark_grp_prop |>
                                            select(-n_tags,
                                                   -starts_with("prop")),
                                          by = c("population",
                                                 "site_code",
                                                 "origin", "tot_tags", "mark_grp")) |>
                                bind_rows(mark_grp_prop |>
                                            filter(tot_tags == 0) |>
                                            crossing(iter = 1:n_iter) |>
                                            select(-prop_se)) |>
                                arrange(population,
                                        site_code,
                                        origin, mark_grp, iter) |>
                                select(iter, any_of(names(mark_grp_prop)))

                              # posterior samples
                              mark_post = escape_post |>
                                inner_join(prop_samps,
                                           by = c("iter", "origin",
                                                  "param" = "site_code"),
                                           multiple = "all") %>%
                                mutate(across(c(prop,
                                                n_tags,
                                                tot_tags),
                                              replace_na,
                                              0)) |>
                                mutate(prop = if_else(value == 0,
                                                      0,
                                                      prop)) |>
                                mutate(n_fish = escp * prop)

                              mark_grp_yr = mark_post |>
                                group_by(site_code = param,
                                         origin,
                                         ad_clip,
                                         cwt,
                                         mark_grp) |>
                                summarise(across(n_fish,
                                                 list(mean = mean,
                                                      median = median,
                                                      # mode = DABOM::estMode,
                                                      se = sd,
                                                      skew = moments::skewness,
                                                      kurtosis = moments::kurtosis,
                                                      lowerCI = ~ coda::HPDinterval(coda::as.mcmc(.x))[,1],
                                                      upperCI = ~ coda::HPDinterval(coda::as.mcmc(.x))[,2]),
                                                 # na.rm = T,
                                                 .names = "{.fn}"),
                                          .groups = "drop") |>
                                # mutate(across(mode,
                                #               ~ if_else(. < 0, 0, .))) |>
                                left_join(mark_grp_prop |>
                                            rename(proportion = prop),
                                          by = c("site_code",
                                                 "origin", "ad_clip", "cwt", "mark_grp")) |>
                                select(population,
                                       site_code,
                                       origin:mark_grp,
                                       n_tags:proportion,
                                       prop_se,
                                       everything()) |>
                                # mutate(across(population,
                                #               factor,
                                #               levels = levels(tag_summ$group)),
                                #        across(population,
                                #               fct_drop)) |>
                                arrange(population,
                                        site_code,
                                        mark_grp) |>
                                mutate(across(origin,
                                              fct_recode,
                                              "Natural" = "W",
                                              "Hatchery" = "H"),
                                       across(population,
                                              factor,
                                              levels = c('Wenatchee',
                                                         'Entiat',
                                                         'Methow',
                                                         'Okanogan'))) |>
                                select(population:prop_se,
                                       estimate = median,
                                       se,
                                       lowerCI,
                                       upperCI)

                              mark_grp_pop_yr <- site_pop_df |>
                                group_by(population) |>
                                filter(str_length(path) == min(str_length(path))) |>
                                ungroup() |>
                                filter(population %in% c('Wenatchee',
                                                         'Entiat',
                                                         'Methow',
                                                         'Okanogan')) |>
                                arrange(population,
                                        site_code) |>
                                select(population,
                                       site_code) |>
                                left_join(mark_grp_yr,
                                          by = c("population",
                                                 "site_code"),
                                          multiple = "all") |>
                                group_by(population,
                                         origin,
                                         ad_clip,
                                         cwt,
                                         mark_grp) %>%
                                summarize(across(c(n_tags, tot_tags,
                                                   estimate,
                                                   ends_with("CI")),
                                                 sum,
                                                 na.rm = T),
                                          across(se,
                                                 ~ sqrt(sum(.^2, na.rm = T))),
                                          .groups = "drop") |>
                                group_by(population,
                                         origin) |>
                                mutate(tot_tags = sum(n_tags),
                                       prop = n_tags / tot_tags,
                                       # using normal approximation
                                       prop_se = sqrt((prop * (1 - prop))/tot_tags)) |>
                                ungroup() |>
                                select(any_of(names(mark_grp_yr)))

                              res_list <- list(tag_summ = tag_summ,
                                               age_prop_yr = age_prop_yr,
                                               pop_escp_yr = pop_escp_yr,
                                               escp_yr = escp_yr,
                                               detect_yr = detect_yr,
                                               sex_yr = sex_yr,
                                               size_yr = size_yr,
                                               mark_grp_yr = mark_grp_yr,
                                               mark_grp_pop_yr = mark_grp_pop_yr)

                              rm(age_prop_yr,
                                 pop_escp_yr,
                                 escp_yr,
                                 detect_yr,
                                 sexed_tags,
                                 sex_yr,
                                 mark_grp_yr,
                                 mark_grp_pop_yr,
                                 mark_grp_prop,
                                 mark_post,
                                 mark_tag_summ,
                                 site_pop_df,
                                 prop_samps)

                              return(res_list)
                            }))

pop_escp <- dabom_est |>
  mutate(res = map(all_results,
                   "pop_escp_yr")) |>
  select(-dam_cnt_name,
         -all_results) |>
  unnest(res) |>
  makeTableNms()

all_escp <- dabom_est |>
  mutate(res = map(all_results,
                   "escp_yr")) |>
  select(-dam_cnt_name,
         -all_results) |>
  unnest(res) |>
  makeTableNms()

detect_df <- dabom_est |>
  mutate(res = map(all_results,
                   "detect_yr")) |>
  select(-dam_cnt_name,
         -all_results) |>
  unnest(res) |>
  makeTableNms()

tag_df <- dabom_est |>
  mutate(res = map(all_results,
                   "tag_summ")) |>
  select(-dam_cnt_name,
         -all_results) |>
  unnest(res) |>
  select(run_year:tag_code,
         species,
         record_id,
         tag_other,
         sex:age,
         ad_clip,
         cwt,
         trap_date,
         final_node:tag_detects,
         path) |>
  makeTableNms()

sex_prop <- dabom_est |>
  mutate(res = map(all_results,
                   "sex_yr")) |>
  select(-dam_cnt_name,
         -all_results) |>
  unnest(res) |>
  makeTableNms()

size_df <- dabom_est |>
  mutate(res = map(all_results,
                   "size_yr")) |>
  select(-dam_cnt_name,
         -all_results) |>
  unnest(res) |>
  mutate(across(
    population,
    ~ factor(.,
             levels = c('Wenatchee',
                        'Entiat',
                        'Methow',
                        'Okanogan',
                        "Below Priest",
                        "Wells Pool",
                        "Priest Rapids Tags"))
  )) |>
  arrange(spawn_year,
          population,
          origin,
          salt_water_age)

size_tab <- expand(size_df,
                   nesting(run_year, spawn_year),
                   population, origin, salt_water_age) |>
  left_join(size_df,
            by = c("run_year",
                   "spawn_year",
                   "population",
                   "origin",
                   "salt_water_age")) |>
  mutate(across(salt_water_age,
                ~ paste0("Salt-", .))) |>
  rename(Mean = mean,
         SD = sd,
         N = n) |>
  group_by(spawn_year,
           population,
           origin) |>
  mutate(sample_size = sum(N, na.rm = T)) |>
  ungroup() |>
  pivot_wider(names_from = salt_water_age,
              values_from = c(Mean, SD, N),
              names_glue = "{salt_water_age}_{.value}",
              names_vary = "slowest") |>
  makeTableNms()



age_prop <- dabom_est %>%
  mutate(res = map(all_results,
                   "age_prop_yr")) %>%
  select(-dam_cnt_name,
         -all_results) %>%
  unnest(res) |>
  mutate(across(
    population,
    ~ factor(.,
             levels = c('Wenatchee',
                        'Entiat',
                        'Methow',
                        'Okanogan',
                        "Below Priest",
                        "Wells Pool",
                        "Priest Rapids Tags"))
  ))

tot_age_prop <- age_prop |>
  filter(age_type == "total_age") |>
  mutate(across(age_label,
                factor,
                levels = 2:11)) |>
  select(-c(age_type, n)) |>
  pivot_wider(names_from = age_label,
              names_expand = T,
              values_from = prop,
              values_fill = 0) |>
  makeTableNms()

sw_age_prop <- age_prop |>
  filter(age_type == "salt_water") |>
  mutate(across(age_label,
                factor,
                levels = 0:4)) |>
  select(-c(age_type, n)) |>
  pivot_wider(names_from = age_label,
              names_expand = T,
              values_from = prop,
              values_fill = 0) |>
  makeTableNms()

fw_age_prop <- age_prop |>
  filter(age_type == "fresh_water") |>
  mutate(across(age_label,
                ~ paste0(., ".X")),
         across(age_label,
                factor,
                levels = paste0(1:5, ".X"))) |>
  select(-c(age_type, n)) |>
  pivot_wider(names_from = age_label,
              names_expand = T,
              values_from = prop,
              values_fill = 0) |>
  makeTableNms()


mark_grp_df <- dabom_est %>%
  mutate(res = map(all_results,
                   "mark_grp_yr")) %>%
  select(-dam_cnt_name,
         -all_results) %>%
  unnest(res) |>
  makeTableNms()

mark_grp_pop_df <- dabom_est %>%
  mutate(res = map(all_results,
                   "mark_grp_pop_yr")) %>%
  select(-dam_cnt_name,
         -all_results) %>%
  unnest(res) |>
  makeTableNms()

#-----------------------------------------------------------------
# pull together estimates of sex call error rates at Priest
# check for duplicate tags
read_excel(here("data","STHD_UC Brood Collections_2011 to current.xlsx"),
           sheet = "Brood Collected_PIT Tagged Only") |>
  clean_names() |>
  distinct() |>
  # one tag has two records; choose the one that matches PTAGIS recapture details
  filter(!(recaptured_pit == "3DD.0077DA4CC6" & sex_final == "M")) |>
  # filter(recaptured_pit %in% recaptured_pit[duplicated(recaptured_pit)]) |>
  # arrange(recaptured_pit) |>
  # group_by(recaptured_pit) |>
  # mutate(n_final_sex = n_distinct(sex_final)) |>
  # ungroup() |>
  # filter(n_final_sex > 1)
  # as.data.frame()
  select(spawn_year,
         tag_code = recaptured_pit,
         sex_final) |>
  distinct() |>
  filter(tag_code %in% tag_code[duplicated(tag_code)]) |>
  arrange(tag_code) #|>
  inner_join(tag_df |>
               clean_names() |>
               select(tag_code))

# estimate error rate for each sex
sex_err_rate <- tag_df |>
  clean_names() |>
  select(spawn_year,
         tag_code,
         sex_field = sex) |>
  inner_join(read_excel(here("data",
                               "STHD_UC Brood Collections_2011 to current.xlsx"),
                        sheet = "Brood Collected_PIT Tagged Only") |>
               clean_names() |>
               # one tag has two records; choose the one that matches PTAGIS recapture details
               filter(!(recaptured_pit == "3DD.0077DA4CC6" & sex_final == "M")) |>
               rename(tag_code = recaptured_pit) |>
               select(spawn_year,
                      tag_code,
                      sex_final) |>
               distinct(),
             by = c("spawn_year",
                    "tag_code")) |>
  filter(!is.na(sex_final),
         !is.na(sex_field)) |>
  mutate(agree = if_else(sex_field == sex_final,
                         T, F)) |>
  group_by(spawn_year,
           sex = sex_field) |>
  summarize(n_tags = n_distinct(tag_code),
            n_true = sum(agree),
            n_false = sum(!agree),
            .groups = "drop") |>
  mutate(binom_ci = map2(n_false,
                         n_tags,
                         .f = function(x, y) {
                           DescTools::BinomCI(x, y) |>
                             as_tibble()
                         })) |>
  unnest(binom_ci) |>
  clean_names() |>
  rename(perc_false = est,
         lowerci = lwr_ci,
         upperci = upr_ci) |>
  mutate(perc_se = sqrt((perc_false * (1 - perc_false)) / n_tags)) |>
  relocate(perc_se,
           .after = "perc_false")

# re-calculate proportion of each sex in each population each year
sex_prop <- dabom_est |>
  mutate(res = map(all_results,
                   "sex_yr")) |>
  select(-dam_cnt_name,
         -all_results) |>
  unnest(res) |>
  mutate(across(origin,
                as_factor)) |>
  select(run_year:total_sexed) |>
  rename(priest_tags = n_tags) |>
  left_join(sex_err_rate |>
              select(spawn_year,
                     sex,
                     starts_with("perc")),
            by = c("spawn_year", "sex")) |>
  pivot_wider(names_from = sex,
              values_from = c(priest_tags,
                              perc_false,
                              perc_se)) |>
  mutate(true_male = priest_tags_M - (priest_tags_M * perc_false_M) + (priest_tags_F * perc_false_F),
         true_female = priest_tags_F - (priest_tags_F * perc_false_F) + (priest_tags_M * perc_false_M)) |>
  # mutate(across(starts_with("true"),
  #               janitor::round_half_up)) |>
  rowwise() |>
  mutate(true_m_se = msm::deltamethod(~ x1 - (x1 * x2) + (x3 * x4),
                                      mean = c(priest_tags_M,
                                               perc_false_M,
                                               priest_tags_F,
                                               perc_false_F),
                                      cov = diag(c(0,
                                                   perc_se_M,
                                                   0,
                                                   perc_se_F)^2)),
         true_f_se = msm::deltamethod(~ x1 - (x1 * x2) + (x3 * x4),
                                      mean = c(priest_tags_F,
                                               perc_false_F,
                                               priest_tags_M,
                                               perc_false_M),
                                      cov = diag(c(0,
                                                   perc_se_F,
                                                   0,
                                                   perc_se_M)^2))) |>
  mutate(n_sexed = true_male + true_female,
         across(n_sexed,
                round_half_up),
         prop_m = true_male / (true_male + true_female),
         prop_se = msm::deltamethod(~ x1 / (x1 + x2),
                                    mean = c(true_male,
                                             true_female),
                                    cov = diag(c(true_m_se,
                                                 true_f_se)^2))) |>
  ungroup() |>
  select(run_year:origin,
         total_sexed,
         contains("priest"),
         n_sexed,
         starts_with("true"),
         starts_with("prop")) |>
  mutate(prop_f = 1 - prop_m) |>
  # filter(total_sexed != n_sexed)
  clean_names() |>
  pivot_longer(cols = c(starts_with("priest"),
                        starts_with("true"),
                        starts_with("prop"))) |>
  mutate(sex = if_else(str_detect(name, "_m"),
                       "M", "F"),
         across(name,
                str_remove,
                "_male"),
         across(name,
                str_remove,
                "_female"),
         across(name,
                str_remove,
                "_m"),
         across(name,
                str_remove,
                "_f"),
         across(name,
                str_remove,
                "_M"),
         across(name,
                str_remove,
                "_F"),
         across(name,
                str_replace,
                "true",
                "n_tags")) |>
  pivot_wider(names_from = "name",
              values_from = "value") |>
  arrange(spawn_year,
          population,
          origin,
          sex) |>
  fill(prop_se,
       .direction = "down") |>
  select(run_year:origin,
         sex,
         priest_tags,
         n_tags,
         total_sexed = n_sexed,
         # total_sexed,
         starts_with("prop")) |>
  mutate(across(c(n_tags,
                  total_sexed),
                round_half_up)) |>
  left_join(pop_escp |>
              clean_names() |>
              rename(total_est = total_escapement) %>%
              rowwise() %>%
              mutate(total_se = sqrt(sum(c(nos_se,
                                           hos_se)^2))) %>%
              ungroup() %>%
              select(run_year,
                     spawn_year,
                     population,
                     ends_with("est"),
                     ends_with("se")) %>%
              pivot_longer(c(ends_with("est"),
                             ends_with("se"))) %>%
              mutate(type = if_else(str_detect(name, "_se$"),
                                    "se",
                                    "est"),
                     origin = str_split(name, "_", simplify = T)[,1],
                     across(origin,
                            recode,
                            "total" = "All",
                            "nos" = "Natural",
                            "hos" = "Hatchery")) %>%
              select(-name) %>%
              pivot_wider(names_from = type,
                          values_from = value),
            by = c("run_year", "spawn_year", "population", "origin")) |>
  rowwise() |>
  mutate(escp = est * prop,
         escp_se = msm::deltamethod(~ x1 * x2,
                                    mean = c(est,
                                             prop),
                                    cov = diag(c(se, prop_se)^2))) |>
  ungroup() |>
  select(-est,
         -se) |>
  makeTableNms()


#-----------------------------------------------------------------
# generate table of escapement estimates at Priest
priest_df = tibble(spawn_year = 2011:max_yr) |>
  mutate(prd_df = map(spawn_year,
                      .f = function(yr) {
                        # load data
                        load(here("analysis/data/derived_data/estimates",
                                  "PriestRapids",
                                  paste0("UC_Sthd_DABOM_", yr, ".rda")))

                        tag_num <- bio_df |>
                          group_by(origin) |>
                          summarize(n_tags = n_distinct(tag_code),
                                    .groups = "drop") |>
                          mutate(total_tags = sum(n_tags))

                        res <- dam_escp_df |>
                          filter(dam == "PriestRapids") |>
                          mutate(run_year = year - 1) |>
                          select(run_year,
                                 # spawn_year = year,
                                 win_cnt,
                                 origin,
                                 contains("reasc_rate"),
                                 contains("priest_cnt")) |>
                          mutate(total = sum(priest_cnt),
                                 total_se = sqrt(sum(priest_cnt_se^2))) |>
                          left_join(tag_num,
                                    by = "origin")
                        return(res)
                      })) |>
  unnest(prd_df) |>
  relocate(run_year,
           .before = 1) |>
  mutate(across(origin,
                recode,
                "H" = "Hatchery",
                "W" = "Natural")) |>
  makeTableNms() |>
  rename(`PRD Ladder Count` = `Win Cnt`,
         `Adult Reascension Adjustment Rate` = `Reasc Rate`,
         `Reascension Rate SE` = `Reasc Rate SE`,
         `Priest Est` = `Priest Cnt`,
         `Priest SE` = `Priest Cnt SE`)


# now a table for years when Rock Island was used
rock_isl_df = priest_df |>
  select(run_year = `Run Year`,
         spawn_year = `Spawn Year`) |>
  filter(!spawn_year %in% c(2011:2015, 2018)) |>
  distinct() |>
  mutate(ria_df = map(spawn_year,
                      .f = function(yr) {
                        # load data
                        load(here("analysis/data/derived_data/estimates",
                                  "RockIsland",
                                  paste0("UC_Sthd_DABOM_", yr, ".rda")))

                        res <- dam_escp_df |>
                          filter(dam == "RockIsland") |>
                          select(origin,
                                 #contains("reasc_rate"),
                                 contains("priest_cnt"),
                                 dam) |>
                          mutate(total = sum(priest_cnt),
                                 total_se = sqrt(sum(priest_cnt_se^2))) |>
                          relocate(dam,
                                   .after = "total_se")
                        return(res)
                      })) |>
  unnest(ria_df) |>
  mutate(across(origin,
                recode,
                "H" = "Hatchery",
                "W" = "Natural"),
         across(dam,
                recode,
                "RockIsland" = "Rock Island")) |>
  makeTableNms() |>
  rename(#`PRD Ladder Count` = `Win Cnt`,
    # `Adult Reascension Adjustment Rate` = `Reasc Rate`,
    # `Reascension Rate SE` = `Reasc Rate_se`,
    `Modelled Priest Est` = `Priest Cnt`,
    `Modelled Priest SE` = `Priest Cnt SE`,
    `Modelled Total` = Total,
    `Modelled Total SE` = `Total SE`,
    `Dam Method` = Dam)


dam_cnt_tab <- priest_df |>
  left_join(rock_isl_df,
            by = c("Run Year",
                   "Spawn Year",
                   "Origin")) |>
  relocate(c(`N Tags`,
             `Total Tags`),
           .after = "Dam Method") |>
  mutate(`Tag Rate` = if_else(is.na(`Dam Method`),
                              `N Tags` / `Priest Est`,
                              `N Tags` / `Modelled Priest Est`),
         `Total Tag Rate` = if_else(is.na(`Dam Method`),
                                    `Total Tags` / Total,
                                    `Total Tags` / `Modelled Total`))

#-----------------------------------------------------------------
# put together in a list to save / write to Excel
# make some decisions about rounding
save_list <- list(
  "Redd Abundance" = redd_spwn |>
    mutate(across(ends_with("Year"),
                  as.integer),
           across(c(`Redd Est`,
                    `NOS Est`,
                    `HOS Est`),
                  round_half_up)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Run Escp by Origin" = pop_escp |>
    mutate(across(ends_with("Year"),
                  as.integer),
           across(c(`Total Escapement`,
                    ends_with("Est")),
                  round_half_up)) |>
    mutate(`Total Escapement` = `NOS Est` + `HOS Est`) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Spawning Escp by Origin" = spwn_df |>
    bind_rows(pop_escp |>
                filter(Population == "Entiat") |>
                rename(`Total Spawners` = `Total Escapement`) |>
                mutate(Stream = as.factor("Total")) |>
                rowwise() |>
                mutate(`pHOS SE` = deltamethod(~ x1 / (x1 + x2),
                                               mean = c(`HOS Est`,
                                                        `NOS Est`),
                                               cov = diag(c(`HOS SE`,
                                                            `NOS SE`)^2))) |>
                ungroup() |>
                mutate(across(Population,
                              factor,
                              levels = levels(pop_escp$Population)))) |>
    arrange(`Spawn Year`,
            Population,
            Stream) |>
    mutate(across(ends_with("Year"),
                  as.integer),
           across(c(`Total Spawners`,
                    ends_with("Est")),
                  round_half_up)) |>
    mutate(`Total Spawners` = `NOS Est` + `HOS Est`) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Run Escp All Locations" = all_escp |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Detection Probability" = detect_df |>
    mutate(across(c(ends_with("Year"),
                    `N Tags`),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Tag Summary" = tag_df |>
    mutate(across(ends_with("Year"),
                  as.integer)),

  "Sex Error Rates" = sex_err_rate |>
    makeTableNms() |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Run Escp Sex" = sex_prop |>
    mutate(across(ends_with("Year"),
                  as.integer),
           across(Escp,
                  ~ as.integer(round_half_up(.)))) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "SaltAge at Maturity" = sw_age_prop |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "FWAge at Maturity" = fw_age_prop |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "TotalAge at Maturity" = tot_age_prop |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Size at Maturity" = size_tab |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Priest Rapids Escapement" = dam_cnt_tab |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(c("PRD Ladder Count",
                    ends_with("Priest Est"),
                    ends_with("Total")),
                  ~ as.integer(round_half_up(.)))) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Mark Group Population" = mark_grp_pop_df |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(Estimate,
                  ~ as.integer(round_half_up(.)))) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Mark Group All Locations" = mark_grp_df |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(Estimate,
                  ~ as.integer(round_half_up(.)))) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3))

)

#-----------------------------------------------------------------
# actually save entire file
write_xlsx(x = save_list,
           path = here("output",
                         "UC_STHD_Model_Output.xlsx"))

#--------------------------------------------------------------
# read in previous estimates, add latest year to them

# what year are we overwriting or adding?
yr = 2022

library(readxl)
output_path <- here("output",
                      "UC_STHD_Model_Output.xlsx")

tab_nms <- excel_sheets(output_path)

previous_output <- as.list(tab_nms) |>
  rlang::set_names() %>%
  map(.f = function(x) {
    read_excel(output_path,
               sheet = x) %>%
      filter(`Spawn Year` != yr)
  })

lst_one_yr <- save_list |>
  map(.f = function(x) {
    x |>
      filter(`Spawn Year` == yr)
  })



identical(names(previous_output), names(lst_one_yr))

save_list = vector("list",
                   length = length(previous_output))
names(save_list) = names(previous_output)

previous_output[[3]]$'HOS CV' <- as.numeric(previous_output[[3]]$'HOS CV')

for(i in 1:length(previous_output)) {
  save_list[[i]] <- previous_output[[i]] |>
    bind_rows(lst_one_yr[[i]]) %>%
    arrange(`Spawn Year`)
}

write_xlsx(x = save_list,
           path = here("output",
                         "UC_STHD_Model_Output_new.xlsx"))
