# ...............................................................................
# ASSEGNO PER IL LAVORO - 04auto Heterogeneity through subgroup analyses
# Author: Álvaro F. Junquera
# ...............................................................................

library(tidyverse)
library(purrr)
library(rdrobust)

# 0. Reading -------------
indi_ns_ss1_c <- readRDS("intermediate/script02/indi_ns_ss1_c.RDS")
indi_ns_ss2_c <- readRDS("intermediate/script02/indi_ns_ss2_c.RDS")

# 1. Base function to run rdrobust -----------
run_rd <- function(y, x) {
  rdrobust(y = y, x = x,
           kernel = "triangular",
           c = 0, p = 1, bwselect = "mserd",
           cluster = NULL)
}

# 2. Function to run rdrobust for all outcomes -----------
run_rd_outcomes <- function(df, score_var, outcomes) {
  outcomes %>%
    set_names() %>%
    map(~ run_rd(y = df[[.x]], x = df[[score_var]]))
}


# 3. Function to run subgroup analyses -----------
analyze_subgroups <- function(data, subgroup_var, score_var, outcomes,
                              continuous_splits = NULL) {

  # Case A: categorical variable
  if (is.factor(data[[subgroup_var]]) || is.character(data[[subgroup_var]])) {
    categories <- unique(na.omit(data[[subgroup_var]]))

    res <- categories %>%
      set_names() %>%
      map(~ data %>%
            filter(.data[[subgroup_var]] == .x) %>%
            run_rd_outcomes(score_var, outcomes))
  }

  # Case B: continuous variable (requires cutoffs)
  else if (is.numeric(data[[subgroup_var]])) {
    if (is.null(continuous_splits)) {
      stop("You need to specify cutoffs for continuous variables")
    }

    res <- continuous_splits %>%
      map(~ {
        above <- data %>%
          filter(.data[[subgroup_var]] >= .x) %>%
          run_rd_outcomes(score_var, outcomes)

        below <- data %>%
          filter(.data[[subgroup_var]] < .x) %>%
          run_rd_outcomes(score_var, outcomes)

        set_names(list(above, below), nm = c(paste0("≥", .x), paste0("<", .x)))
      }) %>%
      flatten()

  }

  else {
    stop("Unknown type of variable")
  }

  return(res)
}

# 4. Extract results from the product of rdrobust ------------
rd_to_df <- function(rd_obj, outcome, subgroup, level) {
  tibble::tibble(
    subgroup_var = subgroup,
    subgroup_lvl = level,
    outcome  = outcome,
    coef = rd_obj$coef[1],
    se = rd_obj$se[1],
    ci_lower = rd_obj$ci[1, 1],
    ci_upper = rd_obj$ci[1, 2],
    p_value = rd_obj$pv[1],
    stars = case_when(
      rd_obj$pv[1] < 0.01 ~ "***",
      rd_obj$pv[1] < 0.05 ~ "**",
      rd_obj$pv[1] < 0.1 ~ "*",
      T ~ ""),
    n_eff = sum(rd_obj$N_h)
  )
}

# Convertir lista jerárquica de resultados en data.frame tidy
rdlist_to_df <- function(rdlist, subgroup_var) {
  # Extraemos los nombres de los niveles
  levels <- names(rdlist)

  # Iteramos sobre cada nivel
  purrr::map2_dfr(rdlist, levels, function(res_outcomes, lvl) {
    # Para cada outcome dentro del nivel
    outcome_names <- names(res_outcomes)

    purrr::map2_dfr(res_outcomes, outcome_names, function(rd_obj, outcome) {
      rd_to_df(rd_obj, outcome, subgroup_var, lvl)
    })
  })
}


# Heterogeneity analyses -------------

# Possible covariates:
# - Participant covariates:
#   - Categorical: province of residence, sex, level of education, disability,
#     last occupation, last sector, last type of contract
#   - Continuous: age, duration of the last contract (days)
# - Provider covariate: provider

outcomes <- c("post_interval6", "post_interval712",
              "post_interval1318", "post_interval1924",
              "jshours", "attiv_form_ore_prev")

## Categorical variables -----------
### By sex -----------
res_genere <- analyze_subgroups(
  data = indi_ns_ss1_c,
  subgroup_var = "genere",
  score_var = "scoringD1_0",
  outcomes = outcomes
)

df_genere <- rdlist_to_df(res_genere, subgroup_var = "genere")
df_genere |>
  filter(outcome %in% c("jshours", "attiv_form_ore_prev"))


### By level of education ---------
indi_ns_ss1_c <- indi_ns_ss1_c %>%
  mutate(educa = case_when(studio2_grouped2 %in% c("1_primary", "2_lowersec") ~ "Edu_1_2",
                           studio2_grouped2 %in% c("344_uppersecG",
                                                   "353_uppersecP_3y",
                                                   "353_uppersecP_4y_prof",
                                                   "353_uppersecP_4y_tec") ~ "Edu_3",
                           studio2_grouped2 %in% c("660_terBachy+_lev1",
                                                   "760+_terMastery+_lev2") ~ "Edu_6_7",
                           TRUE ~ "Edu_NA"))

## Run
res_educa <- analyze_subgroups(
  data = indi_ns_ss1_c,
  subgroup_var = "educa",
  score_var = "scoringD1_0",
  outcomes = outcomes
)

df_educa <- rdlist_to_df(res_educa, subgroup_var = "educa")

df_educa |>
  filter(outcome %in% c("jshours", "attiv_form_ore_prev"))

### By last occupation --------------
indi_ns_ss1_c <- indi_ns_ss1_c %>%
  mutate(lastocc = str_sub(edu_last_contract, 1, 1)) %>%
  mutate(lastoccu = case_when(lastocc == "9" | is.na(lastocc) ~ "NA",
                              TRUE ~ lastocc)) %>%
  select(-lastocc)

res_lastoccu <- analyze_subgroups(
  data = indi_ns_ss1_c,
  subgroup_var = "lastoccu",
  score_var = "scoringD1_0",
  outcomes = outcomes
)

df_lo <- rdlist_to_df(res_lastoccu, subgroup_var = "lastoccu")

df_lo |>
  filter(outcome %in% c("jshours", "attiv_form_ore_prev"))




### By provider -----------
provs_atleast300 <- indi_ns_ss1_c %>%
  janitor::tabyl(axl_ente) %>%
  filter(n > 299) %>%
  arrange(desc(n)) %>%
  pull(axl_ente)

indi_ns_ss1_c <- indi_ns_ss1_c %>%
  mutate(provider = if_else(axl_ente %in% provs_atleast300, axl_ente, "SmallProvider"))

## Run
res_prov <- analyze_subgroups(
  data = indi_ns_ss1_c,
  subgroup_var = "provider",
  score_var = "scoringD1_0",
  outcomes = outcomes
)

df_prov <- rdlist_to_df(res_prov, subgroup_var = "provider")

df_prov %>%
  filter(n_eff > 200) %>%
  filter(subgroup_lvl != "SmallProvider") %>%
  print(n = "all")

df_prov |>
  filter(n_eff > 200) %>%
  filter(subgroup_lvl != "SmallProvider") %>%
  filter(outcome %in% c("jshours", "attiv_form_ore_prev"))


## Continuous variables -------------
### By age ---------------
res_eta <- analyze_subgroups(
  data = indi_ns_ss1_c,
  subgroup_var = "eta",
  score_var = "scoringD1_0",
  outcomes = outcomes,
  continuous_splits = 50
)

df_edad <- rdlist_to_df(res_eta, subgroup_var = "eta")
df_edad
