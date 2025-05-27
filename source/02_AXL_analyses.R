# ...............................................................................
# ASSEGNO PER IL LAVORO - 02 Main analyses
# Author: Álvaro F. Junquera (UAB)
# ...............................................................................

library(tidyverse)
library(data.table)
library(tidytable)
library(lubridate)
library(descr)
library(janitor)

library(readxl)
library(openxlsx)

library(rdrobust)
library(QTE.RD)
library(rddensity)
library(binsreg)

# Place "rd.categorical.R" in the Project file and then run renv::install(normalizePath("rd.categorical"))
# It is a package developed by Ke-Li Xu (2017) and not available in CRAN
library(rd.categorical)

library(modelsummary)
library(glue)

# Table 2 is produced with 01_AXL_wrangling.R
# Table 3 is produced manually by joining those estimates with the highest AMSE estimate


# 1. Reading data -------------
indi_ns_ss1 <- readRDS("intermediate/script01/indi_ns_ss1_190225.RDS")
indi_ns_ss2 <- readRDS("intermediate/script01/indi_ns_ss2_190225.RDS")
longest <- readRDS("intermediate/script01/longest_190225.RDS")

# Some functions -----------

create_ci <- function(x) {
  paste0("(", round(x[1], 3), ", ", round(x[2], 3), ")")
}

# Design matrix

create_ci_vec <- function(x, fila) {
  paste0("(", round(x[fila, 1], 3), ", ", round(x[fila, 2], 3), ")")
}


rdcate_multinom <- function(score, outcome) {

  dm_D <- model.matrix(~ score + outcome - 1)

  colnames(dm_D)[1] <- "scoring"

  ilc_ord <- as.data.frame(dm_D) %>% arrange(scoring) # de forma ascendente
  ilc_ordm <- as.matrix(ilc_ord)

  miJ <- length(ilc_ordm[1, ]) - 2

  rd_cat_est99 <- rd.mnl(
    DAT = ilc_ordm, c = 0,
    H0_t = matrix(data = 0, nrow = miJ, ncol = 1),
    H0_R = diag(miJ),
    H0_r = matrix(0, miJ, 1),
    level = 0.99
  )

  rd_cat_est95 <- rd.mnl(
    DAT = ilc_ordm, c = 0,
    H0_t = matrix(data = 0, nrow = miJ, ncol = 1),
    H0_R = diag(miJ),
    H0_r = matrix(0, miJ, 1),
    level = 0.95
  )

  rd_cat_est90 <- rd.mnl(
    DAT = ilc_ordm, c = 0,
    H0_t = matrix(data = 0, nrow = miJ, ncol = 1),
    H0_R = diag(miJ),
    H0_r = matrix(0, miJ, 1),
    level = 0.90
  )

  significance90 <- rep(NA, miJ)
  significance95 <- rep(NA, miJ)
  significance99 <- rep(NA, miJ)

  for(f in 1:length(significance90)) {
    significance90[f] <- if (data.table::between(0,
                                                 rd_cat_est90$ci_rob[f, 1],
                                                 rd_cat_est90$ci_rob[f, 2]) == T) {
      F
    } else {
      T
    }

  }

  for(f in 1:length(significance95)) {
    significance95[f] <- if (data.table::between(0,
                                                 rd_cat_est95$ci_rob[f, 1],
                                                 rd_cat_est95$ci_rob[f, 2]) == T) {
      F
    } else {
      T
    }
  }

  for(f in 1:length(significance99)) {
    significance99[f] <- if (data.table::between(0,
                                                 rd_cat_est99$ci_rob[f, 1],
                                                 rd_cat_est99$ci_rob[f, 2]) == T) {
      F
    } else {
      T
    }
  }

  confidenceinterval95 <- rep(NA, miJ)

  for(g in 1:length(confidenceinterval95)) {
    confidenceinterval95[g] <- create_ci(x = rd_cat_est95$ci_rob[g,])
  }

  tabla_est <- data.frame(
    pointest = round(rd_cat_est95$ATE, 3),
    ci95 = confidenceinterval95,
    signif90 = significance90,
    signif95 = significance95,
    signif99 = significance99
  )

  tabla_est <- tabla_est %>%
    mutate(stars = case_when(
      signif99 == T ~ "***",
      signif95 == T ~ "**",
      signif90 == T ~ "*",
      T ~ ""
    ))

  tabla_est
}


# 2. ANALYSIS: LATE with 3 treatments (A, B, C) -------------

## Preparing table of estimates
tidy.rdrobust <- function(model, ...) {
  ret <- data.frame(
    term = row.names(model$coef),
    estimate = model$coef[, 1], # Conventional (nor bias-corrected, nor robust-bias-corrected)
    std.error = model$se[, 1],
    p.value = model$pv[3] # Robust-bias-corrected inference
  )
  row.names(ret) <- NULL
  ret
}

glance.rdrobust <- function(model, ...) {
  ret <- data.frame(
    Bandwidth = round(model$bws[1], 3),
    N_h = sum(model$N_h),
    Kernel = model$kernel,
    Polynomial_degree = round(model$p)
  )
  ret
}

cm <- c("Conventional" = "Treatment")


## 2.1. Treatment 1 (D1: B vs. A) ------------
### Quantitative outcomes about employment ----------------

### Outcome 1.1: days worked post months [1, 6]
rddD1ei6_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL, all = T
)
summary(rddD1ei6_p1kT)

rddD1ei6_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(rddD1ei6_p2kT)

# Slight change at the optimal bw
srdd_2_D1ei6_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", h = rddD1ei6_p1kT$bws[1, 1] - 0.01,
  c = 0, p = 1, cluster = NULL
)
summary(srdd_2_D1ei6_p1kT)

srdd_2_D1ei6_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", h = rddD1ei6_p2kT$bws[1, 1] - 0.01,
  c = 0, p = 2, cluster = NULL
)
summary(srdd_2_D1ei6_p2kT)

# Table E3
if (!file.exists("intermediate/script02/E3_cont_D1_post6_nocl.docx")) {

  modelsummary(list(rddD1ei6_p1kT, rddD1ei6_p2kT, srdd_2_D1ei6_p1kT, srdd_2_D1ei6_p2kT),
               statistic = "std.error", coef_map = cm,
               stars = c('*' = .1, '**' = .05, '***' = 0.01),
               output = "intermediate/script02/E3_cont_D1_post6_nocl.docx")
}



### Outcome 1.2: days worked post months [7, 12]
rddD1ei12_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval712, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
  cluster = NULL
)
summary(rddD1ei12_p1kT)

rddD1ei12_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval712, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(rddD1ei12_p2kT)

# Slight change at the bw
srddD1ei12_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval712, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", c = 0, p = 1,
  h = rddD1ei12_p1kT$bws[1, 1] - 0.01,
  cluster = NULL
)
summary(srddD1ei12_p1kT)

srddD1ei12_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval712, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, cluster = NULL,
  h = rddD1ei12_p2kT$bws[1, 1] - 0.01
)
summary(srddD1ei12_p2kT)

# Table E4
if (!file.exists("intermediate/script02/E4_cont_D1_post12_nocl.docx")) {

modelsummary(list(rddD1ei12_p1kT, rddD1ei12_p2kT, srddD1ei12_p1kT, srddD1ei12_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E4_cont_D1_post12_nocl.docx")

}

### Outcome 1.3: days worked post [13, 18] months
rddD1ei18_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval1318, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
  cluster = NULL
)
summary(rddD1ei18_p1kT)

rddD1ei18_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval1318, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", c = 0, p = 2, bwselect = "mserd",
  cluster = NULL
)
summary(rddD1ei18_p2kT)

# Slight change at the bw
srddD1ei18_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval1318, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", c = 0, p = 1,
  h = rddD1ei18_p1kT$bws[1, 1] - 0.01,
  cluster = NULL
)
summary(srddD1ei18_p1kT)

srddD1ei18_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval1318, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", c = 0, p = 2,
  h = rddD1ei18_p2kT$bws[1, 1] - 0.01,
  cluster = NULL
)
summary(srddD1ei18_p2kT)

## Table E5
if (!file.exists("intermediate/script02/E5_cont_D1_post18_nocl.docx")) {

modelsummary(list(rddD1ei18_p1kT, rddD1ei18_p2kT, srddD1ei18_p1kT, srddD1ei18_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E5_cont_D1_post18_nocl.docx")

}

### Outcome 1.4: days worked post [19, 24] months
rdd24_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval1924, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(rdd24_p1kT)

rdd24_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval1924, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(rdd24_p2kT)

# Slight change at the optimal bw
srdd24_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval1924, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, cluster = NULL,
  h = rdd24_p1kT$bws[1, 1] - 0.01
)
summary(srdd24_p1kT)

srdd24_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval1924, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, cluster = NULL,
  h = rdd24_p2kT$bws[1, 1] - 0.01
)
summary(srdd24_p2kT)

# Table E6
if (!file.exists("intermediate/script02/E6_cont_D1_post24_nocl.docx")) {

modelsummary(list(rdd24_p1kT, rdd24_p2kT, srdd24_p1kT, srdd24_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E6_cont_D1_post24_nocl.docx")

}


### Qualitative outcomes about employment -------------------

### Outcome 2: longest type of contract, i.e. declared duration [with Xu's method]
#### Joining datasets
longest$id <- paste0(longest$ide_soggetto, "_", longest$axl_data_attribuzione)
longestforjoin <- longest[, c("id", "longestcontract", "longestcontract_standard")]

indi_ns_ss1$id <- paste0(indi_ns_ss1$ide_soggetto, "_", indi_ns_ss1$axl_data_attribuzione)
indi_ns_ss1_c <- left_join(indi_ns_ss1, longestforjoin, by = "id")

# Italian classification of contracts

indi_ns_ss1_c$lc_AC <- ifelse(indi_ns_ss1_c$longestcontract == "AC", 1, 0)
indi_ns_ss1_c$lc_DC <- ifelse(indi_ns_ss1_c$longestcontract == "DC", 1, 0)
indi_ns_ss1_c$lc_FTC <- ifelse(indi_ns_ss1_c$longestcontract == "FTC", 1, 0)
indi_ns_ss1_c$lc_NE <- ifelse(indi_ns_ss1_c$longestcontract == "NE", 1, 0)
indi_ns_ss1_c$lc_OEC <- ifelse(indi_ns_ss1_c$longestcontract == "OEC", 1, 0)
indi_ns_ss1_c$lc_SC <- ifelse(indi_ns_ss1_c$longestcontract == "SC", 1, 0)

# Standard classification of contracts
indi_ns_ss1_c$lcris_OEC <- ifelse(indi_ns_ss1_c$longestcontract_standard == "OEC", 1, 0)
indi_ns_ss1_c$lcris_FTCm12 <- ifelse(indi_ns_ss1_c$longestcontract_standard == "FTCm12", 1, 0)
indi_ns_ss1_c$lcris_FTC612 <- ifelse(indi_ns_ss1_c$longestcontract_standard == "FTC612", 1, 0)
indi_ns_ss1_c$lcris_FTCme6 <- ifelse(indi_ns_ss1_c$longestcontract_standard == "FTCme6", 1, 0)

# Multinomial

## Order levels of the factor

indi_ns_ss1_c$longestcontract_standardo <- factor(indi_ns_ss1_c$longestcontract_standard,
                                                  levels = c("OEC", "FTCm12", "FTC612", "FTCme6", "NE"))

homod1quali <- rdcate_multinom(score = indi_ns_ss1_c$scoringD1_0,
                               outcome = indi_ns_ss1_c$longestcontract_standardo)

homod1quali

### Outcome 3: awarded contract, i.e. duration checked at month 3
indi_ns_ss1_c$awardedcontract <- factor(indi_ns_ss1_c$prizetype ,
                                                  levels = c("prize1", "prize2", "prize3", "prize0")) # the last category is the base category

awarded_d1quali <- rdcate_multinom(score = indi_ns_ss1_c$scoringD1_0,
                               outcome = indi_ns_ss1_c$awardedcontract)

awarded_d1quali

### Saving indi_ns_ss1_c
#saveRDS(indi_ns_ss1_c, "intermediate/script02/indi_ns_ss1_c.RDS")


### Outcome 4: occupation of the longest employment spell after treatment
indi_ns_ss1_c1 <- indi_ns_ss1_c %>%
  mutate(ide_data = paste(ide_soggetto, axl_data_attribuzione, sep = "_")) %>%
  left_join(., occupations, by = "ide_data")

indi_ns_ss1_c1 <- indi_ns_ss1_c1 %>%
  mutate(occ2 = str_sub(occ_longeste, 1, 3)) %>%
  mutate(occ2 = if_else(is.na(occ2), "Any", occ2))

indi_ns_ss1_c1$occ2 <- factor(indi_ns_ss1_c1$occ2) # the last category is the base category

occ2_d1 <- rdcate_multinom(score = indi_ns_ss1_c1$scoringD1_0,
                                   outcome = indi_ns_ss1_c1$occ2)

View(occ2_d1)




### Mechanisms: hours of treatment ----------
### Outcome HS: actual treatment received of job search measures
actualtr_p1kT <- rdrobust(
  y = indi_ns_ss1$jshours, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(actualtr_p1kT)

actualtr_p2kT <- rdrobust(
  y = indi_ns_ss1$jshours, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(actualtr_p2kT)

## With slight change at the bandwidth
sactualtr_p1kT <- rdrobust(
  y = indi_ns_ss1$jshours, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, cluster = NULL,
  h = actualtr_p1kT$bws[1, 1] - 0.01
)
summary(sactualtr_p1kT)

sactualtr_p2kT <- rdrobust(
  y = indi_ns_ss1$jshours, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, cluster = NULL,
  h = actualtr_p2kT$bws[1, 1] - 0.01
)
summary(sactualtr_p2kT)

# Table E11

if (!file.exists("intermediate/script02/E11_cont_D1_jshours_nocl.docx")) {

modelsummary(list(actualtr_p1kT, actualtr_p2kT, sactualtr_p1kT, sactualtr_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E11_cont_D1_jshours_nocl.docx")

}

### Outcome HT: actual treatment received of training measures
prevtr_p1kT <- rdrobust(
  y = indi_ns_ss1$attiv_form_ore_prev, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(prevtr_p1kT)

prevtr_p2kT <- rdrobust(
  y = indi_ns_ss1$attiv_form_ore_prev, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(prevtr_p2kT)

# Slight changes at the optimal bw
sprevtr_p1kT <- rdrobust(
  y = indi_ns_ss1$attiv_form_ore_prev, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, cluster = NULL,
  h = prevtr_p1kT$bws[1, 1] - 0.01
)
summary(sprevtr_p1kT)

sprevtr_p2kT <- rdrobust(
  y = indi_ns_ss1$attiv_form_ore_prev, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, cluster = NULL,
  h = prevtr_p2kT$bws[1, 1] - 0.01
)
summary(sprevtr_p2kT)

# Table E12

if (!file.exists("intermediate/script02/E12_cont_D1_trhours_nocl.docx")) {

modelsummary(list(prevtr_p1kT, prevtr_p2kT, sprevtr_p1kT, sprevtr_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E12_cont_D1_trhours_nocl.docx")

}

## 2.2. Treatment 2 (D2: C vs. B) ----------

### Quantitative outcomes on employment -----------------

### Outcome 1.1: days worked post months [1, 6]
rddei6_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval6, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(rddei6_p1kT)

rddei6_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval6, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(rddei6_p2kT)

# Slight change at the optimal bw
srddei6_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval6, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, cluster = NULL,
  h = rddei6_p1kT$bws[1, 1] - 0.01
)
summary(srddei6_p1kT)

srddei6_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval6, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, cluster = NULL,
  h = rddei6_p2kT$bws[1, 1] - 0.01
)
summary(srddei6_p2kT)

# Table E7
if (!file.exists("intermediate/script02/E7_cont_D2_post6_nocl.docx")) {

modelsummary(list(rddei6_p1kT, rddei6_p2kT, srddei6_p1kT, srddei6_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E7_cont_D2_post6_nocl.docx")
}

### Outcome 1.2: days worked post [7, 12] months
rddei12_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval712, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(rddei12_p1kT)

rddei12_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval712, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(rddei12_p2kT)

# Slight changes at the optimal bw
srddei12_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval712, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, cluster = NULL,
  h = rddei12_p1kT$bws[1, 1] - 0.01
)
summary(srddei12_p1kT)

srddei12_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval712, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, cluster = NULL,
  h = rddei12_p2kT$bws[1, 1] - 0.01
)
summary(srddei12_p2kT)

# Table E8
if (!file.exists("intermediate/script02/E8_cont_D2_post712_nocl.docx")) {

modelsummary(list(rddei12_p1kT, rddei12_p2kT, srddei12_p1kT, srddei12_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E8_cont_D2_post712_nocl.docx")

}

### Outcome 1.3: days worked post [13, 18] months
rddei18_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval1318, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(rddei18_p1kT)

rddei18_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval1318, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(rddei18_p2kT)

# Slight change at the bw
srddei18_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval1318, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, cluster = NULL,
  h = rddei18_p1kT$bws[1, 1] - 0.01
)
summary(srddei18_p1kT)

srddei18_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval1318, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, cluster = NULL,
  h = rddei18_p2kT$bws[1, 1] - 0.01
)
summary(srddei18_p2kT)

# Table E9
if (!file.exists("intermediate/script02/E9_cont_D2_post1318_nocl.docx")) {

modelsummary(list(rddei18_p1kT, rddei18_p2kT, srddei18_p1kT, srddei18_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E9_cont_D2_post1318_nocl.docx")

}

### Outcome 1.4: days worked post [19, 24] months (RDD)
rddei24_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval1924, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(rddei24_p1kT)

rddei24_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval1924, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(rddei24_p2kT)

# Slight change at the optimal bw
srdd24_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval1924, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, cluster = NULL,
  h = rdd24_p1kT$bws[1, 1] - 0.01
)
summary(srdd24_p1kT)

srdd24_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval1924, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, cluster = NULL,
  h = rdd24_p2kT$bws[1, 1] - 0.01
)
summary(srdd24_p2kT)

cm <- c("Conventional" = "Treatment")

# Table E10
if (!file.exists("intermediate/script02/E10_cont_D2_post1924_nocl.docx")) {

modelsummary(list(rddei24_p1kT, rddei24_p2kT, srdd24_p1kT, srdd24_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E10_cont_D2_post1924_nocl.docx")

}

### Qualitative outcomes on employment ----------------
### Outcome 5: type of contract (longest one) with Xu's method
#### Joining datasets
indi_ns_ss2$id <- paste0(indi_ns_ss2$ide_soggetto, "_", indi_ns_ss2$axl_data_attribuzione)
indi_ns_ss2_c <- left_join(indi_ns_ss2, longestforjoin, by = "id")

# Italian classification of contracts
indi_ns_ss2_c$lc_AC <- ifelse(indi_ns_ss2_c$longestcontract == "AC", 1, 0)
indi_ns_ss2_c$lc_DC <- ifelse(indi_ns_ss2_c$longestcontract == "DC", 1, 0)
indi_ns_ss2_c$lc_FTC <- ifelse(indi_ns_ss2_c$longestcontract == "FTC", 1, 0)
indi_ns_ss2_c$lc_NE <- ifelse(indi_ns_ss2_c$longestcontract == "NE", 1, 0)
indi_ns_ss2_c$lc_OEC <- ifelse(indi_ns_ss2_c$longestcontract == "OEC", 1, 0)
indi_ns_ss2_c$lc_SC <- ifelse(indi_ns_ss2_c$longestcontract == "SC", 1, 0)


# Standard classification of contracts
indi_ns_ss2_c$lcris_OEC <- ifelse(indi_ns_ss2_c$longestcontract_standard == "OEC", 1, 0)
indi_ns_ss2_c$lcris_FTCm12 <- ifelse(indi_ns_ss2_c$longestcontract_standard == "FTCm12", 1, 0)
indi_ns_ss2_c$lcris_FTC612 <- ifelse(indi_ns_ss2_c$longestcontract_standard == "FTC612", 1, 0)
indi_ns_ss2_c$lcris_FTCme6 <- ifelse(indi_ns_ss2_c$longestcontract_standard == "FTCme6", 1, 0)

# Multinomial (automatic)

indi_ns_ss2_c$longestcontract_standardo <- factor(indi_ns_ss2_c$longestcontract_standard,
                                                  levels = c("OEC", "FTCm12", "FTC612", "FTCme6", "NE"))

homoD2quali <- rdcate_multinom(score = indi_ns_ss2_c$scoringD2_0,
                               outcome = indi_ns_ss2_c$longestcontract_standardo)

T5_latej <- data.frame(treat = c("D1", "", "D2", ""),
                       Y_OEC = c(paste0(homod1quali$pointest[1], homod1quali$stars[1]),
                                 homod1quali$ci95[1],
                                 paste0(homoD2quali$pointest[1], homoD2quali$stars[1]),
                                 homoD2quali$ci95[1]),
                       Y_m12 = c(paste0(homod1quali$pointest[2], homod1quali$stars[2]),
                                 homod1quali$ci95[2],
                                 paste0(homoD2quali$pointest[2], homoD2quali$stars[2]),
                                 homoD2quali$ci95[2]),
                       Y_612 = c(paste0(homod1quali$pointest[3], homod1quali$stars[3]),
                                 homod1quali$ci95[3],
                                 paste0(homoD2quali$pointest[3], homoD2quali$stars[3]),
                                 homoD2quali$ci95[3]),
                       Y_05 = c(paste0(homod1quali$pointest[4], homod1quali$stars[4]),
                                homod1quali$ci95[4],
                                paste0(homoD2quali$pointest[4], homoD2quali$stars[4]),
                                homoD2quali$ci95[4]))

T5_latej

if (!file.exists("intermediate/script02/T5_latej.xlsx")) {

  openxlsx::write.xlsx(T5_latej, 'intermediate/script02/T5_latej.xlsx')

}

### Outcome 6:
indi_ns_ss2_c$awardedcontract <- factor(indi_ns_ss2_c$prizetype ,
                                      levels = c("prize1", "prize2", "prize3", "prize0"))

awarded_d2quali <- rdcate_multinom(score = indi_ns_ss2_c$scoringD2_0,
                                   outcome = indi_ns_ss2_c$awardedcontract)

awarded_d2quali

T6_latej <- data.frame(treat = c("D1", "", "D2", ""),
                       P_OEC = c(paste0(awarded_d1quali$pointest[3], awarded_d1quali$stars[3]),
                                 awarded_d1quali$ci95[3],
                                 paste0(awarded_d2quali$pointest[3], awarded_d2quali$stars[3]),
                                 awarded_d2quali$ci95[3]),
                       P_m12 = c(paste0(awarded_d1quali$pointest[2], awarded_d1quali$stars[2]),
                                 awarded_d1quali$ci95[2],
                                 paste0(awarded_d2quali$pointest[2], awarded_d2quali$stars[2]),
                                 awarded_d2quali$ci95[2]),
                       P_612 = c(paste0(awarded_d1quali$pointest[1], awarded_d1quali$stars[1]),
                                 awarded_d1quali$ci95[1],
                                 paste0(awarded_d2quali$pointest[1], awarded_d2quali$stars[1]),
                                 awarded_d2quali$ci95[1]))

T6_latej


if (!file.exists("intermediate/script02/T6_latej.xlsx")) {

  openxlsx::write.xlsx(T6_latej, 'intermediate/script02/T6_latej.xlsx')

}

### Saving indi_ns_ss2_c
#saveRDS(indi_ns_ss2_c, "intermediate/script02/indi_ns_ss2_c.RDS")

### Mechanisms: hours of treatment -----------------
### Outcome HS: actual treatment received of job search measures
actualtr_p1kT_d2 <- rdrobust(
  y = indi_ns_ss2$jshours, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(actualtr_p1kT_d2)

actualtr_p2kT_d2 <- rdrobust(
  y = indi_ns_ss2$jshours, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(actualtr_p2kT_d2)

# Slight changes at the optimal bw
sactualtr_p1kT <- rdrobust(
  y = indi_ns_ss2$jshours, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, cluster = NULL,
  h = actualtr_p1kT$bws[1, 1] - 0.01
)
summary(sactualtr_p1kT)

sactualtr_p2kT <- rdrobust(
  y = indi_ns_ss2$jshours, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, cluster = NULL,
  h = actualtr_p2kT$bws[1, 1] - 0.01
)
summary(sactualtr_p2kT)


# Table E13
if (!file.exists("intermediate/script02/E13_cont_D2_jshours_nocl.docx")) {

modelsummary(list(actualtr_p1kT_d2, actualtr_p2kT_d2, sactualtr_p1kT, sactualtr_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E13_cont_D2_jshours_nocl.docx")

}

#### Outcome HT: actual treatment received of training measures
prevtr_p1kT_d2 <- rdrobust(
  y = indi_ns_ss2$attiv_form_ore_prev, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(prevtr_p1kT_d2)

prevtr_p2kT_d2 <- rdrobust(
  y = indi_ns_ss2$attiv_form_ore_prev, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(prevtr_p2kT_d2)

# Slight changes at the optimal bw
sprevtr_p1kT <- rdrobust(
  y = indi_ns_ss2$attiv_form_ore_prev, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, cluster = NULL,
  h = prevtr_p1kT$bws[1, 1] - 0.01
)
summary(sprevtr_p1kT)

sprevtr_p2kT <- rdrobust(
  y = indi_ns_ss2$attiv_form_ore_prev, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, cluster = NULL,
  h = prevtr_p2kT$bws[1, 1] - 0.01
)
summary(sprevtr_p2kT)

# Table E14
if (!file.exists("intermediate/script02/E14_cont_D2_trhours_nocl.docx")) {

modelsummary(list(prevtr_p1kT_d2, prevtr_p2kT_d2, sprevtr_p1kT, sprevtr_p2kT),
             statistic = "std.error", coef_map = cm,
             stars = c('*' = .1, '**' = .05, '***' = 0.01),
             output = "intermediate/script02/E14_cont_D2_trhours_nocl.docx")

}

## 2.3. Writing the Stata script to estimate AMSE ---------------

amse_script_d1 <- glue("net install rdmse, from(https://raw.githubusercontent.com/peizhuan/rdmse/master) replace
                       use \"C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script01/indi_ns_stata_D1.dta\"

                      matrix res = J(1, 12, .) // creates the matrix to save results


                       * Y1.1
                       rdmse post_interval6 scoringD1_0, deriv(0) c(0) p(1) h({rddD1ei6_p1kT$bws[1,1]}) b({rddD1ei6_p1kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,1] = r(amse_cl)

                       rdmse post_interval6 scoringD1_0, deriv(0) c(0) p(2) h({rddD1ei6_p2kT$bws[1,1]}) b({rddD1ei6_p2kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,2] = r(amse_cl)

                       * Y1.2
                       rdmse post_interval712 scoringD1_0, deriv(0) c(0) p(1) h({rddD1ei12_p1kT$bws[1,1]}) b({rddD1ei12_p1kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,3] = r(amse_cl)

                       rdmse post_interval712 scoringD1_0, deriv(0) c(0) p(2) h({rddD1ei12_p2kT$bws[1,1]}) b({rddD1ei12_p2kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,4] = r(amse_cl)

                       * Y1.3
                       rdmse post_interval1318 scoringD1_0, deriv(0) c(0) p(1) h({rddD1ei18_p1kT$bws[1,1]}) b({rddD1ei18_p1kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,5] = r(amse_cl)

                       rdmse post_interval1318 scoringD1_0, deriv(0) c(0) p(2) h({rddD1ei18_p2kT$bws[1,1]}) b({rddD1ei18_p2kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,6] = r(amse_cl)

                       * Y1.4
                       rdmse post_interval1924 scoringD1_0, deriv(0) c(0) p(1) h({rdd24_p1kT$bws[1,1]}) b({rdd24_p1kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,7] = r(amse_cl)

                       rdmse post_interval1924 scoringD1_0, deriv(0) c(0) p(2) h({rdd24_p2kT$bws[1,1]}) b({rdd24_p2kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,8] = r(amse_cl)

                       * MECHANISMS (js hours)
                       rdmse jshours scoringD1_0, deriv(0) c(0) p(1) h({actualtr_p1kT$bws[1,1]}) b({actualtr_p1kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,9] = r(amse_cl)

                       rdmse jshours scoringD1_0, deriv(0) c(0) p(2) h({actualtr_p2kT$bws[1,1]}) b({actualtr_p2kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,10] = r(amse_cl)

                       * MECHANISMS (tr hours)
                       rdmse attiv_form_ore_prev scoringD1_0, deriv(0) c(0) p(1) h({prevtr_p1kT$bws[1,1]}) b({prevtr_p1kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,11] = r(amse_cl)

                       rdmse attiv_form_ore_prev scoringD1_0, deriv(0) c(0) p(2) h({prevtr_p2kT$bws[1,1]}) b({prevtr_p2kT$bws[2,1]}) kernel(triangular)
                       matrix res[1,12] = r(amse_cl)

                       * Export to Excel
                       // Crear un dataset temporal con solo la matriz
                       clear // Asegurarse de que no haya datos existentes
                       svmat res, names(col)

                       // Exportar el dataset como Excel
                       export excel using \"C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script02/stata/resu_amse_1.xlsx\", firstrow(variables) replace

                       ")

writeLines(amse_script_d1, "intermediate/script02/amse_script_d1.txt")



amse_script_d2 <- glue("net install rdmse, from(https://raw.githubusercontent.com/peizhuan/rdmse/master) replace
                       use \"C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script01/indi_ns_stata_D2.dta\"

                       matrix resu = J(1, 12, .) // creates the matrix to save results


                       * Y1.1
                       rdmse post_interval6 scoringD2_0, deriv(0) c(0) p(1) h({rddei6_p1kT$bws[1,1]}) b({rddei6_p1kT$bws[2,1]}) kernel(triangular)
                       matrix resu[1,1] = r(amse_cl)

                       rdmse post_interval6 scoringD2_0, deriv(0) c(0) p(2) h({rddei6_p2kT$bws[1,1]}) b({rddei6_p2kT$bws[2,1]}) kernel(triangular)
                       matrix resu[1,2] = r(amse_cl)

                       * Y1.2
                       rdmse post_interval712 scoringD2_0, deriv(0) c(0) p(1) h({rddei12_p1kT$bws[1,1]}) b({rddei12_p1kT$bws[2,1]}) kernel(triangular)
                       matrix resu[1,3] = r(amse_cl)

                       rdmse post_interval712 scoringD2_0, deriv(0) c(0) p(2) h({rddei12_p2kT$bws[1,1]}) b({rddei12_p2kT$bws[2,1]}) kernel(triangular)
                       matrix resu[1,4] = r(amse_cl)

                       * Y1.3
                       rdmse post_interval1318 scoringD2_0, deriv(0) c(0) p(1) h({rddei18_p1kT$bws[1,1]}) b({rddei18_p1kT$bws[2,1]}) kernel(triangular)
                       matrix resu[1,5] = r(amse_cl)

                       rdmse post_interval1318 scoringD2_0, deriv(0) c(0) p(2) h({rddei18_p2kT$bws[1,1]}) b({rddei18_p2kT$bws[2,1]}) kernel(triangular)
                       matrix resu[1,6] = r(amse_cl)

                       * Y1.4
                       rdmse post_interval1924 scoringD2_0, deriv(0) c(0) p(1) h({rddei24_p1kT$bws[1,1]}) b({rddei24_p1kT$bws[2,1]}) kernel(triangular)
                       matrix resu[1,7] = r(amse_cl)

                       rdmse post_interval1924 scoringD2_0, deriv(0) c(0) p(2) h({rddei24_p2kT$bws[1,1]}) b({rddei24_p2kT$bws[2,1]}) kernel(triangular)
                       matrix resu[1,8] = r(amse_cl)

                       * MECHANISMS (js hours)
                       rdmse jshours scoringD2_0, deriv(0) c(0) p(1) h({actualtr_p1kT_d2$bws[1,1]}) b({actualtr_p1kT_d2$bws[2,1]}) kernel(triangular)
                       matrix resu[1,9] = r(amse_cl)

                       rdmse jshours scoringD2_0, deriv(0) c(0) p(2) h({actualtr_p2kT_d2$bws[1,1]}) b({actualtr_p2kT_d2$bws[2,1]}) kernel(triangular)
                       matrix resu[1,10] = r(amse_cl)

                       * MECHANISMS (tr hours)
                       rdmse attiv_form_ore_prev scoringD2_0, deriv(0) c(0) p(1) h({prevtr_p1kT_d2$bws[1,1]}) b({prevtr_p1kT_d2$bws[2,1]}) kernel(triangular)
                       matrix resu[1,11] = r(amse_cl)

                       rdmse attiv_form_ore_prev scoringD2_0, deriv(0) c(0) p(2) h({prevtr_p2kT_d2$bws[1,1]}) b({prevtr_p2kT_d2$bws[2,1]}) kernel(triangular)
                       matrix resu[1,12] = r(amse_cl)

                       * Export to Excel
                       // Crear un dataset temporal con solo la matriz
                       clear // Asegurarse de que no haya datos existentes
                       svmat resu, names(col)

                       // Exportar el dataset como Excel
                       export excel using \"C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script02/stata/resu_amse_2.xlsx\", firstrow(variables) replace

                       ")

writeLines(amse_script_d2, "intermediate/script02/amse_script_d2.txt")


# 3. ANALYSIS: LQTE ------------
indi_ns_ss1$treatedD1 <- ifelse(indi_ns_ss1$scoringD1_0 > 0, TRUE, FALSE)
indi_ns_ss2$treatedD2 <- ifelse(indi_ns_ss2$scoringD2_0 > 0, TRUE, FALSE)

running1 <- indi_ns_ss1$scoringD1_0

outcome0 <- indi_ns_ss1$unemplength
outcome1 <- indi_ns_ss1$post_interval6
outcome2 <- indi_ns_ss1$post_interval712
outcome3 <- indi_ns_ss1$post_interval1318
outcome4 <- indi_ns_ss1$post_interval1924
d1 <- indi_ns_ss1$treatedD1
x0 <- 0

running2 <- indi_ns_ss2$scoringD2_0
outcome0_d2 <- indi_ns_ss2$unemplength
outcome1_d2 <- indi_ns_ss2$post_interval6
outcome2_d2 <- indi_ns_ss2$post_interval712
outcome3_d2 <- indi_ns_ss2$post_interval1318
outcome4_d2 <- indi_ns_ss2$post_interval1924
d2 <- indi_ns_ss2$treatedD2



quantile(indi_ns_ss2$post_interval6[which(indi_ns_ss2$scoringD2_0 > -0.1 & indi_ns_ss2$scoringD2_0 < 0.1)],
         prob=seq(0, 1, length = 21)) # para obtener deciles



container_lqte_d1 <- data.frame(object = c("q1_pe", "q2_pe", "q3_pe", "nonsign_test_stars", "homo_test_stars"),
                                y1 = rep(NA, 5),
                                y2 = rep(NA, 5),
                                y3 = rep(NA, 5),
                                y4 = rep(NA, 5)) # pe = point estimate

container_lqte_d2 <- data.frame(object = c("q1_pe", "q2_pe", "q3_pe", "nonsign_test_stars", "homo_test_stars"),
                                y1 = rep(NA, 5),
                                y2 = rep(NA, 5),
                                y3 = rep(NA, 5),
                                y4 = rep(NA, 5)) # pe = point estimate

# Estimation and inference following Qu and Yoon (2019)
#source("inputcode/2018quyoon_code/qte_rdd.R")
## Following the same logic as Cattaneo et al. (2020, p.65), the point estimate obtained from a MSE-optimal bw is not bias-corrected.


## 3.1.  Treatment 1 (D1: B vs. A) ----------

### Outcome 1.1: days worked post months [1, 6] -------------
# 1) Estimate bandwidth
bwd1o1_bdy <- rdq.bandwidth(y = outcome1, x = running1, d = d1, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 1, val = seq(0.01, 0.1, by = 0.01))


bwd1o1_int <- rdq.bandwidth(y = outcome1, x = running1, d = d1, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 0, val = seq(0.01, 0.1, by = 0.01))


bwd1o1_bdyint <- data.frame(cv = c(bwd1o1_bdy$cv, bwd1o1_int$cv),
                            bdy = c(bwd1o1_bdy$opt.m, bwd1o1_bdy$opt.p),
                            int = c(bwd1o1_int$opt.m, bwd1o1_int$opt.p))

#saveRDS(bwd1o1_bdyint, "intermediate/script02/lqte/bwd1o1_bdyint.RDS")
bwd1o1_bdyint <- readRDS("intermediate/script02/lqte/bwd1o1_bdyint.RDS")


bwd1o1_min <- min(bwd1o1_bdyint[, 2:3])

# 2) Estimate quantiles
## a. Detecting quantile effects == 0
qte_d1o1_05 <- rd.qte(y = outcome1, x = running1, d = d1, x0 = 0,
                        cov = 0, bias = 0, bdw = bwd1o1_min,
                        tau = seq(0.1, 0.9, by = 0.05))
qte_d1o1_05

qte_d1o1 <- rd.qte(y = outcome1, x = running1, d = d1, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd1o1_min,
                      tau = seq(0.25, 0.75, by = 0.25))
qte_d1o1


## b. Estimable quantile effects
qte_d1o1_05_est <- rdq.band(y = outcome1, x = running1, d = d1, x0 = 0,
                        cov = 0, bdw = bwd1o1_min, alpha = 0.1,
                        tau = seq(0.35, 0.9, by = 0.05))

#saveRDS(qte_d1o1, "intermediate/script02/lqte/qte_d1o1.RDS")
#saveRDS(qte_d1o1_05_est, "intermediate/script02/lqte/qte_d1o1_05_est.RDS")

qte_d1o1 <- readRDS("intermediate/script02/lqte/qte_d1o1.RDS")
qte_d1o1_05_est <- readRDS("intermediate/script02/lqte/qte_d1o1_05_est.RDS")

summary(qte_d1o1, alpha=0.1)
summary(qte_d1o1, alpha=0.05)
summary(qte_d1o1, alpha=0.01)

## c. Plot
df_d1o1 <- data.frame(
  cuantil = qte_d1o1_05_est$tau,
  efecto_estimado = qte_d1o1_05_est$qte,
  lower_band = qte_d1o1_05_est[["uband.robust"]][, 1, ],
  upper_band = qte_d1o1_05_est[["uband.robust"]][, 2, ])

# Gráfico
p_lqte_d1y1 <- ggplot(df_d1o1, aes(x = cuantil, y = efecto_estimado)) +
  geom_ribbon(aes(ymin = lower_band, ymax = upper_band), fill = "grey", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 0) +
  theme_light() +
  ylim(-180, +180) +
  labs(x = expression(paste("Quantile (", tau, ")")),
       y = expression(paste("Estimate of LQTE(", tau, ") of ", D[1]," on ", Y[1]))) # Saved at 395x300

ggsave(filename = "intermediate/script02/plots/lqte_d1_y1.svg",  plot = p_lqte_d1y1,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


# 3) Run hypothesis tests

# Running the following chunk with rdq.test we find (using traceback()) the following trace of the error:

# Error en solve.default(crossprod((kx * xj), (d * fxp[, k] * xj))/en):
#  Lapack routine dgesv: system is exactly singular: U[1,1] = 0

# 4: solve.default(crossprod((kx * xj), (d * fxp[, k] * xj))/en)
# 3: solve(crossprod((kx * xj), (d * fxp[, k] * xj))/en)
# 2: rdq.sim(x, d, x0, z0, dz, cov, tt = tt, hh2[ind], hh2[ind], fxp = fp$ff,
#           fxm = fm$ff, n.sim) at #81
# 1: rdq.test2(y = outcome1, x = running1, d = d1, x0 = 0, cov = 0,
#             alpha = 0.1, tau = seq(0.1, 0.9, by = 0.05), bdw = bwd1o1_min,
#             bias = 1, type = c(1, 2))

# Denotemos Ñ a la matriz fruto de crossprod((kx * xj), (d * fxp[, k] * xj))
# No podemos resolver Ñ/en porque Ñ contiene vectores cero!

Wd1o1_signhomo <- rdq.test2(y = outcome1, x = running1, d = d1, x0 = 0,
                           cov = 0,
                           alpha = 0.1, tau = seq(0.1, 0.9, by = 0.05),
                           bdw = bwd1o1_min, bias = 1, type = c(1,2)) # to debug, check ?rdq.sim

fxp_prob <- as.data.frame(Wd1o1_signhomo)
all(fxp_prob$V1 == fxp_prob$V2)
all(fxp_prob$V2 == fxp_prob$V3)


Wd1o1_signhomo_nz <- rdq.test(y = outcome1, x = running1, d = d1, x0 = 0,
                           cov = 0,
                           alpha = 0.95, tau = seq(0.35, 0.9, by = 0.05),
                           bdw = bwd1o1_min, bias = 1, type = c(1,2))

#saveRDS(Wd1o1_signhomo, "intermediate/script02/lqte/Wd1o1_signhomo.RDS")


# 4) Writing the table
container_lqte_d1$y1[1:3] <- round(qte_d1o1$qte, 3)

container_lqte_d1$y1[4] <- if(Wd1o1_signhomo_nz$p.value$significance[1] < 0.01) {"***"} else{
  if(Wd1o1_signhomo_nz$p.value$significance[1] < 0.05) {"**"} else{
    if(Wd1o1_signhomo_nz$p.value$significance[1] < 0.1) {"*"} else{""}
  }
}

container_lqte_d1$y1[5] <- if(Wd1o1_signhomo_nz$p.value$homogeneity[1] < 0.01) {"***"} else{
  if(Wd1o1_signhomo_nz$p.value$homogeneity[1] < 0.05) {"**"} else{
    if(Wd1o1_signhomo_nz$p.value$homogeneity[1] < 0.1) {"*"} else{""}
  }
}






### Outcome 1.2: days worked post months [7, 12] ---------------
# 1) Estimate bandwidth
bwd1o2_bdy <- rdq.bandwidth(y = outcome2, x = running1, d = d1, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 1, val = seq(0.01, 0.1, by = 0.01))


bwd1o2_int <- rdq.bandwidth(y = outcome2, x = running1, d = d1, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 0, val = seq(0.01, 0.1, by = 0.01))


bwd1o2_bdyint <- data.frame(cv = c(bwd1o2_bdy$cv, bwd1o2_int$cv),
                            bdy = c(bwd1o2_bdy$opt.m, bwd1o2_bdy$opt.p),
                            int = c(bwd1o2_int$opt.m, bwd1o2_int$opt.p))

#saveRDS(bwd1o2_bdyint, "intermediate/script02/lqte/bwd1o2_bdyint.RDS")
bwd1o2_bdyint <- readRDS("intermediate/script02/lqte/bwd1o2_bdyint.RDS")

bwd1o2_min <- min(bwd1o2_bdyint)

# 2) Estimate quantiles
## a. Detecting quantile effects == 0
qte_d1o2_05 <- rd.qte(y = outcome2, x = running1, d = d1, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd1o2_min,
                   tau = seq(0.1, 0.9, by = 0.05))

qte_d1o2_05

qte_d1o2 <- rd.qte(y = outcome2, x = running1, d = d1, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd1o2_min,
                   tau = seq(0.25, 0.75, by = 0.25))


saveRDS(qte_d1o2, "intermediate/script02/lqte/qte_d1o2.RDS")

## b. Non-null quantile effects
qte_d1o2_05_est <- rdq.band(y = outcome2, x = running1, d = d1, x0 = 0,
                            cov = 0, bdw = bwd1o2_min, alpha = 0.1,
                            tau = seq(0.3, 0.75, by = 0.05))

saveRDS(qte_d1o2_05_est, "intermediate/script02/lqte/qte_d1o2_05_est.RDS")

## c. Plot
df_d1o2 <- data.frame(
  cuantil = qte_d1o2_05_est$tau,
  efecto_estimado = qte_d1o2_05_est$qte,
  lower_band = qte_d1o2_05_est[["uband.robust"]][, 1, ],
  upper_band = qte_d1o2_05_est[["uband.robust"]][, 2, ])

# Gráfico
p_lqte_d1y2 <- ggplot(df_d1o2, aes(x = cuantil, y = efecto_estimado)) +
  geom_ribbon(aes(ymin = lower_band, ymax = upper_band), fill = "grey", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 0) +
  theme_light() +
  ylim(-180, +180) +
  labs(x = expression(paste("Quantile (", tau, ")")),
       y = expression(paste("Estimate of LQTE(", tau, ") of ", D[1]," on ", Y[2]))) # Saved at 395x300

ggsave(filename = "intermediate/script02/plots/lqte_d1_y2.svg",  plot = p_lqte_d1y2,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


# 3) Run hypothesis tests
Wd1o2_signhomo_nz <- rdq.test(y = outcome2, x = running1, d = d1, x0 = 0,
                           cov = 0,
                           alpha = 0.95, tau = seq(0.3, 0.75, by = 0.05),
                           bdw = bwd1o2_min, bias = 1, type = c(1,2))

#saveRDS(Wd1o2_signhomo_nz, "intermediate/script02/lqte/Wd1o2_signhomo_nz.RDS")


# 4) Writing the table
container_lqte_d1$y2[1:3] <- round(qte_d1o2$qte, 3)

container_lqte_d1$y2[4] <- if(Wd1o2_signhomo_nz$p.value$significance[1] < 0.01) {"***"} else{
  if(Wd1o2_signhomo_nz$p.value$significance[1] < 0.05) {"**"} else{
    if(Wd1o2_signhomo_nz$p.value$significance[1] < 0.1) {"*"} else{""}
  }
}

container_lqte_d1$y2[5] <- if(Wd1o2_signhomo_nz$p.value$homogeneity[1] < 0.01) {"***"} else{
  if(Wd1o2_signhomo_nz$p.value$homogeneity[1] < 0.05) {"**"} else{
    if(Wd1o2_signhomo_nz$p.value$homogeneity[1] < 0.1) {"*"} else{""}
  }
}



### Outcome 1.3: days worked post months [13, 18] ------------
# 1) Estimate bandwidth
bwd1o3_bdy <- rdq.bandwidth(y = outcome3, x = running1, d = d1, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 1, val = seq(0.01, 0.1, by = 0.01))


bwd1o3_int <- rdq.bandwidth(y = outcome3, x = running1, d = d1, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 0, val = seq(0.01, 0.1, by = 0.01))


bwd1o3_bdyint <- data.frame(cv = c(bwd1o3_bdy$cv, bwd1o3_int$cv),
                            bdy = c(bwd1o3_bdy$opt.m, bwd1o3_bdy$opt.p),
                            int = c(bwd1o3_int$opt.m, bwd1o3_int$opt.p))

#saveRDS(bwd1o3_bdyint, "intermediate/script02/lqte/bwd1o3_bdyint.RDS")
bwd1o3_bdyint <- readRDS("intermediate/script02/lqte/bwd1o3_bdyint.RDS")

bwd1o3_min <- min(bwd1o3_bdyint)


# 2) Estimate quantiles
## a. Detecting quantile effects == 0
qte_d1o3_05 <- rd.qte(y = outcome3, x = running1, d = d1, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd1o3_min,
                   tau = seq(0.1, 0.9, by = 0.05))

qte_d1o3 <- rd.qte(y = outcome3, x = running1, d = d1, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd1o3_min,
                   tau = seq(0.25, 0.75, by = 0.25))

saveRDS(qte_d1o3, "intermediate/script02/lqte/qte_d1o3.RDS")

## b. Non-null quantile effects
qte_d1o3_05_est <- rdq.band(y = outcome3, x = running1, d = d1, x0 = 0,
                            cov = 0, bdw = bwd1o3_min, alpha = 0.1,
                            tau = seq(0.3, 0.75, by = 0.05))

saveRDS(qte_d1o3_05_est, "intermediate/script02/lqte/qte_d1o3_05_est.RDS")

## c. Plot
df_d1o3 <- data.frame(
  cuantil = qte_d1o3_05_est$tau,
  efecto_estimado = qte_d1o3_05_est$qte,
  lower_band = qte_d1o3_05_est[["uband.robust"]][, 1, ],
  upper_band = qte_d1o3_05_est[["uband.robust"]][, 2, ])

# Gráfico
p_lqte_d1y3 <- ggplot(df_d1o3, aes(x = cuantil, y = efecto_estimado)) +
  geom_ribbon(aes(ymin = lower_band, ymax = upper_band), fill = "grey", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 0) +
  theme_light() +
  ylim(-180, +180) +
  labs(x = expression(paste("Quantile (", tau, ")")),
       y = expression(paste("Estimate of LQTE(", tau, ") of ", D[1]," on ", Y[3]))) # Saved at 395x300

ggsave(filename = "intermediate/script02/plots/lqte_d1_y3.svg",  plot = p_lqte_d1y3,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


# 3) Run hypothesis tests
Wd1o3_signhomo_nz <- rdq.test(y = outcome3, x = running1, d = d1, x0 = 0,
                            cov = 0,
                            alpha = 0.95, tau = seq(0.3, 0.75, by = 0.05),
                            bdw = bwd1o3_min, bias = 1, type = c(1,2))



# Automatically writing the table
container_lqte_d1$y3[1:3] <- round(qte_d1o3$qte, 3)

container_lqte_d1$y3[4] <- if(Wd1o3_signhomo_nz$p.value$significance[1] < 0.01) {"***"} else{
  if(Wd1o3_signhomo_nz$p.value$significance[1] < 0.05) {"**"} else{
    if(Wd1o3_signhomo_nz$p.value$significance[1] < 0.1) {"*"} else{""}
  }
}

container_lqte_d1$y3[5] <- if(Wd1o3_signhomo_nz$p.value$homogeneity[1] < 0.01) {"***"} else{
  if(Wd1o3_signhomo_nz$p.value$homogeneity[1] < 0.05) {"**"} else{
    if(Wd1o3_signhomo_nz$p.value$homogeneity[1] < 0.1) {"*"} else{""}
  }
}



### Outcome 1.4: days worked post months [19, 24] -------------
# 1) Estimate bandwidth
bwd1o4_bdy <- rdq.bandwidth(y = outcome4, x = running1, d = d1, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 1, val = seq(0.01, 0.1, by = 0.01))

bwd1o4_bdy_min <- min(bwd1o4_bdy$opt.m, bwd1o4_bdy$opt.p)

bwd1o4_int <- rdq.bandwidth(y = outcome4, x = running1, d = d1, x0 = 0,
                        cv = 1, cov = 0, pm.each = 0,
                        bdy = 0, val = seq(0.01, 0.1, by = 0.01))

bwd1o4_int_min <- min(bwd1o4_int$opt.m, bwd1o4_int$opt.p)


bwd1o4_bdyint <- data.frame(cv = c(bwd1o4_bdy$cv, bwd1o4_int$cv),
                            bdy = c(bwd1o4_bdy$opt.m, bwd1o4_bdy$opt.p),
                            int = c(bwd1o4_int$opt.m, bwd1o4_int$opt.p))

#saveRDS(bwd1o4_bdyint, "intermediate/script02/lqte/bwd1o4_bdyint.RDS")
bwd1o4_bdyint <- readRDS("intermediate/script02/lqte/bwd1o4_bdyint.RDS")

bwd1o4_min <- min(bwd1o4_bdyint)


# 2) Estimate quantiles
## a. Detecting quantile effects == 0
qte_d1o4_05 <- rd.qte(y = outcome4, x = running1, d = d1, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd1o4_min,
                   tau = seq(0.1, 0.9, by = 0.05))

qte_d1o4 <- rd.qte(y = outcome4, x = running1, d = d1, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd1o4_min,
                   tau = seq(0.25, 0.75, by = 0.25))

#saveRDS(qte_d1o4, "intermediate/script02/lqte/qte_d1o4.RDS")
qte_d1o4 <- readRDS("intermediate/script02/lqte/qte_d1o4.RDS")

## b. Non-null quantile effects
qte_d1o4_05_est <- rdq.band(y = outcome4, x = running1, d = d1, x0 = 0,
                            cov = 0, bdw = bwd1o4_min, alpha = 0.1,
                            tau = seq(0.3, 0.75, by = 0.05))

saveRDS(qte_d1o4_05_est, "intermediate/script02/lqte/qte_d1o4_05_est.RDS")

## c. Plot
df_d1o4 <- data.frame(
  cuantil = qte_d1o4_05_est$tau,
  efecto_estimado = qte_d1o4_05_est$qte,
  lower_band = qte_d1o4_05_est[["uband.robust"]][, 1, ],
  upper_band = qte_d1o4_05_est[["uband.robust"]][, 2, ])

# Gráfico
p_lqte_d1y4 <- ggplot(df_d1o4, aes(x = cuantil, y = efecto_estimado)) +
  geom_ribbon(aes(ymin = lower_band, ymax = upper_band), fill = "grey", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 0) +
  theme_light() +
  ylim(-180, +180) +
  labs(x = expression(paste("Quantile (", tau, ")")),
       y = expression(paste("Estimate of LQTE(", tau, ") of ", D[1]," on ", Y[4]))) # Saved at 395x300

ggsave(filename = "intermediate/script02/plots/lqte_d1_y4.svg",  plot = p_lqte_d1y4,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


# 3) Run hypothesis tests
Wd1o4_signhomo_nz <- rdq.test(y = outcome4, x = running1, d = d1, x0 = 0,
                           cov = 0,
                           alpha = 0.95, tau = seq(0.3, 0.75, by = 0.05),
                           bdw = bwd1o4_min, bias = 1, type = c(1,2))

saveRDS(Wd1o4_signhomo_nz, "intermediate/script02/lqte/Wd1o4_signhomo_nz.RDS")

# Automatically writing the table
container_lqte_d1$y4[1:3] <- round(qte_d1o4$qte, 3)

container_lqte_d1$y4[4] <- if(Wd1o4_signhomo_nz$p.value$significance[1] < 0.01) {"***"} else{
  if(Wd1o4_signhomo_nz$p.value$significance[1] < 0.05) {"**"} else{
    if(Wd1o4_signhomo_nz$p.value$significance[1] < 0.1) {"*"} else{""}
  }
}

container_lqte_d1$y4[5] <- if(Wd1o4_signhomo_nz$p.value$homogeneity[1] < 0.01) {"***"} else{
  if(Wd1o4_signhomo_nz$p.value$homogeneity[1] < 0.05) {"**"} else{
    if(Wd1o4_signhomo_nz$p.value$homogeneity[1] < 0.1) {"*"} else{""}
  }
}


### Saving and exporting --------

#saveRDS(container_lqte_d1, "intermediate/script02/lqte/container_lqte_d1.RDS")
#openxlsx::write.xlsx(container_lqte_d1, 'intermediate/script02/lqte/container_lqte_d1.xlsx')

container_lqte_d1 <- readRDS("intermediate/script02/lqte/container_lqte_d1.RDS")

## 3.2. Treatment 2 (D2: C vs. B) ---------

### Outcome 1.1: days worked post months [1, 6] ------------
# 1- Estimate bandwidth
bwd2o1_bdy <- rdq.bandwidth(y = outcome1_d2, x = running2, d = d2, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 1, val = seq(0.01, 0.1, by = 0.01))


bwd2o1_int <- rdq.bandwidth(y = outcome1_d2, x = running2, d = d2, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 0, val = seq(0.01, 0.1, by = 0.01))


bwd2o1_bdyint <- data.frame(cv = c(bwd2o1_bdy$cv, bwd2o1_int$cv),
                            bdy = c(bwd2o1_bdy$opt.m, bwd2o1_bdy$opt.p),
                            int = c(bwd2o1_int$opt.m, bwd2o1_int$opt.p))

#saveRDS(bwd2o1_bdyint, "intermediate/script02/lqte/bwd2o1_bdyint.RDS")
bwd2o1_bdyint <- readRDS("intermediate/script02/lqte/bwd2o1_bdyint.RDS")


bwd2o1_min <- min(bwd2o1_bdyint[, 2:3]) # MSE-optimal


# 2- Estimate quantiles
## a. Detecting non-zero quantile effects
qte_d2o1_05 <- rd.qte(y = outcome1_d2, x = running2, d = d2, x0 = 0,
                       cov = 0, bias = 0, bdw = bwd2o1_min,
                       tau = seq(0.1, 0.9, by = 0.05))

qte_d2o1 <- rd.qte(y = outcome1_d2, x = running2, d = d2, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd2o1_min,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Estimating non-zero quantile effects
qte_d2o1_005b <- rdq.band(y = outcome1_d2, x = running2, d = d2, x0 = 0,
                       cov = 0, bdw = bwd2o1_min, alpha = 0.1,
                       tau = seq(0.45, 0.9, by = 0.05))
## c. Plot
df_d2o1 <- data.frame(
  cuantil = qte_d2o1_005b$tau,
  efecto_estimado = qte_d2o1_005b$qte,
  lower_band = qte_d2o1_005b[["uband.robust"]][, 1, ],
  upper_band = qte_d2o1_005b[["uband.robust"]][, 2, ])

p_lqte_d2y1 <- ggplot(df_d2o1, aes(x = cuantil, y = efecto_estimado)) +
  geom_ribbon(aes(ymin = lower_band, ymax = upper_band), fill = "grey", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 0) +
  theme_light() +
  ylim(-180, +180) +
  labs(x = expression(paste("Quantile (", tau, ")")),
       y = expression(paste("Estimate of LQTE(", tau, ") of ", D[2]," on ", Y[1]))) # Saved at 395x300

ggsave(filename = "intermediate/script02/plots/lqte_d2_y1.svg",  plot = p_lqte_d2y1,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# 3) Hypothesis testing
Wd2o1_signhomo_nz <- rdq.test(y = outcome1_d2, x = running2, d = d2, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.45, 0.9, by = 0.05),
                              bdw = bwd2o1_min, bias = 1, type = c(1,2))

saveRDS(Wd2o1_signhomo_nz, "intermediate/script02/lqte/Wd2o1_signhomo_nz.RDS")


# 4) Writing the table
container_lqte_d2$y1[1:3] <- round(qte_d2o1$qte, 3)

container_lqte_d2$y1[4] <- if(Wd2o1_signhomo_nz$p.value$significance[1] < 0.01) {"***"} else{
  if(Wd2o1_signhomo_nz$p.value$significance[1] < 0.05) {"**"} else{
    if(Wd2o1_signhomo_nz$p.value$significance[1] < 0.1) {"*"} else{""}
  }
}

container_lqte_d2$y1[5] <- if(Wd2o1_signhomo_nz$p.value$homogeneity[1] < 0.01) {"***"} else{
  if(Wd2o1_signhomo_nz$p.value$homogeneity[1] < 0.05) {"**"} else{
    if(Wd2o1_signhomo_nz$p.value$homogeneity[1] < 0.1) {"*"} else{""}
  }
}


### Outcome 1.2: days worked post months [7, 12] ----------------
# 1- Estimate bandwidth
bwd2o2_bdy <- rdq.bandwidth(y = outcome2_d2, x = running2, d = d2, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 1, val = seq(0.01, 0.1, by = 0.01))


bwd2o2_int <- rdq.bandwidth(y = outcome2_d2, x = running2, d = d2, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 0, val = seq(0.01, 0.1, by = 0.01))


bwd2o2_bdyint <- data.frame(cv = c(bwd2o1_bdy$cv, bwd2o2_int$cv),
                            bdy = c(bwd2o2_bdy$opt.m, bwd2o2_bdy$opt.p),
                            int = c(bwd2o2_int$opt.m, bwd2o2_int$opt.p))

#saveRDS(bwd2o2_bdyint, "intermediate/script02/lqte/bwd2o2_bdyint.RDS")
bwd2o2_bdyint <- readRDS("intermediate/script02/lqte/bwd2o2_bdyint.RDS")

bwd2o2_min <- min(bwd2o2_bdyint[, 2:3])


# 2- Estimate quantiles
## a. Detecting non-zero quantile effects
qte_d2o2_05 <- rd.qte(y = outcome2_d2, x = running2, d = d2, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd2o2_min,
                      tau = seq(0.1, 0.9, by = 0.05))

qte_d2o2 <- rd.qte(y = outcome2_d2, x = running2, d = d2, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd2o2_min,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Estimating non-zero quantile effects
quantile(indi_ns_ss2$post_interval712[which(indi_ns_ss2$scoringD2_0 > 0 & indi_ns_ss2$scoringD2_0 < bwd2o2_min)],
         prob=seq(0, 1, length = 21)) # para obtener percentiles en tratados

quantile(indi_ns_ss2$post_interval712[which(indi_ns_ss2$scoringD2_0 > -bwd2o2_min & indi_ns_ss2$scoringD2_0 < 0)],
         prob=seq(0, 1, length = 21)) # para obtener deciles en no-tratados (el percentil 40 es cero)

qte_d2o2_05b <- rdq.band(y = outcome2_d2, x = running2, d = d2, x0 = 0,
                          cov = 0, bdw = bwd2o2_min, alpha = 0.1,
                          tau = seq(0.45, 0.9, by = 0.05))

## c. Plot
df_d2o2 <- data.frame(
  cuantil = qte_d2o2_05b$tau,
  efecto_estimado = qte_d2o2_05b$qte,
  lower_band = qte_d2o2_05b[["uband.robust"]][, 1, ],
  upper_band = qte_d2o2_05b[["uband.robust"]][, 2, ])

p_lqte_d2y2 <- ggplot(df_d2o2, aes(x = cuantil, y = efecto_estimado)) +
  geom_ribbon(aes(ymin = lower_band, ymax = upper_band), fill = "grey", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 0) +
  theme_light() +
  ylim(-180, +180) +
  labs(x = expression(paste("Quantile (", tau, ")")),
       y = expression(paste("Estimate of LQTE(", tau, ") of ", D[2]," on ", Y[2]))) # Saved at 395x300

ggsave(filename = "intermediate/script02/plots/lqte_d2_y2.svg",  plot = p_lqte_d2y2,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")



# 3- Hypothesis testing
Wd2o2_signhomo_nz <- rdq.test(y = outcome2_d2, x = running2, d = d2, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.45, 0.9, by = 0.05),
                              bdw = bwd2o2_min, bias = 1, type = c(1,2))

saveRDS(Wd2o2_signhomo_nz, "intermediate/script02/lqte/Wd2o2_signhomo_nz.RDS")


# 4) Writing the table
container_lqte_d2$y2[1:3] <- round(qte_d2o2$qte, 3)

container_lqte_d2$y2[4] <- if(Wd2o2_signhomo_nz$p.value$significance[1] < 0.01) {"***"} else{
  if(Wd2o2_signhomo_nz$p.value$significance[1] < 0.05) {"**"} else{
    if(Wd2o2_signhomo_nz$p.value$significance[1] < 0.1) {"*"} else{""}
  }
}

container_lqte_d2$y2[5] <- if(Wd2o2_signhomo_nz$p.value$homogeneity[1] < 0.01) {"***"} else{
  if(Wd2o2_signhomo_nz$p.value$homogeneity[1] < 0.05) {"**"} else{
    if(Wd2o2_signhomo_nz$p.value$homogeneity[1] < 0.1) {"*"} else{""}
  }
}


### Outcome 1.3: days worked post months [13, 18] ----------------
# 1- Estimate bandwidth
bwd2o3_bdy <- rdq.bandwidth(y = outcome3_d2, x = running2, d = d2, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 1, val = seq(0.01, 0.1, by = 0.01))


bwd2o3_int <- rdq.bandwidth(y = outcome3_d2, x = running2, d = d2, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 0, val = seq(0.01, 0.1, by = 0.01))


bwd2o3_bdyint <- data.frame(cv = c(bwd2o3_bdy$cv, bwd2o3_int$cv),
                            bdy = c(bwd2o3_bdy$opt.m, bwd2o3_bdy$opt.p),
                            int = c(bwd2o3_int$opt.m, bwd2o3_int$opt.p))

#saveRDS(bwd2o3_bdyint, "intermediate/script02/lqte/bwd2o3_bdyint.RDS")
bwd2o3_bdyint <- readRDS("intermediate/script02/lqte/bwd2o3_bdyint.RDS")

bwd2o3_min <- min(bwd2o3_bdyint[, 2:3]) # MSE-optimal


# 2- Estimate quantiles
## a- Detecting non-zero quantile effects
qte_d2o3_05 <- rd.qte(y = outcome3_d2, x = running2, d = d2, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd2o3_min,
                      tau = seq(0.1, 0.9, by = 0.05))
qte_d2o3_05

qte_d2o3 <- rd.qte(y = outcome3_d2, x = running2, d = d2, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd2o3_min,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Estimating non-zero quantile effects
quantile(indi_ns_ss2$post_interval1318[which(indi_ns_ss2$scoringD2_0 > 0 & indi_ns_ss2$scoringD2_0 < bwd2o3_min)],
         prob=seq(0, 1, length = 21)) # para obtener percentiles en tratados

quantile(indi_ns_ss2$post_interval1318[which(indi_ns_ss2$scoringD2_0 > -bwd2o3_min & indi_ns_ss2$scoringD2_0 < 0)],
         prob=seq(0, 1, length = 21)) # para obtener deciles en no-tratados (el percentil 40 es cero)

qte_d2o3_05b <- rdq.band(y = outcome3_d2, x = running2, d = d2, x0 = 0,
                         cov = 0, bdw = bwd2o3_min, alpha = 0.1,
                         tau = seq(0.45, 0.9, by = 0.05))

## c. Plot
df_d2o3 <- data.frame(
  cuantil = qte_d2o3_05b$tau,
  efecto_estimado = qte_d2o3_05b$qte,
  lower_band = qte_d2o3_05b[["uband.robust"]][, 1, ],
  upper_band = qte_d2o3_05b[["uband.robust"]][, 2, ])

p_lqte_d2y3 <- ggplot(df_d2o3, aes(x = cuantil, y = efecto_estimado)) +
  geom_ribbon(aes(ymin = lower_band, ymax = upper_band), fill = "grey", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 0) +
  theme_light() +
  ylim(-180, +180) +
  labs(x = expression(paste("Quantile (", tau, ")")),
       y = expression(paste("Estimate of LQTE(", tau, ") of ", D[2]," on ", Y[3]))) # Saved at 395x300

ggsave(filename = "intermediate/script02/plots/lqte_d2_y3.svg",  plot = p_lqte_d2y3,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# 3- Hypothesis testing
Wd2o3_signhomo_nz <- rdq.test(y = outcome3_d2, x = running2, d = d2, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.45, 0.9, by = 0.05),
                              bdw = bwd2o3_min, bias = 1, type = c(1,2))

saveRDS(Wd2o3_signhomo_nz, "intermediate/script02/lqte/Wd2o3_signhomo_nz.RDS")


# 4) Writing the table
container_lqte_d2$y3[1:3] <- round(qte_d2o3$qte, 3)

container_lqte_d2$y3[4] <- if(Wd2o3_signhomo_nz$p.value$significance[1] < 0.01) {"***"} else{
  if(Wd2o3_signhomo_nz$p.value$significance[1] < 0.05) {"**"} else{
    if(Wd2o3_signhomo_nz$p.value$significance[1] < 0.1) {"*"} else{""}
  }
}

container_lqte_d2$y3[5] <- if(Wd2o3_signhomo_nz$p.value$homogeneity[1] < 0.01) {"***"} else{
  if(Wd2o3_signhomo_nz$p.value$homogeneity[1] < 0.05) {"**"} else{
    if(Wd2o3_signhomo_nz$p.value$homogeneity[1] < 0.1) {"*"} else{""}
  }
}



### Outcome 1.4: days worked post months [19, 24] ---------------
# 1- Estimate bandwidth
bwd2o4_bdy <- rdq.bandwidth(y = outcome4_d2, x = running2, d = d2, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 1, val = seq(0.01, 0.1, by = 0.01))


bwd2o4_int <- rdq.bandwidth(y = outcome4_d2, x = running2, d = d2, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 0, val = seq(0.01, 0.1, by = 0.01))


bwd2o4_bdyint <- data.frame(cv = c(bwd2o4_bdy$cv, bwd2o4_int$cv),
                            bdy = c(bwd2o4_bdy$opt.m, bwd2o4_bdy$opt.p),
                            int = c(bwd2o4_int$opt.m, bwd2o4_int$opt.p))

#saveRDS(bwd2o4_bdyint, "intermediate/script02/lqte/bwd2o4_bdyint.RDS")
bwd2o4_bdyint <- readRDS("intermediate/script02/lqte/bwd2o4_bdyint.RDS")


bwd2o4_min <- min(bwd2o4_bdyint[, 2:3])

# 2- Estimate quantiles
## a- Detecting non-zero quantile effects
qte_d2o4_05 <- rd.qte(y = outcome4_d2, x = running2, d = d2, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd2o4_min,
                      tau = seq(0.1, 0.9, by = 0.05))

qte_d2o4 <- rd.qte(y = outcome4_d2, x = running2, d = d2, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd2o4_min,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Estimating non-zero quantile effects
quantile(indi_ns_ss2$post_interval1924[which(indi_ns_ss2$scoringD2_0 > 0 & indi_ns_ss2$scoringD2_0 < bwd2o4_min)],
         prob=seq(0, 1, length = 21)) # para obtener percentiles en tratados

quantile(indi_ns_ss2$post_interval1924[which(indi_ns_ss2$scoringD2_0 > -bwd2o4_min & indi_ns_ss2$scoringD2_0 < 0)],
         prob=seq(0, 1, length = 21)) # para obtener deciles en no-tratados (el percentil 40 es cero)

qte_d2o4_05b <- rdq.band(y = outcome4_d2, x = running2, d = d2, x0 = 0,
                         cov = 0, bdw = bwd2o4_min, alpha = 0.1,
                         tau = seq(0.45, 0.9, by = 0.05))

## c. Plot
df_d2o4 <- data.frame(
  cuantil = qte_d2o4_05b$tau,
  efecto_estimado = qte_d2o4_05b$qte,
  lower_band = qte_d2o4_05b[["uband.robust"]][, 1, ],
  upper_band = qte_d2o4_05b[["uband.robust"]][, 2, ])

p_lqte_d2y4 <- ggplot(df_d2o4, aes(x = cuantil, y = efecto_estimado)) +
  geom_ribbon(aes(ymin = lower_band, ymax = upper_band), fill = "grey", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 0.5) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 0) +
  theme_light() +
  ylim(-180, +180) +
  labs(x = expression(paste("Quantile (", tau, ")")),
       y = expression(paste("Estimate of LQTE(", tau, ") of ", D[2]," on ", Y[4]))) # Saved at 395x300

ggsave(filename = "intermediate/script02/plots/lqte_d2_y4.svg",  plot = p_lqte_d2y4,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# 3- Hypothesis testing
Wd2o4_signhomo_nz <- rdq.test(y = outcome4_d2, x = running2, d = d2, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.45, 0.9, by = 0.05),
                              bdw = bwd2o4_min, bias = 1, type = c(1,2))

saveRDS(Wd2o4_signhomo_nz, "intermediate/script02/lqte/Wd2o4_signhomo_nz.RDS")


# 4) Writing the table
container_lqte_d2$y4[1:3] <- round(qte_d2o4$qte, 3)

container_lqte_d2$y4[4] <- if(Wd2o4_signhomo_nz$p.value$significance[1] < 0.01) {"***"} else{
  if(Wd2o4_signhomo_nz$p.value$significance[1] < 0.05) {"**"} else{
    if(Wd2o4_signhomo_nz$p.value$significance[1] < 0.1) {"*"} else{""}
  }
}

container_lqte_d2$y4[5] <- if(Wd2o4_signhomo_nz$p.value$homogeneity[1] < 0.01) {"***"} else{
  if(Wd2o4_signhomo_nz$p.value$homogeneity[1] < 0.05) {"**"} else{
    if(Wd2o4_signhomo_nz$p.value$homogeneity[1] < 0.1) {"*"} else{""}
  }
}

### Saving and exporting --------

#saveRDS(container_lqte_d2, "intermediate/script02/lqte/container_lqte_d2.RDS")
#openxlsx::write.xlsx(container_lqte_d2, 'intermediate/script02/lqte/container_lqte_d2.xlsx')

container_lqte_d2 <- readRDS("intermediate/script02/lqte/container_lqte_d2.RDS")

