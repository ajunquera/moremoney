# ...............................................................................
# ASSEGNO PER IL LAVORO - 03 Plausibility and robustness
# Author: Álvaro F. Junquera (UAB)
# ...............................................................................

library(tidyverse)
library(rdrobust)
library(rddensity)
library(modelsummary)
library(openxlsx)
library(tidytable)
library(descr)

library(QTE.RD)
library(rd.categorical) # package by Ke-Li Xu (2017), follow "installingXuR.R" to install it

library(glue)

# 0. Reading data -------------
indi_ns_ss1 <- readRDS("intermediate/script01/indi_ns_ss1_190225.RDS")
indi_ns_ss2 <- readRDS("intermediate/script01/indi_ns_ss2_190225.RDS")
longest <- readRDS("intermediate/script01/longest_190225.RDS")

indi_ns_ss1_c <- readRDS("intermediate/script02/indi_ns_ss1_c.RDS")
indi_ns_ss2_c <- readRDS("intermediate/script02/indi_ns_ss2_c.RDS")


indi_ns_ss1 <- as.data.frame(indi_ns_ss1)
indi_ns_ss2 <- as.data.frame(indi_ns_ss2)
longest <- as.data.frame(longest)

indi_ns_ss1_c <- as.data.frame(indi_ns_ss1_c)
indi_ns_ss2_c <- as.data.frame(indi_ns_ss2_c)


# 1. Some necessary functions -----------

create_ci <- function(x) {
  paste0("(", round(x[1], 3), ", ", round(x[2], 3), ")")
}


rdcate <- function(score, outcome) {
  # iD1lc <- indi_ns_ss1_c[, c("scoringD1_0", "lcris_OEC")]
  # setnames(iD1lc, "lcris_OEC", "lcris_OEC_yes")
  # iD1lc$lcris_OEC_no <- ifelse(iD1lc$lcris_OEC_yes == 1, 0, 1)
  outcome_no <- ifelse(outcome == 1, 0, 1)
  ilc <- cbind(score, outcome, outcome_no)
  colnames(ilc) <- c("scoring", "outcome_yes", "outcome_no")

  ilc_ord <- as.data.frame(ilc) %>% arrange(scoring) # de forma ascendente
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

  tabla_est <- data.frame(
    pointest = round(rd_cat_est95$ATE, 3),
    ci95 = create_ci(rd_cat_est95$ci_rob),
    signif90 = if (data.table::between(0, rd_cat_est90$ci_rob[1], rd_cat_est90$ci_rob[2]) == T) {
      F
    } else {
      T
    },
    signif95 = if (data.table::between(0, rd_cat_est95$ci_rob[1], rd_cat_est95$ci_rob[2]) == T) {
      F
    } else {
      T
    },
    signif99 = if (data.table::between(0, rd_cat_est99$ci_rob[1], rd_cat_est99$ci_rob[2]) == T) {
      F
    } else {
      T
    }
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


create_ci_vec <- function(x, fila) {
  paste0("(", round(x[fila, 1], 3), ", ", round(x[fila, 2], 3), ")")
}


rdcate_multinom <- function(score, outcome, dataset) {

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


tostata <- function(model) {
  code <- paste0(
    "rdmse ", model$call$y[[3]], " ", model$call$x[[3]], ", ",
    "deriv(0) c(0) twosided pl(", model$call$p, ") pr(", model$call$p,
    ") hl(", model$bws[1, 1], ") hr(", model$bws[1, 2], ") bl(",
    model$bws[2, 1], ") br(", model$bws[2, 2], ") kernel(",
    str_to_lower(model$kernel), ")",
    "\n", "scalar ", model$call$y[[3]], model$call$p, model$kernel, "_l = ", "r(amse_l_cl)",
    "\n", "scalar ", model$call$y[[3]], model$call$p, model$kernel, "_r = ", "r(amse_r_cl)"
  )
  # code2 <- paste0( "scalar ", model$call$y[[3]], "_l = ", "r(amse_l_cl)")

  # cat(code)
  # cat(code2)
}


# 2. PLAUSIBILITY ------------------

## Plausibility checks for LATE ----------------
### D1 ----------------
#### Pseudo-treatments --------------
indi_ns_ss1_rmT1 <- indi_ns_ss1 %>%
  filter(ppa_data_avvio_d > ymd("2017-06-30"))

# haven::write_dta(indi_ns_ss1[, c("scoringD1_0", "prewd1_6", "prewd7_12", "prewd13_18", "prewd19_24")], "indi_ns_ss1_pse_n2023_rmT1.dta")


le_val1_prewd16 <- rdrobust(
  y = indi_ns_ss1_rmT1$prewd1_6, x = indi_ns_ss1_rmT1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(le_val1_prewd16)

le_val1_prewd16_2 <- rdrobust(
  y = indi_ns_ss1_rmT1$prewd1_6, x = indi_ns_ss1_rmT1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(le_val1_prewd16_2)



#### Balance in pretreatment variables ----------------
## Continuous outcome variables (age)
val1_age <- rdrobust(
  y = indi_ns_ss1$eta, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(val1_age)

val1_age2 <- rdrobust(
  y = indi_ns_ss1$eta, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(val1_age2)

val1u_age <- rdrobust(
  y = indi_ns_ss1$eta, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(val1_age)

val1u_age2 <- rdrobust(
  y = indi_ns_ss1$eta, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(val1_age2)


val1tabc <- modelsummary(list(val1_age, val1_age2),
                         statistic = "{p.value} [{conf.low}, {conf.high}]", coef_map = cm,
                         stars = c("*" = .1, "**" = .05, "***" = 0.01),
                         output = "data.frame")

val1tabc <- val1tabc[, -c(1:3)]
val1tabct <- as.data.frame(t(val1tabc))
val1tabct$pvalue <- str_sub(val1tabct$V2, 1, 5)
val1tabct$ci <- str_sub(val1tabct$V2, 7)
val1tabct$V2 <- NULL
colnames(val1tabct) <- c("point", "bw", "n_bw", "pvalue", "ci")
val1tabct$outcome <- c("age", "age")
val1tabct$polydegreeX <- rep(c(1, 2))


## Categorical outcome variables (sex, foreign, studies, disability, lastoccupation, last sector)
# Data preparation
### Sex
indi_ns_ss1$sex_f <- ifelse(indi_ns_ss1$genere == "F", 1, 0)
# indi_ns_ss1$sex_f_f <- as.factor(indi_ns_ss1$sex_f)

v1_sexF <- rdrobust(
  y = indi_ns_ss1$sex_f, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_sexF)

v1_sexF_2 <- rdrobust(
  y = indi_ns_ss1$sex_f, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_sexF_2)


v1u_sexF <- rdrobust(
  y = indi_ns_ss1$sex_f, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_sexF)

v1u_sexF_2 <- rdrobust(
  y = indi_ns_ss1$sex_f, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_sexF_2)

### Foreign
indi_ns_ss1$foreign <- ifelse(indi_ns_ss1$flg_italiano == "no", 1, 0)

v1_foreignF <- rdrobust(
  y = indi_ns_ss1$foreign, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_foreignF)

v1_foreignF_2 <- rdrobust(
  y = indi_ns_ss1$foreign, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_foreignF_2)


v1u_foreignF <- rdrobust(
  y = indi_ns_ss1$foreign, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_foreignF)

v1u_foreignF_2 <- rdrobust(
  y = indi_ns_ss1$foreign, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_foreignF_2)


### Education
indi_ns_ss1 <- indi_ns_ss1 %>% mutate(studio2_groupedX = case_when(
  studio2 == "Licenza elementare" ~ "1_LicElementare",
  studio2 == "Licenza media" ~ "2_LicMedia",
  studio2 %in% c("Diploma di qualifica professionale") ~ "3_QualificaProf",
  studio2 %in% c(
    "Diploma Tecnico", "Diploma Professionale", "Diploma Liceale", "Diploma universitario",
    "Diploma Conservatorio musicale",
    "Maestro d'arte", "Scuola magistrale (triennale)", "Diploma di istruzione artistica", "Diploma interprete, traduttore, archivista"
  ) ~ "4_Diploma",
  studio2 %in% c("Laurea I livello (triennale)", "Laurea - vecchio o nuovo ordinamento", "Post laurea") ~ "5_Laurea",
  TRUE ~ "Nondetermined"
))
freq(indi_ns_ss1$studio2_groupedX)

indi_ns_ss1$elementare <- ifelse(indi_ns_ss1$studio2_groupedX == "1_LicElementare", 1, 0)
indi_ns_ss1$media <- ifelse(indi_ns_ss1$studio2_groupedX == "2_LicMedia", 1, 0)
indi_ns_ss1$qualifica <- ifelse(indi_ns_ss1$studio2_groupedX == "3_QualificaProf", 1, 0)
indi_ns_ss1$diploma <- ifelse(indi_ns_ss1$studio2_groupedX == "4_Diploma", 1, 0)
indi_ns_ss1$laurea <- ifelse(indi_ns_ss1$studio2_groupedX == "5_Laurea", 1, 0)

#### elementare

v1_elementareF <- rdrobust(
  y = indi_ns_ss1$elementare, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_elementareF)

v1_elementareF_2 <- rdrobust(
  y = indi_ns_ss1$elementare, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_elementareF_2)


v1u_elementareF <- rdrobust(
  y = indi_ns_ss1$elementare, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_elementareF)

v1u_elementareF_2 <- rdrobust(
  y = indi_ns_ss1$elementare, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_elementareF_2)

#### media

v1_mediaF <- rdrobust(
  y = indi_ns_ss1$media, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_mediaF)

v1_mediaF_2 <- rdrobust(
  y = indi_ns_ss1$media, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_mediaF_2)


v1u_mediaF <- rdrobust(
  y = indi_ns_ss1$media, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_mediaF)

v1u_mediaF_2 <- rdrobust(
  y = indi_ns_ss1$media, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_mediaF_2)


#### qualifica

v1_qualificaF <- rdrobust(
  y = indi_ns_ss1$qualifica, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_qualificaF)

v1_qualificaF_2 <- rdrobust(
  y = indi_ns_ss1$qualifica, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_qualificaF_2)


v1u_qualificaF <- rdrobust(
  y = indi_ns_ss1$qualifica, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_qualificaF)

v1u_qualificaF_2 <- rdrobust(
  y = indi_ns_ss1$qualifica, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_qualificaF_2)


#### diploma

v1_diplomaF <- rdrobust(
  y = indi_ns_ss1$diploma, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_diplomaF)

v1_diplomaF_2 <- rdrobust(
  y = indi_ns_ss1$diploma, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_diplomaF_2)


v1u_diplomaF <- rdrobust(
  y = indi_ns_ss1$diploma, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_diplomaF)

v1u_diplomaF_2 <- rdrobust(
  y = indi_ns_ss1$diploma, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_diplomaF_2)



#### laurea

v1_laureaF <- rdrobust(
  y = indi_ns_ss1$laurea, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_laureaF)

v1_laureaF_2 <- rdrobust(
  y = indi_ns_ss1$laurea, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_laureaF_2)


v1u_laureaF <- rdrobust(
  y = indi_ns_ss1$laurea, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_laureaF)

v1u_laureaF_2 <- rdrobust(
  y = indi_ns_ss1$laurea, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_laureaF_2)


### Disability
indi_ns_ss1$disa <- ifelse(indi_ns_ss1$disability == "Yes", 1, 0)

v1_disaF <- rdrobust(
  y = indi_ns_ss1$disa, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_disaF)

v1_disaF_2 <- rdrobust(
  y = indi_ns_ss1$disa, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_disaF_2)


v1u_disaF <- rdrobust(
  y = indi_ns_ss1$disa, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1_disaF)

v1u_disaF_2 <- rdrobust(
  y = indi_ns_ss1$disa, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd"
)
summary(v1_disaF_2)

### Occupation of last employment
indi_ns_ss1$qualifica1d <- str_sub(indi_ns_ss1$edu_last_contract, 1, 1)

indi_ns_ss1 <- indi_ns_ss1 %>%
  mutate(qualifica1dm = case_when(
    is.na(qualifica1d) | qualifica1d %in% c("1", "9") ~ "Pre36_Other",
    qualifica1d == "2" ~ "2_Intellectual",
    qualifica1d == "3" ~ "3_Technical",
    qualifica1d == "4" ~ "4_WhiteLowSkilled",
    qualifica1d == "5" ~ "5_WhiteHighSkilled",
    qualifica1d == "6" ~ "6_BlueHighSkilled",
    qualifica1d == "7" ~ "7_BlueMediumSkilled",
    qualifica1d == "8" ~ "8_NonSkilled"
  ))

indi_ns_ss1$o2I <- ifelse(indi_ns_ss1$qualifica1dm == "2_Intellectual", 1, 0)
indi_ns_ss1$o3T <- ifelse(indi_ns_ss1$qualifica1dm == "3_Technical", 1, 0)
indi_ns_ss1$o4WLS <- ifelse(indi_ns_ss1$qualifica1dm == "4_WhiteLowSkilled", 1, 0)
indi_ns_ss1$o5WHS <- ifelse(indi_ns_ss1$qualifica1dm == "5_WhiteHighSkilled", 1, 0)
indi_ns_ss1$o6BHS <- ifelse(indi_ns_ss1$qualifica1dm == "6_BlueHighSkilled", 1, 0)
indi_ns_ss1$o7BMS <- ifelse(indi_ns_ss1$qualifica1dm == "7_BlueMediumSkilled", 1, 0)
indi_ns_ss1$o8NS <- ifelse(indi_ns_ss1$qualifica1dm == "8_NonSkilled", 1, 0)
indi_ns_ss1$oX <- ifelse(indi_ns_ss1$qualifica1dm == "Pre36_Other", 1, 0)

v1u_o2 <- rdrobust(
  y = indi_ns_ss1$o2I, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1u_o2)


v1u_o3 <- rdrobust(
  y = indi_ns_ss1$o3T, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_o3)


v1u_o4 <- rdrobust(
  y = indi_ns_ss1$o4WLS, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_o4)


v1u_o5 <- rdrobust(
  y = indi_ns_ss1$o5WHS, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_o5)


v1u_o6 <- rdrobust(
  y = indi_ns_ss1$o6BHS, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_o6)


v1u_o7 <- rdrobust(
  y = indi_ns_ss1$o7BMS, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_o7)


v1u_o8 <- rdrobust(
  y = indi_ns_ss1$o8NS, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_o8)


v1u_oX <- rdrobust(
  y = indi_ns_ss1$oX, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_oX)

### Sector of last employment
indi_ns_ss1 <- indi_ns_ss1 %>%
  mutate(sectorVLm = case_when(
    sectorVL %in% c("Agricoltura", "Unobserved") ~ "Pre36_Other",
    TRUE ~ sectorVL
  ))

indi_ns_ss1$sAE <- ifelse(indi_ns_ss1$sectorVLm == "AltreIndustrie", 1, 0)
indi_ns_ss1$sAS <- ifelse(indi_ns_ss1$sectorVLm == "AltriServizi", 1, 0)
indi_ns_ss1$sCTL <- ifelse(indi_ns_ss1$sectorVLm == "CommercioTempoLib", 1, 0)
indi_ns_ss1$sC <- ifelse(indi_ns_ss1$sectorVLm == "Costruzione", 1, 0)
indi_ns_ss1$sIL <- ifelse(indi_ns_ss1$sectorVLm == "IngrossoLogistica", 1, 0)
indi_ns_ss1$sMII <- ifelse(indi_ns_ss1$sectorVLm == "MadeInItaly", 1, 0)
indi_ns_ss1$sMM <- ifelse(indi_ns_ss1$sectorVLm == "MetalMeccanico", 1, 0)
indi_ns_ss1$sX <- ifelse(indi_ns_ss1$sectorVLm == "Pre36_Other", 1, 0)
indi_ns_ss1$sSP <- ifelse(indi_ns_ss1$sectorVLm == "ServiziPersona", 1, 0)
indi_ns_ss1$sTA <- ifelse(indi_ns_ss1$sectorVLm == "TerziarioAvanzato", 1, 0)


v1u_sAE <- rdrobust(
  y = indi_ns_ss1$sAE, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1u_sAE)


v1u_sAS <- rdrobust(
  y = indi_ns_ss1$sAS, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_sAS)


v1u_sCTL <- rdrobust(
  y = indi_ns_ss1$sCTL, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_sCTL)


v1u_sC <- rdrobust(
  y = indi_ns_ss1$sC, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_sC)


v1u_sIL <- rdrobust(
  y = indi_ns_ss1$sIL, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_sIL)


v1u_sMII <- rdrobust(
  y = indi_ns_ss1$sMII, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_sMII)

v1u_sMM <- rdrobust(
  y = indi_ns_ss1$sMM, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(v1u_sMM)

v1u_sX <- rdrobust(
  y = indi_ns_ss1$sX, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_sX)


v1u_sSP <- rdrobust(
  y = indi_ns_ss1$sSP, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_sSP)

v1u_sTA <- rdrobust(
  y = indi_ns_ss1$sTA, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(v1u_sTA)

##### Export ---------------
tidy.rdrobust <- function(model, ...) {
  ret <- data.frame(
    term = row.names(model$coef),
    estimate = model$coef[, 1],
    std.error = model$se[, 1],
    p.value = model$pv[3]
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

# Table A6
pretable_E2_d1 <- modelsummary(list(val1u_age, v1u_sexF, v1u_foreignF, v1u_disaF,
                  v1u_elementareF, v1u_mediaF, v1u_qualificaF, v1u_diplomaF, v1u_laureaF,
                  v1u_o2, v1u_o3, v1u_o4, v1u_o5, v1u_o6, v1u_o7, v1u_o8, v1u_oX,
                  v1u_sMII, v1u_sMM, v1u_sAE, v1u_sC, v1u_sCTL, v1u_sIL, v1u_sTA, v1u_sSP, v1u_sAS, v1u_sX),
             statistic = "std.error", coef_map = cm,
             stars = F,
             output = "data.frame")

pretable_E2_d1_est <- t(pretable_E2_d1[1, 4:30])

pretable_E2_d1_vest <- data.frame(var_cat = c("Age", "Sex (Woman)", "Foreign (Yes)", "Disability (Yes)",
                                              "Level of educational attainment",
                                              "  Primary", "  Lower secondary", "  Upper secondary of 3 years",
                                              "  Upper secondary of 4 years", "  University or equivalent",
                                              "Last occupation",
                                              "  Intellectual (2)", "  Technical (3)", "  White collar low-skilled (4)",
                                              "  White collar high-skilled (5)", "  Blue collar high-skilled (6)",
                                              "  Blue collar intermediate-skilled (7)", "  Blue collar low-skilled (8)",
                                              "  Other occupations or unemployed (1,9)",
                                              "Sector of last employment",
                                              "  Made in Italy", "  Metal and mechanical",
                                              "  Other industries", "  Construction", "  Retail and leisure",
                                              "  Wholesale and logistics", "  Advanced third sector",
                                              "  Personal services", "  Other services", "Agriculture or unobserved"),
                                  D1 = c(pretable_E2_d1_est[1:4], NA, pretable_E2_d1_est[5:9], NA,
                                         pretable_E2_d1_est[10:17], NA, pretable_E2_d1_est[18:27]))



## Size of the problem: standardized mean difference for variable AGE
# see specification of the Absolute Standardized Difference (ASD) in
# Zhou, T., Tong, G., Li, F., & Thomas, L. E. (2022). PSweight: An R Package for Propensity Score Weighting Analysis. R Journal, 14(1).

denom_inside <- (var(indi_ns_ss1$eta[indi_ns_ss1$scoringD1_0 > 0 & indi_ns_ss1$scoringD1_0 < val1u_age$bws[1,1]]) +
                   var(indi_ns_ss1$eta[indi_ns_ss1$scoringD1_0 < 0 & indi_ns_ss1$scoringD1_0 > -val1u_age$bws[1,1]])) / 2
asd_age_d1 <- round(abs(val1u_age$Estimate[1]) / sqrt(denom_inside), 4)

pretable_E2_d1_vest$D1[1] <- paste0(pretable_E2_d1_vest$D1[1], " (", asd_age_d1, ")")


### D2 ---------------

#### Pseudo-treatments ----------------
indi_ns_ss2_rmT1 <- indi_ns_ss2 %>%
  filter(ppa_data_avvio_d > ymd("2017-06-30"))

# haven::write_dta(indi_ns_ss2[, c("scoringD2_0", "prewd1_6", "prewd7_12", "prewd13_18", "prewd19_24")], "indi_ns_ss2_pse_n2023_rmT1.dta")


le_val2_prewd16 <- rdrobust(
  y = indi_ns_ss2_rmT1$prewd1_6, x = indi_ns_ss2_rmT1$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(le_val2_prewd16)

le_val2_prewd16_2 <- rdrobust(
  y = indi_ns_ss2_rmT1$prewd1_6, x = indi_ns_ss2_rmT1$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "mserd", cluster = NULL
)
summary(le_val2_prewd16_2)


tableE1a <- data.frame(parameter = c("D_1", "SE", "P-value", "Kernel", "Poly. degree",
                                     "D_2", "SE", "P-value", "Kernel", "Poly. degree"),
                       LATE_Ym1 = c(round(le_val1_prewd16$Estimate[1], 3),
                                    round(le_val1_prewd16$Estimate[3], 3),
                                    round(le_val1_prewd16$pv[3], 3),
                                    "Triangular",
                                    "1",
                                    round(le_val2_prewd16$Estimate[1], 3),
                                    round(le_val2_prewd16$Estimate[3], 3),
                                    round(le_val2_prewd16$pv[3], 3),
                                    "Triangular",
                                    "1"))


#### Balance in pretreatment variables ----------------
tidy.rdrobust <- function(model, ...) {
  ret <- data.frame(
    term = row.names(model$coef),
    estimate = model$coef[, 1],
    p.value = model$pv[3],
    conf.low = model$ci[3, 1],
    conf.high = model$ci[3, 2]
  )
  row.names(ret) <- NULL
  ret
}

# glance.rdrobust <- function(model, ...) {
#  ret <- data.frame(
#    Bandwidth = model$bws[1],
#    N_h = sum(model$N_h)
#  )
#  ret
# }

## Continuous outcome variables (age)
vd2u_age <- rdrobust(
  y = indi_ns_ss2$eta, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_age)



## Categorical outcome variables (sex, foreign, studies, disability, lastoccupation, last sector)
# Data preparation
### Sex
indi_ns_ss2$sex_f <- ifelse(indi_ns_ss2$genere == "F", 1, 0)
# indi_ns_ss2$sex_f_f <- as.factor(indi_ns_ss2$sex_f)


vd2u_sexF <- rdrobust(
  y = indi_ns_ss2$sex_f, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(vd2u_sexF)

### Foreign
indi_ns_ss2$foreign <- ifelse(indi_ns_ss2$flg_italiano == "no", 1, 0)

vd2u_foreignF <- rdrobust(
  y = indi_ns_ss2$foreign, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(vd2u_foreignF)

### Education
indi_ns_ss2 <- indi_ns_ss2 %>% mutate(studio2_groupedX = case_when(
  studio2 == "Licenza elementare" ~ "1_LicElementare",
  studio2 == "Licenza media" ~ "2_LicMedia",
  studio2 %in% c("Diploma di qualifica professionale") ~ "3_QualificaProf",
  studio2 %in% c(
    "Diploma Tecnico", "Diploma Professionale", "Diploma Liceale", "Diploma universitario",
    "Diploma Conservatorio musicale",
    "Maestro d'arte", "Scuola magistrale (triennale)", "Diploma di istruzione artistica", "Diploma interprete, traduttore, archivista"
  ) ~ "4_Diploma",
  studio2 %in% c("Laurea I livello (triennale)", "Laurea - vecchio o nuovo ordinamento", "Post laurea") ~ "5_Laurea",
  TRUE ~ "Nondetermined"
))
freq(indi_ns_ss2$studio2_groupedX)

indi_ns_ss2$elementare <- ifelse(indi_ns_ss2$studio2_groupedX == "1_LicElementare", 1, 0)
indi_ns_ss2$media <- ifelse(indi_ns_ss2$studio2_groupedX == "2_LicMedia", 1, 0)
indi_ns_ss2$qualifica <- ifelse(indi_ns_ss2$studio2_groupedX == "3_QualificaProf", 1, 0)
indi_ns_ss2$diploma <- ifelse(indi_ns_ss2$studio2_groupedX == "4_Diploma", 1, 0)
indi_ns_ss2$laurea <- ifelse(indi_ns_ss2$studio2_groupedX == "5_Laurea", 1, 0)

#### elementare


vd2u_elementareF <- rdrobust(
  y = indi_ns_ss2$elementare, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(vd2u_elementareF)

#### media

vd2u_mediaF <- rdrobust(
  y = indi_ns_ss2$media, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(vd2u_mediaF)


#### qualifica


vd2u_qualificaF <- rdrobust(
  y = indi_ns_ss2$qualifica, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(vd2u_qualificaF)



#### diploma


vd2u_diplomaF <- rdrobust(
  y = indi_ns_ss2$diploma, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(vd2u_diplomaF)



#### laurea

vd2u_laureaF <- rdrobust(
  y = indi_ns_ss2$laurea, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(vd2u_laureaF)



### Disability
indi_ns_ss2$disa <- ifelse(indi_ns_ss2$disability == "Yes", 1, 0)


vd2u_disaF <- rdrobust(
  y = indi_ns_ss2$disa, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(vd2u_disaF)

### Occupation of last employment
indi_ns_ss2$qualifica1d <- str_sub(indi_ns_ss2$edu_last_contract, 1, 1)

indi_ns_ss2 <- indi_ns_ss2 %>%
  mutate(qualifica1dm = case_when(
    is.na(qualifica1d) | qualifica1d %in% c("1", "9") ~ "Pre36_Other",
    qualifica1d == "2" ~ "2_Intellectual",
    qualifica1d == "3" ~ "3_Technical",
    qualifica1d == "4" ~ "4_WhiteLowSkilled",
    qualifica1d == "5" ~ "5_WhiteHighSkilled",
    qualifica1d == "6" ~ "6_BlueHighSkilled",
    qualifica1d == "7" ~ "7_BlueMediumSkilled",
    qualifica1d == "8" ~ "8_NonSkilled"
  ))

indi_ns_ss2$o2I <- ifelse(indi_ns_ss2$qualifica1dm == "2_Intellectual", 1, 0)
indi_ns_ss2$o3T <- ifelse(indi_ns_ss2$qualifica1dm == "3_Technical", 1, 0)
indi_ns_ss2$o4WLS <- ifelse(indi_ns_ss2$qualifica1dm == "4_WhiteLowSkilled", 1, 0)
indi_ns_ss2$o5WHS <- ifelse(indi_ns_ss2$qualifica1dm == "5_WhiteHighSkilled", 1, 0)
indi_ns_ss2$o6BHS <- ifelse(indi_ns_ss2$qualifica1dm == "6_BlueHighSkilled", 1, 0)
indi_ns_ss2$o7BMS <- ifelse(indi_ns_ss2$qualifica1dm == "7_BlueMediumSkilled", 1, 0)
indi_ns_ss2$o8NS <- ifelse(indi_ns_ss2$qualifica1dm == "8_NonSkilled", 1, 0)
indi_ns_ss2$oX <- ifelse(indi_ns_ss2$qualifica1dm == "Pre36_Other", 1, 0)

vd2u_o2 <- rdrobust(
  y = indi_ns_ss2$o2I, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(vd2u_o2)


vd2u_o3 <- rdrobust(
  y = indi_ns_ss2$o3T, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_o3)


vd2u_o4 <- rdrobust(
  y = indi_ns_ss2$o4WLS, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_o4)


vd2u_o5 <- rdrobust(
  y = indi_ns_ss2$o5WHS, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_o5)


vd2u_o6 <- rdrobust(
  y = indi_ns_ss2$o6BHS, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_o6)


vd2u_o7 <- rdrobust(
  y = indi_ns_ss2$o7BMS, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_o7)


vd2u_o8 <- rdrobust(
  y = indi_ns_ss2$o8NS, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_o8)


vd2u_oX <- rdrobust(
  y = indi_ns_ss2$oX, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_oX)

### Sector of last employment
indi_ns_ss2 <- indi_ns_ss2 %>%
  mutate(sectorVLm = case_when(
    sectorVL %in% c("Agricoltura", "Unobserved") ~ "Pre36_Other",
    TRUE ~ sectorVL
  ))

indi_ns_ss2$sAE <- ifelse(indi_ns_ss2$sectorVLm == "AltreIndustrie", 1, 0)
indi_ns_ss2$sAS <- ifelse(indi_ns_ss2$sectorVLm == "AltriServizi", 1, 0)
indi_ns_ss2$sCTL <- ifelse(indi_ns_ss2$sectorVLm == "CommercioTempoLib", 1, 0)
indi_ns_ss2$sC <- ifelse(indi_ns_ss2$sectorVLm == "Costruzione", 1, 0)
indi_ns_ss2$sIL <- ifelse(indi_ns_ss2$sectorVLm == "IngrossoLogistica", 1, 0)
indi_ns_ss2$sMII <- ifelse(indi_ns_ss2$sectorVLm == "MadeInItaly", 1, 0)
indi_ns_ss2$sMM <- ifelse(indi_ns_ss2$sectorVLm == "MetalMeccanico", 1, 0)
indi_ns_ss2$sX <- ifelse(indi_ns_ss2$sectorVLm == "Pre36_Other", 1, 0)
indi_ns_ss2$sSP <- ifelse(indi_ns_ss2$sectorVLm == "ServiziPersona", 1, 0)
indi_ns_ss2$sTA <- ifelse(indi_ns_ss2$sectorVLm == "TerziarioAvanzato", 1, 0)


vd2u_sAE <- rdrobust(
  y = indi_ns_ss2$sAE, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)

summary(vd2u_sAE)


vd2u_sAS <- rdrobust(
  y = indi_ns_ss2$sAS, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_sAS)


vd2u_sCTL <- rdrobust(
  y = indi_ns_ss2$sCTL, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_sCTL)


vd2u_sC <- rdrobust(
  y = indi_ns_ss2$sC, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_sC)


vd2u_sIL <- rdrobust(
  y = indi_ns_ss2$sIL, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_sIL)


vd2u_sMII <- rdrobust(
  y = indi_ns_ss2$sMII, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_sMII)


vd2u_sMM <- rdrobust(
  y = indi_ns_ss2$sMM, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_sMM)


vd2u_sX <- rdrobust(
  y = indi_ns_ss2$sX, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_sX)


vd2u_sSP <- rdrobust(
  y = indi_ns_ss2$sSP, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_sSP)

vd2u_sTA <- rdrobust(
  y = indi_ns_ss2$sTA, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd"
)
summary(vd2u_sTA)

##### Export -------------
pretable_E2_d2 <- modelsummary(list(vd2u_age, vd2u_sexF, vd2u_foreignF, vd2u_disaF,
                                    vd2u_elementareF, vd2u_mediaF, vd2u_qualificaF, vd2u_diplomaF, vd2u_laureaF,
                                    vd2u_o2, vd2u_o3, vd2u_o4, vd2u_o5, vd2u_o6, vd2u_o7, vd2u_o8, vd2u_oX,
                                    vd2u_sMII, vd2u_sMM, vd2u_sAE, vd2u_sC, vd2u_sCTL, vd2u_sIL, vd2u_sTA, vd2u_sSP, vd2u_sAS, vd2u_sX),
                               statistic = "std.error", coef_map = cm,
                               stars = F,
                               output = "data.frame")

pretable_E2_d2_est <- t(pretable_E2_d2[1, 4:30])

pretable_E2_d2_vest <- data.frame(var_cat = c("Age", "Sex (Woman)", "Foreign (Yes)", "Disability (Yes)",
                                              "Level of educational attainment",
                                              "  Primary", "  Lower secondary", "  Upper secondary of 3 years",
                                              "  Upper secondary of 4 years", "  University or equivalent",
                                              "Last occupation",
                                              "  Intellectual (2)", "  Technical (3)", "  White collar low-skilled (4)",
                                              "  White collar high-skilled (5)", "  Blue collar high-skilled (6)",
                                              "  Blue collar intermediate-skilled (7)", "  Blue collar low-skilled (8)",
                                              "  Other occupations or unemployed (1,9)",
                                              "Sector of last employment",
                                              "  Made in Italy", "  Metal and mechanical",
                                              "  Other industries", "  Construction", "  Retail and leisure",
                                              "  Wholesale and logistics", "  Advanced third sector",
                                              "  Personal services", "  Other services", "Agriculture or unobserved"),
                                  D2 = c(pretable_E2_d2_est[1:4], NA, pretable_E2_d2_est[5:9], NA,
                                         pretable_E2_d2_est[10:17], NA, pretable_E2_d2_est[18:27]))

## Size of the problem
denom_inside2 <- (var(indi_ns_ss2$eta[indi_ns_ss2$scoringD2_0 > 0 & indi_ns_ss2$scoringD2_0 < vd2u_age$bws[1,1]]) +
                    var(indi_ns_ss2$eta[indi_ns_ss2$scoringD2_0 < 0 & indi_ns_ss2$scoringD2_0 > -vd2u_age$bws[1,1]])) / 2

asd_age_d2 <- round(abs(vd2u_age$Estimate[1]) / sqrt(denom_inside2), 4)

pretable_E2_d2_vest$D2[1] <- paste0(pretable_E2_d2_vest$D2[1], " (", asd_age_d2, ")")


## Export all
table_E2 <- left_join(pretable_E2_d1_vest, pretable_E2_d2_vest, by = "var_cat")

if (!file.exists("intermediate/script03/T_E2_balance.xlsx")) {

  openxlsx::write.xlsx(table_E2, 'intermediate/script03/T_E2_balance.xlsx')

}

## Plausibility checks for LQTE -----------
### D1 --------------
indi_ns_ss1_rmT1$treatedD1 <- ifelse(indi_ns_ss1_rmT1$scoringD1_0 > 0, TRUE, FALSE)

running1rmt1 <- indi_ns_ss1_rmT1$scoringD1_0
outcomepre1 <- indi_ns_ss1_rmT1$prewd1_6
d1rmt1 <- indi_ns_ss1_rmT1$treatedD1
# x0 <- 0


# 1) Estimate bandwidth
bwd1preo1_bdy <- rdq.bandwidth(y = outcomepre1, x = running1rmt1, d = d1rmt1, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 1, val = seq(0.01, 0.1, by = 0.01))

bwd1preo1_int <- rdq.bandwidth(y = outcomepre1, x = running1rmt1, d = d1rmt1, x0 = 0,
                            cv = 1, cov = 0, pm.each = 0,
                            bdy = 0, val = seq(0.01, 0.1, by = 0.01))

bwd1preo1_bdyint <- data.frame(cv = c(bwd1preo1_bdy$cv, bwd1preo1_int$cv),
                            bdy = c(bwd1preo1_bdy$opt.m, bwd1preo1_bdy$opt.p),
                            int = c(bwd1preo1_int$opt.m, bwd1preo1_int$opt.p))

bwd1preo1_min <- min(bwd1preo1_bdyint[, 2:3])



# 2) Estimate quantiles
## a. Detecting quantile effects == 0
qte_d1preo1_05 <- rd.qte(y = outcomepre1, x = running1rmt1, d = d1rmt1, x0 = 0,
                         cov = 0, bias = 0, bdw = bwd1preo1_min,
                         tau = seq(0.1, 0.9, by = 0.05))
qte_d1preo1_05


## b. Estimable quantile effects
qte_d1preo1_est <- rdq.band(y = outcomepre1, x = running1rmt1, d = d1rmt1, x0 = 0,
                            cov = 0, bdw = bwd1preo1_min, alpha = 0.1,
                            tau = seq(0.35, 0.9, by = 0.05))


# 3) Run hypothesis tests
Wd1preo1_sh <- rdq.test(y = outcomepre1, x = running1rmt1, d = d1rmt1, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.3, 0.75, by = 0.05),
                              bdw = bwd1preo1_min, bias = 1, type = c(1,2))



### D2 --------------
indi_ns_ss2_rmT1$treatedD2 <- ifelse(indi_ns_ss2_rmT1$scoringD2_0 > 0, TRUE, FALSE)

running2rmt1 <- indi_ns_ss2_rmT1$scoringD2_0
outcomepre1_d2 <- indi_ns_ss2_rmT1$prewd1_6
d2rmt1 <- indi_ns_ss2_rmT1$treatedD2

# 1) Estimate bandwidth
bwd2preo1_bdy <- rdq.bandwidth(y = outcomepre1_d2, x = running2rmt1, d = d2rmt1, x0 = 0,
                               cv = 1, cov = 0, pm.each = 0,
                               bdy = 1, val = seq(0.01, 0.1, by = 0.01))

bwd2preo1_int <- rdq.bandwidth(y = outcomepre1_d2, x = running2rmt1, d = d2rmt1, x0 = 0,
                               cv = 1, cov = 0, pm.each = 0,
                               bdy = 0, val = seq(0.01, 0.1, by = 0.01))

bwd2preo1_bdyint <- data.frame(cv = c(bwd2preo1_bdy$cv, bwd2preo1_int$cv),
                               bdy = c(bwd2preo1_bdy$opt.m, bwd2preo1_bdy$opt.p),
                               int = c(bwd2preo1_int$opt.m, bwd2preo1_int$opt.p))

bwd2preo1_min <- min(bwd2preo1_bdyint[, 2:3])




# 2) Estimate quantiles
## a. Detecting quantile effects == 0
qte_d2preo1 <- rd.qte(y = outcomepre1_d2, x = running2rmt1, d = d2rmt1, x0 = 0,
                         cov = 0, bias = 0, bdw = bwd2preo1_min,
                         tau = seq(0.1, 0.9, by = 0.05))
qte_d2preo1


## b. Estimable quantile effects
qte_d2preo1_est1 <- rd.qte(y = outcomepre1_d2, x = running2rmt1, d = d2rmt1, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd2preo1_min,
                      tau = seq(0.5, 0.9, by = 0.05))
qte_d2preo1_est1

qte_d2preo1_est <- rdq.band(y = outcomepre1_d2, x = running2rmt1, d = d2rmt1, x0 = 0,
                            cov = 0, bdw = bwd2preo1_min, alpha = 0.1,
                            tau = seq(0.4, 0.9, by = 0.05)) # it doesn't work


# 3) Run hypothesis tests
Wd2preo1_sh <- rdq.test(y = outcomepre1_d2, x = running2rmt1, d = d2rmt1, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.4, 0.9, by = 0.05),
                              bdw = 0.05, bias = 1, type = c(1,2)) # it doesn't work

indi_ns_ss2_rmT1 %>%
  filter(scoringD2_0 > -0.1 & scoringD2_0 < 0) %>%
  tabyl(prewd1_6) # 49 % have 0 days of employment in the control group

indi_ns_ss2_rmT1 %>%
  filter(scoringD2_0 > 0 & scoringD2_0 < 0.1) %>%
  tabyl(prewd1_6) # 57 % have 0 days of employment in the treatment group

## We cannot run hypothesis tests for this case because the bandwidth is estimated
## with the difference at the median. However, the median at the treatment group is 0,
## so there aren't enough data to estimate the optimal bandwidth.

tableE1b <- data.frame(parameter = c("D_1", "SE", "P-value", "Kernel", "Poly. degree",
                                     "D_2", "SE", "P-value", "Kernel", "Poly. degree"),
                       LQTE05_Ym1 = c(round(qte_d1preo1_est$qte[4], 3),
                                    "",
                                    paste0("TS: ", Wd1preo1_sh$p.value[1], "\nTH: ", Wd1preo1_sh$p.value[2]),
                                    "Epanechnikov",
                                    "1",
                                    round(qte_d2preo1_est1$qte[1], 3),
                                    "",
                                    "",
                                    "Epanechnikov",
                                    "1"))


### Export ----------

tableE1 <- bind_cols(tableE1a, LQTE05_Ym1 = tableE1b$LQTE05_Ym1)

if (!file.exists("intermediate/script03/T_E1_pseudotreatments.xlsx")) {

  openxlsx::write.xlsx(tableE1, 'intermediate/script03/T_E1_pseudotreatments.xlsx')

}

# 3. ROBUSTNESS checks ----------------

## In the article, we vary four elements of the model as robustness checks:
### polynomial degree, slight changes of the optimal common bandwidth and different optimal bandwidths at each side of the cutoff

## For LATE --------------

glance.rdrobust <- function(model, ...) {
  ret <- data.frame(
    Bandwidths = paste(round(model$bws[1, 1], 3), round(model$bws[1, 2], 3)),
    N_h = sum(model$N_h),
    Kernel = model$kernel,
    Polynomial_degree = round(model$p)
  )
  ret
}

### D1 -----------
#### Changing the bandwidth selector (one selector for each side) ------------

### Outcome 1: days worked post months [1, 6]
rdd_2_D1ei6_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_D1ei6_p1kT)

rdd_2_D1ei6_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_D1ei6_p2kT)


if (!file.exists("intermediate/script03/E15_2bws_D1_Y1.docx")) {

  modelsummary(list(rdd_2_D1ei6_p1kT, rdd_2_D1ei6_p2kT),
               statistic = "std.error", coef_map = cm,
               stars = c('*' = .1, '**' = .05, '***' = 0.01),
               output = "intermediate/script03/E15_2bws_D1_Y1.docx")
}



### Outcome 2: days worked post months [7, 12]
rdd_2_D1ei12_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval712, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", c = 0, p = 1, bwselect = "msetwo",
  cluster = NULL
)
summary(rdd_2_D1ei12_p1kT)

rdd_2_D1ei12_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval712, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_D1ei12_p2kT)

if (!file.exists("intermediate/script03/E16_2bws_D1_Y2.docx")) {

  modelsummary(list(rdd_2_D1ei12_p1kT, rdd_2_D1ei12_p2kT),
               statistic = "std.error", coef_map = cm,
               stars = c('*' = .1, '**' = .05, '***' = 0.01),
               output = "intermediate/script03/E16_2bws_D1_Y2.docx")
}


### Outcome 3: days worked post [13, 18] months
rdd_2_D1ei18_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval1318, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", c = 0, p = 1, bwselect = "msetwo",
  cluster = NULL
)
summary(rdd_2_D1ei18_p1kT)

rdd_2_D1ei18_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval1318, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular", c = 0, p = 2, bwselect = "msetwo",
  cluster = NULL
)
summary(rdd_2_D1ei18_p2kT)

if (!file.exists("intermediate/script03/E17_2bws_D1_Y3.docx")) {

  modelsummary(list(rdd_2_D1ei18_p1kT, rdd_2_D1ei18_p2kT),
               statistic = "std.error", coef_map = cm,
               stars = c('*' = .1, '**' = .05, '***' = 0.01),
               output = "intermediate/script03/E17_2bws_D1_Y3.docx")
}


### Outcome 4: days worked post [19, 24] months
rdd_2_24_p1kT <- rdrobust(
  y = indi_ns_ss1$post_interval1924, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_24_p1kT)

rdd_2_24_p2kT <- rdrobust(
  y = indi_ns_ss1$post_interval1924, x = indi_ns_ss1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_24_p2kT)


if (!file.exists("intermediate/script03/E18_2bws_D1_Y4.docx")) {

  modelsummary(list(rdd_2_24_p1kT, rdd_2_24_p2kT),
               statistic = "std.error", coef_map = cm,
               stars = c('*' = .1, '**' = .05, '***' = 0.01),
               output = "intermediate/script03/E18_2bws_D1_Y4.docx")
}


##### Stata code to estimate AMSE --------

# Function to write Stata code for AMSE estimation
paste0(
  "rmdse ", rdd_2_D1ei6_p1kT$call$y[[3]], " ", rdd_2_D1ei6_p1kT$call$x[[3]], ", ",
  "deriv(0) c(0) pl(", rdd_2_D1ei6_p1kT$call$p, ") pr(", rdd_2_D1ei6_p1kT$call$p,
  ") hl(", rdd_2_D1ei6_p1kT$bws[1, 1], ") hr(", rdd_2_D1ei6_p1kT$bws[1, 2], ") bl(",
  rdd_2_D1ei6_p1kT$bws[2, 1], ") br(", rdd_2_D1ei6_p1kT$bws[2, 2], ") kernel(",
  str_to_lower(rdd_2_D1ei6_p1kT$kernel), ")"
)


# Y1
robust_stata_d1_1 <- data.frame(stata = c(tostata(rdd_2_D1ei6_p1kT), tostata(rdd_2_D1ei6_p2kT)))

# Y2
robust_stata_d1_2 <- data.frame(stata = c(tostata(rdd_2_D1ei12_p1kT), tostata(rdd_2_D1ei12_p2kT)))

# Y3
robust_stata_d1_3 <- data.frame(stata = c(tostata(rdd_2_D1ei18_p1kT), tostata(rdd_2_D1ei18_p2kT)))

# Y4
robust_stata_d1_4 <- data.frame(stata = c(tostata(rdd_2_24_p1kT), tostata(rdd_2_24_p2kT)))

amse_script_2bws_d1 <- glue("net install rdmse, from(https://raw.githubusercontent.com/peizhuan/rdmse/master) replace
                       use \"C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script01/indi_ns_stata_D1.dta\"

                      matrix res = J(1, 16, .) // creates the matrix to save results


                       * Y1.1
                       {robust_stata_d1_1$stata[1]}
                       matrix res[1,1] = r(amse_l_cl)
                       matrix res[1,2] = r(amse_r_cl)

                       {robust_stata_d1_1$stata[2]}
                       matrix res[1,3] = r(amse_l_cl)
                       matrix res[1,4] = r(amse_r_cl)


                       * Y1.2
                       {robust_stata_d1_2$stata[1]}
                       matrix res[1,5] = r(amse_l_cl)
                       matrix res[1,6] = r(amse_r_cl)

                       {robust_stata_d1_2$stata[2]}
                       matrix res[1,7] = r(amse_l_cl)
                       matrix res[1,8] = r(amse_r_cl)

                       * Y1.3
                       {robust_stata_d1_3$stata[1]}
                       matrix res[1,9] = r(amse_l_cl)
                       matrix res[1,10] = r(amse_r_cl)

                       {robust_stata_d1_3$stata[2]}
                       matrix res[1,11] = r(amse_l_cl)
                       matrix res[1,12] = r(amse_r_cl)

                       * Y1.4
                       {robust_stata_d1_4$stata[1]}
                       matrix res[1,13] = r(amse_l_cl)
                       matrix res[1,14] = r(amse_r_cl)

                       {robust_stata_d1_4$stata[2]}
                       matrix res[1,15] = r(amse_l_cl)
                       matrix res[1,16] = r(amse_r_cl)

                       * Export to Excel
                       // Crear un dataset temporal con solo la matriz
                       clear // Asegurarse de que no haya datos existentes
                       svmat res, names(col)

                       // Exportar el dataset como Excel
                       export excel using \"C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script03/stata/resu_amse_1_2bws.xlsx\", firstrow(variables) replace

                       ")

writeLines(amse_script_2bws_d1, "intermediate/script03/amse_script_2bws_d1.txt")



#### Robustness to residual imbalance by adding controls ---------------
indi_ns_ss1 <- indi_ns_ss1 %>%
  mutate(eta2 = eta^2)

indi_ns_ss1_covs <- model.matrix(~indi_ns_ss1$eta +
                                   indi_ns_ss1$eta2 +
                                   indi_ns_ss1$studio2_groupedX +
                                   indi_ns_ss1$sectorVLm - 1) %>%
  as.data.frame()

### Outcome 1.1: days worked post months [1, 6]
rddD1ei6_p1kT_ri <- rdrobust(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0,
  covs = indi_ns_ss1_covs,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL)

summary(rddD1ei6_p1kT_ri)

### Outcome 1.2: days worked post months [7, 12]
rddD1ei12_p1kT_ri <- rdrobust(
  y = indi_ns_ss1$post_interval712, x = indi_ns_ss1$scoringD1_0,
  covs = indi_ns_ss1_covs,
  kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
  cluster = NULL
)

summary(rddD1ei12_p1kT_ri)

### Outcome 1.3: days worked post [13, 18] months
rddD1ei18_p1kT_ri <- rdrobust(
  y = indi_ns_ss1$post_interval1318, x = indi_ns_ss1$scoringD1_0,
  covs = indi_ns_ss1_covs,
  kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
  cluster = NULL
)
summary(rddD1ei18_p1kT_ri)

### Outcome 1.4: days worked post [19, 24] months
rdd24_p1kT_ri <- rdrobust(
  y = indi_ns_ss1$post_interval1924, x = indi_ns_ss1$scoringD1_0,
  covs = indi_ns_ss1_covs,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(rdd24_p1kT_ri)

### D2 -------------
####  Changing the bandwidth selector (one selector for each side) -----------

### Outcome 1: days worked post 6 months (rdd_2_d2_)
rdd_2_d2_ei6_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval6, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_d2_ei6_p1kT)

rdd_2_d2_ei6_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval6, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_d2_ei6_p2kT)

if (!file.exists("intermediate/script03/E19_2bws_D2_Y1.docx")) {

  modelsummary(list(rdd_2_d2_ei6_p1kT, rdd_2_d2_ei6_p2kT),
               statistic = "std.error", coef_map = cm,
               stars = c('*' = .1, '**' = .05, '***' = 0.01),
               output = "intermediate/script03/E19_2bws_D2_Y1.docx")
}


### Outcome 2: days worked post [7, 12] months
rdd_2_d2_ei12_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval712, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_d2_ei12_p1kT)

rdd_2_d2_ei12_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval712, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_d2_ei12_p2kT)

if (!file.exists("intermediate/script03/E20_2bws_D2_Y2.docx")) {

  modelsummary(list(rdd_2_d2_ei12_p1kT, rdd_2_d2_ei12_p2kT),
               statistic = "std.error", coef_map = cm,
               stars = c('*' = .1, '**' = .05, '***' = 0.01),
               output = "intermediate/script03/E20_2bws_D2_Y2.docx")
}


### Outcome 3: days worked post [13, 18] months
rdd_2_d2_ei18_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval1318, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_d2_ei18_p1kT)

rdd_2_d2_ei18_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval1318, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_d2_ei18_p2kT)

if (!file.exists("intermediate/script03/E21_2bws_D2_Y3.docx")) {

  modelsummary(list(rdd_2_d2_ei18_p1kT, rdd_2_d2_ei18_p2kT),
               statistic = "std.error", coef_map = cm,
               stars = c('*' = .1, '**' = .05, '***' = 0.01),
               output = "intermediate/script03/E21_2bws_D2_Y3.docx")
}


### Outcome 4: days worked post [19, 24] months (rdd_2_d2_)
rdd_2_d2_24_p1kT <- rdrobust(
  y = indi_ns_ss2$post_interval1924, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_d2_24_p1kT)

rdd_2_d2_24_p2kT <- rdrobust(
  y = indi_ns_ss2$post_interval1924, x = indi_ns_ss2$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 2, bwselect = "msetwo", cluster = NULL
)
summary(rdd_2_d2_24_p2kT)

if (!file.exists("intermediate/script03/E22_2bws_D2_Y4.docx")) {

  modelsummary(list(rdd_2_d2_24_p1kT, rdd_2_d2_24_p2kT),
               statistic = "std.error", coef_map = cm,
               stars = c('*' = .1, '**' = .05, '***' = 0.01),
               output = "intermediate/script03/E22_2bws_D2_Y4.docx")
}


##### Stata code to estimate AMSE -----------
# Y1
robust_stata_d2_1 <- data.frame(stata = c(tostata(rdd_2_d2_ei6_p1kT), tostata(rdd_2_d2_ei6_p2kT)))

# Y2
robust_stata_d2_2 <- data.frame(stata = c(tostata(rdd_2_d2_ei12_p1kT), tostata(rdd_2_d2_ei12_p2kT)))

# Y3
robust_stata_d2_3 <- data.frame(stata = c(tostata(rdd_2_d2_ei18_p1kT), tostata(rdd_2_d2_ei18_p2kT)))

# Y5
robust_stata_d2_4 <- data.frame(stata = c(tostata(rdd_2_d2_24_p1kT), tostata(rdd_2_d2_24_p2kT)))

amse_script_2bws_d2 <- glue("net install rdmse, from(https://raw.githubusercontent.com/peizhuan/rdmse/master) replace
                       use \"C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script01/indi_ns_stata_d2.dta\"

                      matrix res = J(1, 16, .) // creates the matrix to save results


                       * Y1.1
                       {robust_stata_d2_1$stata[1]}
                       matrix res[1,1] = r(amse_l_cl)
                       matrix res[1,2] = r(amse_r_cl)

                       {robust_stata_d2_1$stata[2]}
                       matrix res[1,3] = r(amse_l_cl)
                       matrix res[1,4] = r(amse_r_cl)


                       * Y1.2
                       {robust_stata_d2_2$stata[1]}
                       matrix res[1,5] = r(amse_l_cl)
                       matrix res[1,6] = r(amse_r_cl)

                       {robust_stata_d2_2$stata[2]}
                       matrix res[1,7] = r(amse_l_cl)
                       matrix res[1,8] = r(amse_r_cl)

                       * Y1.3
                       {robust_stata_d2_3$stata[1]}
                       matrix res[1,9] = r(amse_l_cl)
                       matrix res[1,10] = r(amse_r_cl)

                       {robust_stata_d2_3$stata[2]}
                       matrix res[1,11] = r(amse_l_cl)
                       matrix res[1,12] = r(amse_r_cl)

                       * Y1.4
                       {robust_stata_d2_4$stata[1]}
                       matrix res[1,13] = r(amse_l_cl)
                       matrix res[1,14] = r(amse_r_cl)

                       {robust_stata_d2_4$stata[2]}
                       matrix res[1,15] = r(amse_l_cl)
                       matrix res[1,16] = r(amse_r_cl)

                       * Export to Excel
                       // Crear un dataset temporal con solo la matriz
                       clear // Asegurarse de que no haya datos existentes
                       svmat res, names(col)

                       // Exportar el dataset como Excel
                       export excel using \"C:/Users/afernan5/Nextcloud/AFJ/AXL/intermediate/script03/stata/resu_amse_2_2bws.xlsx\", firstrow(variables) replace

                       ")

writeLines(amse_script_2bws_d2, "intermediate/script03/amse_script_2bws_d2.txt")




#### Robustness to residual imbalance by adding controls ---------
indi_ns_ss2 <- indi_ns_ss2 %>%
  mutate(eta2 = eta^2)

indi_ns_ss2_covs <- model.matrix(~indi_ns_ss2$eta +
                                   indi_ns_ss2$eta2 +
                                   indi_ns_ss2$qualifica1dm - 1) %>%
  as.data.frame()

### Outcome 1.1: days worked post months [1, 6]
rddei6_p1kT_ri <- rdrobust(
  y = indi_ns_ss2$post_interval6, x = indi_ns_ss2$scoringD2_0,
  covs = indi_ns_ss2_covs,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)

summary(rddei6_p1kT_ri)

### Outcome 1.2: days worked post [7, 12] months
rddei12_p1kT_ri <- rdrobust(
  y = indi_ns_ss2$post_interval712, x = indi_ns_ss2$scoringD2_0,
  covs = indi_ns_ss2_covs,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(rddei12_p1kT_ri)

### Outcome 1.3: days worked post [13, 18] months
rddei18_p1kT_ri <- rdrobust(
  y = indi_ns_ss2$post_interval1318, x = indi_ns_ss2$scoringD2_0,
  covs = indi_ns_ss2_covs,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(rddei18_p1kT_ri)

### Outcome 1.4: days worked post [19, 24] months (RDD)
rddei24_p1kT_ri <- rdrobust(
  y = indi_ns_ss2$post_interval1924, x = indi_ns_ss2$scoringD2_0,
  covs = indi_ns_ss2_covs,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL
)
summary(rddei24_p1kT_ri)


## For LQTE --------------

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


container_lqte_d1 <- data.frame(object = c("q1_pe", "q2_pe", "q3_pe", "nonsign_test_stars", "homo_test_stars", "Bandwidth"),
                                y1 = rep(NA, 6),
                                y2 = rep(NA, 6),
                                y3 = rep(NA, 6),
                                y4 = rep(NA, 6)) # pe = point estimate

container_lqte_d2 <- data.frame(object = c("q1_pe", "q2_pe", "q3_pe", "nonsign_test_stars", "homo_test_stars", "Bandwidth"),
                                y1 = rep(NA, 6),
                                y2 = rep(NA, 6),
                                y3 = rep(NA, 6),
                                y4 = rep(NA, 6)) # pe = point estimate



### D1 -------------

#### Changes of the optimal common bandwidth (selecting second smallest bw) ------------

##### Outcome 1.1: days worked post months [1, 6]
# 1) Estimate bandwidth
bwd1o1_bdyint <- readRDS("intermediate/script02/lqte/bwd1o1_bdyint.RDS")

bwd1o1_min2 <- sort(c(bwd1o1_bdyint$bdy, bwd1o1_bdyint$int))[2]


# 2) Estimate quantiles
## a. Detecting quantile effects == 0
qte_d1o1_05 <- rd.qte(y = outcome1, x = running1, d = d1, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd1o1_min2,
                      tau = seq(0.1, 0.9, by = 0.05))
qte_d1o1_05

qte_d1o1 <- rd.qte(y = outcome1, x = running1, d = d1, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd1o1_min2,
                   tau = seq(0.25, 0.75, by = 0.25))
qte_d1o1


## b. Estimable quantile effects
qte_d1o1_05_est <- rdq.band(y = outcome1, x = running1, d = d1, x0 = 0,
                            cov = 0, bdw = bwd1o1_min2, alpha = 0.1,
                            tau = seq(0.35, 0.85, by = 0.05))


# 3) Run hypothesis tests

Wd1o1_signhomo_nz <- rdq.test(y = outcome1, x = running1, d = d1, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.35, 0.85, by = 0.05),
                              bdw = bwd1o1_min2, bias = 1, type = c(1,2))



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

container_lqte_d1$y1[6] <- round(bwd1o1_min2, 3)

### Outcome 1.2: days worked post months [7, 12]
# 1) Estimate bandwidth
bwd1o2_bdyint <- readRDS("intermediate/script02/lqte/bwd1o2_bdyint.RDS")

bwd1o2_min2 <- sort(c(bwd1o2_bdyint$bdy, bwd1o2_bdyint$int))[2]

# 2) Estimate quantiles
## a. Detecting quantile effects == 0
qte_d1o2_05 <- rd.qte(y = outcome2, x = running1, d = d1, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd1o2_min2,
                      tau = seq(0.1, 0.9, by = 0.05))

qte_d1o2_05

qte_d1o2 <- rd.qte(y = outcome2, x = running1, d = d1, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd1o2_min2,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Non-null quantile effects
qte_d1o2_05_est <- rdq.band(y = outcome2, x = running1, d = d1, x0 = 0,
                            cov = 0, bdw = bwd1o2_min2, alpha = 0.1,
                            tau = seq(0.3, 0.75, by = 0.05))

# 3) Run hypothesis tests
Wd1o2_signhomo_nz <- rdq.test(y = outcome2, x = running1, d = d1, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.3, 0.75, by = 0.05),
                              bdw = bwd1o2_min2, bias = 1, type = c(1,2))


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

container_lqte_d1$y2[6] <- round(bwd1o2_min2, 3)

### Outcome 1.3: days worked post months [13, 18]
# 1) Estimate bandwidth
bwd1o3_bdyint <- readRDS("intermediate/script02/lqte/bwd1o3_bdyint.RDS")

bwd1o3_min2 <- sort(c(bwd1o3_bdyint$bdy, bwd1o3_bdyint$int))[2]


# 2) Estimate quantiles
## a. Detecting quantile effects == 0
qte_d1o3_05 <- rd.qte(y = outcome3, x = running1, d = d1, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd1o3_min2,
                      tau = seq(0.1, 0.9, by = 0.05))

qte_d1o3 <- rd.qte(y = outcome3, x = running1, d = d1, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd1o3_min2,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Non-null quantile effects
qte_d1o3_05_est <- rdq.band(y = outcome3, x = running1, d = d1, x0 = 0,
                            cov = 0, bdw = bwd1o3_min2, alpha = 0.1,
                            tau = seq(0.3, 0.65, by = 0.05))


# 3) Run hypothesis tests
Wd1o3_signhomo_nz <- rdq.test(y = outcome3, x = running1, d = d1, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.3, 0.65, by = 0.05),
                              bdw = bwd1o3_min2, bias = 1, type = c(1,2))



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

container_lqte_d1$y3[6] <- round(bwd1o3_min2, 3)

### Outcome 1.4: days worked post months [19, 24]
# 1) Estimate bandwidth
bwd1o4_bdyint <- readRDS("intermediate/script02/lqte/bwd1o4_bdyint.RDS")

bwd1o4_min2 <- sort(c(bwd1o4_bdyint$bdy, bwd1o4_bdyint$int))[2]


# 2) Estimate quantiles
## a. Detecting quantile effects == 0
qte_d1o4_05 <- rd.qte(y = outcome4, x = running1, d = d1, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd1o4_min2,
                      tau = seq(0.1, 0.9, by = 0.05))

qte_d1o4 <- rd.qte(y = outcome4, x = running1, d = d1, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd1o4_min2,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Non-null quantile effects
qte_d1o4_05_est <- rdq.band(y = outcome4, x = running1, d = d1, x0 = 0,
                            cov = 0, bdw = bwd1o4_min2, alpha = 0.1,
                            tau = seq(0.3, 0.75, by = 0.05))

# 3) Run hypothesis tests
Wd1o4_signhomo_nz <- rdq.test(y = outcome4, x = running1, d = d1, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.3, 0.75, by = 0.05),
                              bdw = bwd1o4_min2, bias = 1, type = c(1,2))

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

container_lqte_d1$y4[6] <- round(bwd1o4_min2, 3)

#### Saving and exporting

if (!file.exists("intermediate/script03/E23a_LQTE_secondbw.xlsx")) {

  openxlsx::write.xlsx(container_lqte_d1, 'intermediate/script03/E23a_LQTE_secondbw.xlsx')

}


### D2 -------------

#### Changes of the optimal common bandwidth (selecting second smallest bw) ------------
### Outcome 1.1: days worked post months [1, 6]
# 1- Estimate bandwidth

bwd2o1_bdyint <- readRDS("intermediate/script02/lqte/bwd2o1_bdyint.RDS")


bwd2o1_min2 <- sort(c(bwd2o1_bdyint$bdy, bwd2o1_bdyint$int))[2]


# 2- Estimate quantiles
## a. Detecting non-zero quantile effects
qte_d2o1_05 <- rd.qte(y = outcome1_d2, x = running2, d = d2, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd2o1_min2,
                      tau = seq(0.1, 0.9, by = 0.05))

qte_d2o1 <- rd.qte(y = outcome1_d2, x = running2, d = d2, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd2o1_min2,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Estimating non-zero quantile effects
qte_d2o1_005b <- rdq.band(y = outcome1_d2, x = running2, d = d2, x0 = 0,
                          cov = 0, bdw = bwd2o1_min2, alpha = 0.1,
                          tau = seq(0.45, 0.9, by = 0.05))

# 3) Hypothesis testing
Wd2o1_signhomo_nz <- rdq.test(y = outcome1_d2, x = running2, d = d2, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.45, 0.9, by = 0.05),
                              bdw = bwd2o1_min2, bias = 1, type = c(1,2))


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

container_lqte_d2$y1[6] <- round(bwd2o1_min2, 3)

### Outcome 1.2: days worked post months [7, 12]
# 1- Estimate bandwidth
bwd2o2_bdyint <- readRDS("intermediate/script02/lqte/bwd2o2_bdyint.RDS")

bwd2o2_min2 <- sort(c(bwd2o2_bdyint$bdy, bwd2o2_bdyint$int))[2]


# 2- Estimate quantiles
## a. Detecting non-zero quantile effects
qte_d2o2_05 <- rd.qte(y = outcome2_d2, x = running2, d = d2, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd2o2_min2,
                      tau = seq(0.1, 0.9, by = 0.05))

qte_d2o2 <- rd.qte(y = outcome2_d2, x = running2, d = d2, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd2o2_min2,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Estimating non-zero quantile effects
quantile(indi_ns_ss2$post_interval712[which(indi_ns_ss2$scoringD2_0 > 0 & indi_ns_ss2$scoringD2_0 < bwd2o2_min)],
         prob=seq(0, 1, length = 21)) # para obtener percentiles en tratados

quantile(indi_ns_ss2$post_interval712[which(indi_ns_ss2$scoringD2_0 > -bwd2o2_min & indi_ns_ss2$scoringD2_0 < 0)],
         prob=seq(0, 1, length = 21)) # para obtener deciles en no-tratados (el percentil 40 es cero)

qte_d2o2_05b <- rdq.band(y = outcome2_d2, x = running2, d = d2, x0 = 0,
                         cov = 0, bdw = bwd2o2_min2, alpha = 0.1,
                         tau = seq(0.45, 0.85, by = 0.05))

# 3- Hypothesis testing
Wd2o2_signhomo_nz <- rdq.test(y = outcome2_d2, x = running2, d = d2, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.45, 0.85, by = 0.05),
                              bdw = bwd2o2_min2, bias = 1, type = c(1,2))



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

container_lqte_d2$y2[6] <- round(bwd2o2_min2, 3)

### Outcome 1.3: days worked post months [13, 18]
# 1- Estimate bandwidth
bwd2o3_bdyint <- readRDS("intermediate/script02/lqte/bwd2o3_bdyint.RDS")

bwd2o3_min2 <- sort(c(bwd2o3_bdyint$bdy, bwd2o3_bdyint$int))[2]


# 2- Estimate quantiles
## a- Detecting non-zero quantile effects
qte_d2o3_05 <- rd.qte(y = outcome3_d2, x = running2, d = d2, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd2o3_min2,
                      tau = seq(0.1, 0.9, by = 0.05))
qte_d2o3_05

qte_d2o3 <- rd.qte(y = outcome3_d2, x = running2, d = d2, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd2o3_min2,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Estimating non-zero quantile effects
quantile(indi_ns_ss2$post_interval1318[which(indi_ns_ss2$scoringD2_0 > 0 & indi_ns_ss2$scoringD2_0 < bwd2o3_min2)],
         prob=seq(0, 1, length = 21)) # para obtener percentiles en tratados (el percentil 40 es cero)

quantile(indi_ns_ss2$post_interval1318[which(indi_ns_ss2$scoringD2_0 > -bwd2o3_min2 & indi_ns_ss2$scoringD2_0 < 0)],
         prob=seq(0, 1, length = 21)) # para obtener deciles en no-tratados

qte_d2o3_05b <- rdq.band(y = outcome3_d2, x = running2, d = d2, x0 = 0,
                         cov = 0, bdw = bwd2o3_min2, alpha = 0.1,
                         tau = seq(0.45, 0.8, by = 0.05))

# 3- Hypothesis testing
Wd2o3_signhomo_nz <- rdq.test(y = outcome3_d2, x = running2, d = d2, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.45, 0.8, by = 0.05),
                              bdw = bwd2o3_min2, bias = 1, type = c(1,2))


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

container_lqte_d2$y3[6] <- round(bwd2o3_min2, 3)

### Outcome 1.4: days worked post months [19, 24]
# 1- Estimate bandwidth

bwd2o4_bdyint <- readRDS("intermediate/script02/lqte/bwd2o4_bdyint.RDS")

bwd2o4_min2 <- sort(c(bwd2o4_bdyint$bdy, bwd2o4_bdyint$int))[2]

# 2- Estimate quantiles
## a- Detecting non-zero quantile effects
qte_d2o4_05 <- rd.qte(y = outcome4_d2, x = running2, d = d2, x0 = 0,
                      cov = 0, bias = 0, bdw = bwd2o4_min2,
                      tau = seq(0.1, 0.9, by = 0.05))

qte_d2o4 <- rd.qte(y = outcome4_d2, x = running2, d = d2, x0 = 0,
                   cov = 0, bias = 0, bdw = bwd2o4_min2,
                   tau = seq(0.25, 0.75, by = 0.25))

## b. Estimating non-zero quantile effects
quantile(indi_ns_ss2$post_interval1924[which(indi_ns_ss2$scoringD2_0 > 0 & indi_ns_ss2$scoringD2_0 < bwd2o4_min2)],
         prob=seq(0, 1, length = 21)) # para obtener percentiles en tratados (el percentil 40 es cero)

quantile(indi_ns_ss2$post_interval1924[which(indi_ns_ss2$scoringD2_0 > -bwd2o4_min2 & indi_ns_ss2$scoringD2_0 < 0)],
         prob=seq(0, 1, length = 21)) # para obtener deciles en no-tratados

qte_d2o4_05b <- rdq.band(y = outcome4_d2, x = running2, d = d2, x0 = 0,
                         cov = 0, bdw = bwd2o4_min2, alpha = 0.1,
                         tau = seq(0.45, 0.8, by = 0.05))

# 3- Hypothesis testing
Wd2o4_signhomo_nz <- rdq.test(y = outcome4_d2, x = running2, d = d2, x0 = 0,
                              cov = 0,
                              alpha = 0.95, tau = seq(0.45, 0.8, by = 0.05),
                              bdw = bwd2o4_min2, bias = 1, type = c(1,2))


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

container_lqte_d2$y4[6] <- round(bwd2o4_min2, 3)

### Saving and exporting

if (!file.exists("intermediate/script03/E23b_LQTE_secondbw.xlsx")) {

  openxlsx::write.xlsx(container_lqte_d2, 'intermediate/script03/E23b_LQTE_secondbw.xlsx')

}


## For LATE_j --------------
### D1 ---------------
#### Binarization -------------

#### Type of contract (binary specifications) in 1-24 time interval

## For D1
oec_d1 <- rdcate(
  score = indi_ns_ss1_c$scoringD1_0,
  outcome = indi_ns_ss1_c$lcris_OEC
)

ftcm12_d1 <- rdcate(
  score = indi_ns_ss1_c$scoringD1_0,
  outcome = indi_ns_ss1_c$lcris_FTCm12
)

ftc612_d1 <- rdcate(
  score = indi_ns_ss1_c$scoringD1_0,
  outcome = indi_ns_ss1_c$lcris_FTC612
)

ftcme6_d1 <- rdcate(
  score = indi_ns_ss1_c$scoringD1_0,
  outcome = indi_ns_ss1_c$lcris_FTCme6
)

table_qualityD1 <- data.frame(
  treatment = c(rep("D1", 2)),
  prOEC = c(
    paste0(oec_d1$pointest, oec_d1$stars),
    oec_d1$ci95
  ),
  prFTCm12 = c(
    paste0(ftcm12_d1$pointest, ftcm12_d1$stars),
    ftcm12_d1$ci95
  ),
  prFTC612 = c(
    paste0(ftc612_d1$pointest, ftc612_d1$stars),
    ftc612_d1$ci95
  ),
  prFTCme6 = c(
    paste0(ftcme6_d1$pointest, ftcme6_d1$stars),
    ftcme6_d1$ci95
  )
)



#### Robustness to alternative outcome variable -----------------

# Multinomial
indi_ns_ss1$workedt1cat <- factor(indi_ns_ss1$workedt1,
                                  levels = c("0", "1"))
indi_ns_ss1$workedt2cat <- factor(indi_ns_ss1$workedt2,
                                  levels = c("0", "1"))
indi_ns_ss1$workedt3cat <- factor(indi_ns_ss1$workedt3,
                                  levels = c("0", "1"))
indi_ns_ss1$workedt4cat <- factor(indi_ns_ss1$workedt4,
                                  levels = c("0", "1"))

rddD1_worked1_cat <- rdcate_multinom(score = indi_ns_ss1$scoringD1_0,
                                     outcome = indi_ns_ss1$workedt1cat)

rddD1_worked2_cat <- rdcate_multinom(score = indi_ns_ss1$scoringD1_0,
                                     outcome = indi_ns_ss1$workedt2cat)

rddD1_worked3_cat <- rdcate_multinom(score = indi_ns_ss1$scoringD1_0,
                                     outcome = indi_ns_ss1$workedt3cat)

rddD1_worked4_cat <- rdcate_multinom(score = indi_ns_ss1$scoringD1_0,
                                     outcome = indi_ns_ss1$workedt4cat)

### D2 ---------------

# OEC
rddD2_p1kT_OEC <- rdrobust(
  y = indi_ns_ss2_c$lcris_OEC, x = indi_ns_ss2_c$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL)

summary(rddD2_p1kT_OEC)


# FTCm12
rddD2_p1kT_FTCm12 <- rdrobust(
  y = indi_ns_ss2_c$lcris_FTCm12, x = indi_ns_ss2_c$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL)

summary(rddD2_p1kT_FTCm12)


# FTC612
rddD2_p1kT_FTC612 <- rdrobust(
  y = indi_ns_ss2_c$lcris_FTC612, x = indi_ns_ss2_c$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL)

summary(rddD2_p1kT_FTC612)

# FTC05
rddD2_p1kT_FTC05 <- rdrobust(
  y = indi_ns_ss2_c$lcris_FTCme6, x = indi_ns_ss2_c$scoringD2_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL)

summary(rddD2_p1kT_FTC05)


#### Binarization ------------


## Type of contract (binary specifications) in 1-24 time interval
oec_d2 <- rdcate(
  score = indi_ns_ss2_c$scoringD2_0,
  outcome = indi_ns_ss2_c$lcris_OEC
)

ftcm12_d2 <- rdcate(
  score = indi_ns_ss2_c$scoringD2_0,
  outcome = indi_ns_ss2_c$lcris_FTCm12
)

ftc612_d2 <- rdcate(
  score = indi_ns_ss2_c$scoringD2_0,
  outcome = indi_ns_ss2_c$lcris_FTC612
)

ftcme6_d2 <- rdcate(
  score = indi_ns_ss2_c$scoringD2_0,
  outcome = indi_ns_ss2_c$lcris_FTCme6
)



table_qualityD2 <- data.frame(
  treatment = c(rep("D2", 2)),
  prOEC = c(
    paste0(oec_d2$pointest, oec_d2$stars),
    oec_d2$ci95
  ),
  prFTCm12 = c(
    paste0(ftcm12_d2$pointest, ftcm12_d2$stars),
    ftcm12_d2$ci95
  ),
  prFTC612 = c(
    paste0(ftc612_d2$pointest, ftc612_d2$stars),
    ftc612_d2$ci95
  ),
  prFTCme6 = c(
    paste0(ftcme6_d2$pointest, ftcme6_d2$stars),
    ftcme6_d2$ci95
  )
)

tablequality <- rbind(table_qualityD1, table_qualityD2)

if (!file.exists("intermediate/script03/E24_LTEj_binarized.xlsx")) {

  openxlsx::write.xlsx(tablequality, 'intermediate/script03/E24_LTEj_binarized.xlsx')

}


#### Robustness to alternative outcome variable -----------------

# Multinomial
indi_ns_ss2$workedt1cat <- factor(indi_ns_ss2$workedt1,
                                  levels = c("0", "1"))
indi_ns_ss2$workedt2cat <- factor(indi_ns_ss2$workedt2,
                                  levels = c("0", "1"))
indi_ns_ss2$workedt3cat <- factor(indi_ns_ss2$workedt3,
                                  levels = c("0", "1"))
indi_ns_ss2$workedt4cat <- factor(indi_ns_ss2$workedt4,
                                  levels = c("0", "1"))

rddD2_worked1_cat <- rdcate_multinom(score = indi_ns_ss2$scoringD2_0,
                                     outcome = indi_ns_ss2$workedt1cat)

rddD2_worked2_cat <- rdcate_multinom(score = indi_ns_ss2$scoringD2_0,
                                     outcome = indi_ns_ss2$workedt2cat)

rddD2_worked3_cat <- rdcate_multinom(score = indi_ns_ss2$scoringD2_0,
                                     outcome = indi_ns_ss2$workedt3cat)

rddD2_worked4_cat <- rdcate_multinom(score = indi_ns_ss2$scoringD2_0,
                                     outcome = indi_ns_ss2$workedt4cat)

rddD2_worked1_cat
rddD2_worked2_cat
rddD2_worked3_cat
rddD2_worked4_cat
