# ...............................................................................
# ASSEGNO PER IL LAVORO - 01 Data wrangling
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
library(rddensity)
library(binsreg)

library(modelsummary)

# 0. Some functions -------
"%out%" <- Negate("%in%")

checknas <- function(data_frame) {
  na_counts <- sapply(data_frame, function(x) sum(is.na(x)))
  return(na_counts)
}

# 1. Reading data ---------
## First delivery (nov. 2022)
data <- read_excel("inputdata/unipd_axl.xlsx",
    sheet = "Dati")

## Second delivery (nov. 2022)
profiling <- read_excel("inputdata/indice profiling continuo.xlsx")

## Third delivery (dec. 2022)
longest <- read_excel("inputdata/axl_unipd_post_max.xlsx")

## Fourth delivery (jan. 2023)
provcomuni <- read_excel("inputdata/provcomuniVEN.xlsx")

ISCOCPitalia <- read_excel("inputdata/raccordo_Isco08_CP2011.xls",
    sheet = "Foglio1", col_names = T)

formadomi <- read_excel("inputdata/unipd_axl estrazione 012023.xlsx",
    sheet = "Foglio1", col_names = T)

disab <- read_excel("inputdata/unipd_axl estrazione 022023.xlsx")

data_notp <- read_excel("inputdata/unipd_axl.xlsx",
    sheet = "Dati")


# 2. Data wrangling ---------

## Sharp RDD? Are there problematic cases? ---------------------------

## Joining *data* (extraction 1) and *profiling* (extraction 2)
data <- inner_join(data, profiling, by = c("ide_soggetto", "axl_data_attribuzione"))

### Differences between the original score columns of extraction 1 and extraction 2 for the same participation (i,t)?
data %>%
  mutate(orig_score_1 = as.numeric(axl_profil_punteggio),
         orig_score_2 = as.numeric(indice)) %>%
  mutate(same_orig_score = if_else(orig_score_1 == orig_score_2, TRUE, FALSE)) %>%
  tabyl(same_orig_score) # This may be indicating that the database recorded *all* the profilings made for each individual at each day, since it seems there are more than 1 profiling by (i,t)

data <- data %>%
  mutate(orig_score_1 = as.numeric(axl_profil_punteggio),
         orig_score_2 = as.numeric(indice)) %>%
  mutate(same_orig_score = if_else(orig_score_1 == orig_score_2, TRUE, FALSE)) %>%
  filter(same_orig_score == T) # Since we don't know which profiling was the last one, we need to remove all cases in which there is more than one profiling for each (i,t)


## Possible problems: (a) discrepancy btw original [used] & estimated [recalculated], or (b) discrepancy btw score and treatment received

### Problem (a)

# Variable "indice" was registered at the moment the individual claimed the AXL.
# Variable "indice_previsto" was estimated at November 2022, possibly with changes at the values of the ingredients of the score function

data$orig_score_2_round <- round_half_up(data$orig_score_2, digits = 2) # implemented score [before called indice_num3]
data$recal_score_round <- round_half_up(data$indice_previsto, digits = 2) # recalculated score as it would be implemented [before called roundhu_3]

data$dif_ass_prev_hu3 <- ifelse(data$orig_score_2_round != data$recal_score_round, TRUE, FALSE) # This marks the problematic cases!

freq(data$dif_ass_prev_hu3)

data_withproblematic <- data
data <- data %>% filter(dif_ass_prev_hu3 == FALSE)

### Problem (b)
data$checking_join <- ifelse(data$axl_profil_fascia == data$intensita, TRUE, FALSE)
freq(data$checking_join) # There are 0 problematic cases!

data_withproblematic %>%
  #filter(dif_ass_prev_hu3 == FALSE) %>%
  mutate(treatment = factor(intensita,
                            levels = c("BASSA", "MEDIA", "ALTA"))) %>%
  mutate(treat = case_when(treatment == "BASSA" ~ "A",
                           treatment == "MEDIA" ~ "B",
                           treatment == "ALTA" ~ "C")) %>%
  ggplot(., aes(x = 1 - recal_score_round, y = treat)) +
  geom_point(size = 0.075) +
  geom_vline(xintercept = c(0.41, 0.61), color = "red") +
  theme_light() +
  labs(y = "Treament received",
       x = "Reversed  score") # saved SVG 395x300


data_withproblematic %>%
  #filter(dif_ass_prev_hu3 == FALSE) %>%
  mutate(treatment = factor(intensita,
                            levels = c("BASSA", "MEDIA", "ALTA"))) %>%
  mutate(treat = case_when(treatment == "BASSA" ~ "A",
                           treatment == "MEDIA" ~ "B",
                           treatment == "ALTA" ~ "C")) %>%
ggplot(., aes(x = 1 - recal_score_round, y = treat)) +
  geom_jitter(size = 0.075) +
  geom_vline(xintercept = c(0.41, 0.61), color = "red") +
  theme_light() +
  labs(y = "Treament received",
       x = "Reversed  score") # saved SVG 395x300


data %>%
  mutate(treatment = factor(intensita,
                            levels = c("BASSA", "MEDIA", "ALTA"))) %>%
ggplot(., aes(x = 1 - as.numeric(axl_profil_punteggio), y = treatment)) +
  geom_jitter() +
  theme_light()


## Declaring proper classes --------------
## Declaring factor variables as factor
data$axl_fascia <- factor(data$axl_profil_fascia, levels = c("BASSA", "MEDIA", "ALTA"))

## Declaring numeric variables as numeric
class(data$indice_previsto) # estimated score
class(data$axl_profil_punteggio) # reported score


# data %>% filter(axl_profil_fascia == "BASSA" & ppa_attivata == "attivata") %>% tabyl(indice_previsto) # group A Rounded with 2 digits:  (0.60 <= s) // Rounded with 3 digits: (0.595 <= s)
# data %>% filter(axl_profil_fascia == "BASSA" & ppa_attivata == "attivata") %>% tabyl(flg_da_eliminare)
# imedia <- data %>% filter(axl_profil_fascia == "MEDIA" & ppa_attivata == "attivata") %>% tabyl(indice_previsto) # group B  Rounded with 2 digits: (0.40 <= s < 0.60) // Rounded with 3 digits: (0.395 <= s < 0.595)
# data %>% filter(axl_profil_fascia == "MEDIA" & ppa_attivata == "attivata") %>% tabyl(flg_da_eliminare)
# ialta <- data %>% filter(axl_profil_fascia == "ALTA" & ppa_attivata == "attivata") %>% tabyl(indice_previsto) # group C Rounded with 2 digits: (s < 0.40) // Rounded with 3 digits: (s < 0.395)
# data %>% filter(axl_profil_fascia == "ALTA" & ppa_attivata == "attivata") %>% tabyl(flg_da_eliminare)

## Declaring date variables as dates
data$axl_data_attribuzione_d <- ymd(data$axl_data_attribuzione)
data$axl_year_attribuzione <- year(data$axl_data_attribuzione_d)

data$ppa_data_avvio_d <- ymd(data$ppa_data_avvio) # Fecha inicio acciones # 5919 failed. 5457 due to non activation of the PPA.
freq(data$ppa_attivata)

data$data_fine_ppa_d <- ymd(data$data_fine_ppa)
data$data_fine_assegno_d <- ymd(data$data_fine_assegno)

data$axl_duration <- data$data_fine_assegno_d - data$ppa_data_avvio_d
data$axl_duration_actions <- ifelse(is.na(data$data_fine_ppa_d), NA,
  data$data_fine_ppa_d - data$ppa_data_avvio_d
)


## Scoring variable in the usual direction ----------
# In the program design, the scoring variable is defined with 3 fractional digits. (A if 0.590 <= s <= 1.000) (B if 0.396 <= s <= 0.589) (C if 0.000 <= s <= 0.395)
# Program implementation was executed with a scoring variable with 2 fractional digits.

### Discrete scoring variable ----------
data_notp$indice_num <- as.numeric(data_notp$axl_profil_punteggio)
data_notp$reversed_indice <- 1 - data_notp$indice_num

data_notp$group <- car::recode(data_notp$axl_profil_fascia, "'BASSA' = 'A'; 'MEDIA' = 'B'; 'ALTA' = 'C'",
  levels = c("A", "B", "C")
)

data_notp %>%
  ggplot(aes(x = reversed_indice, y = group)) +
  geom_point() +
  theme_light() +
  geom_vline(xintercept = 0.41, linewidth = 0.5) +
  geom_vline(xintercept = 0.61, linewidth = 0.5) +
  labs(x = "Scoring variable (S)", y = "Treatment group") # Figure A1 (saved 395 x 300)

rm(data_notp)

### Continuous scoring variable ----------
# freq(data$indice_previsto) # Prob[SÍ estar empleado en un plazo de 24 meses desde el registro].
data$axl_scoring <- 1 - data$indice_previsto
# freq(data$axl_scoring) # Prob[NO estar empleado en un plazo de 24 meses desde el registro]

# group A: Rounded with 3 digits (s < 0.405)
iA <- data %>%
  filter(axl_profil_fascia == "BASSA" & ppa_attivata == "attivata") %>%
  tabyl(axl_scoring)

# group B: Rounded with 2 digits (s >= 0.41  &  s < 0.61) // Rounded with 3 digits: (s >= 0.405 & s < 0.605)
iB <- data %>%
  filter(axl_profil_fascia == "MEDIA" & ppa_attivata == "attivata") %>%
  tabyl(axl_scoring)

# group C: Rounded with 2 digits (s >= 0.61) // Rounded with 3 digits (s >= 0.605)
iC <- data %>%
  filter(axl_profil_fascia == "ALTA" & ppa_attivata == "attivata") %>%
  tabyl(axl_scoring)

# And centered at cero...
data$scoringD1_0 <- data$axl_scoring - 0.405 # para el tratamiento 1 (de A a B)
data$scoringD2_0 <- data$axl_scoring - 0.605 # para el tratamiento 2 (de B a C)

## Outcome variables ------------
### Quantitative outcomes --------------
#### Worked days by semester ------------
data$post_interval6 <- data$m1_gg_lav + data$m2_gg_lav + data$m3_gg_lav + data$m4_gg_lav + data$m5_gg_lav + data$m6_gg_lav

data$post_interval712 <- data$m7_gg_lav + data$m8_gg_lav + data$m9_gg_lav + data$m10_gg_lav + data$m11_gg_lav + data$m12_gg_lav

data$post_interval12 <- data$m1_gg_lav + data$m2_gg_lav + data$m3_gg_lav + data$m4_gg_lav + data$m5_gg_lav + data$m6_gg_lav +
  data$m7_gg_lav + data$m8_gg_lav + data$m9_gg_lav + data$m10_gg_lav + data$m11_gg_lav + data$m12_gg_lav

data$post_interval1318 <- data$m13_gg_lav + data$m14_gg_lav + data$m15_gg_lav + data$m16_gg_lav + data$m17_gg_lav + data$m18_gg_lav

data$post_interval18 <- data$m1_gg_lav + data$m2_gg_lav + data$m3_gg_lav + data$m4_gg_lav + data$m5_gg_lav + data$m6_gg_lav +
  data$m7_gg_lav + data$m8_gg_lav + data$m9_gg_lav + data$m10_gg_lav + data$m11_gg_lav + data$m12_gg_lav +
  data$m13_gg_lav + data$m14_gg_lav + data$m15_gg_lav + data$m16_gg_lav + data$m17_gg_lav + data$m18_gg_lav

data$post_interval24 <- data$m1_gg_lav + data$m2_gg_lav + data$m3_gg_lav + data$m4_gg_lav + data$m5_gg_lav + data$m6_gg_lav +
  data$m7_gg_lav + data$m8_gg_lav + data$m9_gg_lav + data$m10_gg_lav + data$m11_gg_lav + data$m12_gg_lav +
  data$m13_gg_lav + data$m14_gg_lav + data$m15_gg_lav + data$m16_gg_lav + data$m17_gg_lav + data$m18_gg_lav +
  data$m19_gg_lav + data$m20_gg_lav + data$m21_gg_lav + data$m22_gg_lav + data$m23_gg_lav + data$m24_gg_lav

data$post_interval1924 <- data$m19_gg_lav + data$m20_gg_lav + data$m21_gg_lav + data$m22_gg_lav + data$m23_gg_lav + data$m24_gg_lav

data$jshours <- data$attiv_indiv_ore_eff


## Costs (for the training cost, we only need to consider the cost of the first AXL, so we need to use the dataset at individual level)
data$ente_ris_occ_importo_r <- car::recode(data$ente_ris_occ_importo, '"no risultato occupaz" = "0"',
  as.factor = F
)

data$spent_incentive <- as.numeric(data$ente_ris_occ_importo_r)
data$spent_js <- data$jshours * 27
data$spent_incentivejs <- data$spent_incentive + data$spent_js


data$axlduration <- as.numeric(data$axl_duration)
data$axldurationactions <- as.numeric(data$axl_duration_actions)


#### Total days in unemployment --------------

data <- data %>%
  mutate(pre_oss_ultimo_data_cessazioneN = if_else(pre_oss_ultimo_data_cessazione == "0000-00-00", NA,
                                                  pre_oss_ultimo_data_cessazione),
         pre_no_oss_ultimo_data_cessazioneN = if_else(pre_no_oss_ultimo_data_cessazione == "0000-00-00", NA,
                                                      pre_no_oss_ultimo_data_cessazione),
         post_oss_primo_data_inizio_rapportoN = if_else(post_oss_primo_data_inizio_rapporto == "0000-00-00", NA,
                                                      post_oss_primo_data_inizio_rapporto),
         post_no_oss_primo_data_inizio_rapportoN = if_else(post_no_oss_primo_data_inizio_rapporto == "0000-00-00", NA,
                                                           post_no_oss_primo_data_inizio_rapporto)) %>%

  mutate(pre_oss_ultimo_data_cessazioneD = ymd(pre_oss_ultimo_data_cessazioneN),
         pre_no_oss_ultimo_data_cessazioneD = ymd(pre_no_oss_ultimo_data_cessazioneN),
         post_oss_primo_data_inizio_rapportoD = ymd(post_oss_primo_data_inizio_rapportoN),
         post_no_oss_primo_data_inizio_rapportoD = ymd(post_no_oss_primo_data_inizio_rapportoN)) %>%

  mutate(pre_ultimo_data_cessazioneD = case_when(is.na(pre_oss_ultimo_data_cessazioneD) & is.na(pre_no_oss_ultimo_data_cessazioneD) ~ ymd(axl_data_attribuzione) - months(36),
                                                 is.na(pre_oss_ultimo_data_cessazioneD) ~ pre_no_oss_ultimo_data_cessazioneD,
                                                 is.na(pre_no_oss_ultimo_data_cessazioneD) ~ pre_oss_ultimo_data_cessazioneD,
                                                 pre_oss_ultimo_data_cessazioneD > pre_no_oss_ultimo_data_cessazioneD ~ pre_oss_ultimo_data_cessazioneD,
                                                 TRUE ~ pre_no_oss_ultimo_data_cessazioneD),
         post_primo_data_inizioD = case_when(is.na(post_oss_primo_data_inizio_rapportoD) & is.na(post_no_oss_primo_data_inizio_rapportoD) ~ ymd(axl_data_attribuzione) + months(24),
                                             is.na(post_oss_primo_data_inizio_rapportoD) ~ post_no_oss_primo_data_inizio_rapportoD,
                                             is.na(post_no_oss_primo_data_inizio_rapportoD) ~ post_oss_primo_data_inizio_rapportoD,
                                             post_oss_primo_data_inizio_rapportoD < post_no_oss_primo_data_inizio_rapportoD ~ post_oss_primo_data_inizio_rapportoD,
                                             TRUE ~ post_no_oss_primo_data_inizio_rapportoD)) %>%

  mutate(unemplength = interval(pre_ultimo_data_cessazioneD, post_primo_data_inizioD) / ddays(1),
         warning_unempl = if_else(post_primo_data_inizioD < pre_ultimo_data_cessazioneD | year(pre_ultimo_data_cessazioneD) == "2099", T, F))

#data2 <- data %>%
#  filter(year(pre_ultimo_data_cessazioneD) != "2099") %>%
#  filter(post_primo_data_inizioD > pre_ultimo_data_cessazioneD)


totalfreq <- data %>% tabyl(unemplength) %>% mutate(cumfreq = cumsum(percent))

datafreqpre <- data %>% tabyl(pre_ultimo_data_cessazioneD)
datafreqpost <- data %>% tabyl(post_primo_data_inizioD)

### Qualitative outcomes ----------

#### Longest contract (expected duration) --------------
class(longest$post_oss_max_durata)
longest$post_oss_max_durata_n <- as.numeric(longest$post_oss_max_durata)

class(longest$post_no_oss_max_durata)
longest$post_no_oss_max_durata_n <- as.numeric(longest$post_no_oss_max_durata)

#### This gets the type of the longest labour contract, considering both collected and not collected by the Observatory
longest <- longest %>% mutate(longestcpost = case_when(
  post_oss_max_durata_n > post_no_oss_max_durata_n ~ post_oss_max_contratto2_ass,
  post_oss_max_durata_n < post_no_oss_max_durata_n ~ post_no_oss_max_contratto2_ass,
  post_oss_max_durata_n == post_no_oss_max_durata_n ~ post_oss_max_contratto2_ass,
  TRUE ~ post_no_oss_max_contratto2_ass
))

freq(longest$longestcpost)
sum(is.na(longest$post_oss_max_contratto2_ass) & is.na(longest$post_no_oss_max_contratto2_ass)) # The number of NAs is ok

#### This gets the duration of the longest labour contract, considering both collected and not collected by the Observatory
longest <- longest %>% mutate(durata_longestcpost = case_when(
  post_oss_max_durata_n > post_no_oss_max_durata_n ~ post_oss_max_durata_n,
  post_oss_max_durata_n < post_no_oss_max_durata_n ~ post_no_oss_max_durata_n,
  post_oss_max_durata_n == post_no_oss_max_durata_n ~ post_oss_max_durata_n,
  TRUE ~ post_no_oss_max_durata_n
))

longest$durata_longestcpost_ld <- days(longest$durata_longestcpost)
longest$durata_longestcpost_ld_months <- longest$durata_longestcpost_ld %/% months(1)

# Option A. Italian classification of labour contracts
longest <- longest %>%
  mutate(longestcontract = case_when(
    longestcpost %in% c("b- Cap", "f- Dom", "g- Par") ~ "AC", # Alternative Contract
    longestcpost == "a- Cti" ~ "OEC", # Open-Ended Contract
    longestcpost == "c- Ctd" ~ "FTC", # Fixed-Term Contract
    longestcpost == "d- Som" ~ "SC", # temporary agency work (called "contratto di somministrazione di lavoro")
    longestcpost == "e- Int" ~ "DC", # Discontinuous Contract (called "contratto di lavoro intermittente"). In Italy, it may be open-ended or temporary.
    is.na(longestcpost) ~ "NE" # Not Employed (there isn't labour contract)
  ))
freq(longest$longestcontract)

# Option B. Standard classification of labour contracts

## longestcontract_ris registers the longest contract signed in the [+1, +24] time interval after *starting* the treatment
longest <- longest %>% mutate(
  longestcontract_standard = case_when(
    durata_longestcpost_ld == "10000d 0H 0M 0S" ~ "OEC", # this includes contracts that are not a-Cti but have open-ended duration
    durata_longestcpost_ld_months > 12 ~ "FTCm12",
    durata_longestcpost_ld_months > 5 & durata_longestcpost_ld_months < 13 ~ "FTC612",
    durata_longestcpost_ld_months < 6 & durata_longestcpost_ld > 0 ~ "FTCme6",
    durata_longestcpost_ld_months == 0 ~ "NE"
  ),
  longestcontract_ris = case_when(
    durata_longestcpost_ld == "10000d 0H 0M 0S" ~ "OEC",
    durata_longestcpost_ld_months > 12 & longestcpost %in% c("b- Cap", "c- Ctd", "d- Som") ~ "FTCm12", # this limits FTC contracts for these 3 types of contract
    durata_longestcpost_ld_months > 5 & durata_longestcpost_ld_months < 13 & longestcpost %in% c("b- Cap", "c- Ctd", "d- Som") ~ "FTC612",
    durata_longestcpost_ld_months < 6 & durata_longestcpost_ld > 0 & longestcpost %in% c("b- Cap", "c- Ctd", "d- Som") ~ "FTCme6",
    TRUE ~ "Not_employed_or_other_contract"
  )
)

#### Awarded contract (checked duration) -----------------
data <- data %>%
  mutate(prizetype = case_when(ente_ris_occ_contratto == "no risultato occupaz" ~ "prize0",
                               ente_ris_occ_contratto == "TD6" ~ "prize1",
                               ente_ris_occ_contratto == "TD12" ~ "prize2",
                               ente_ris_occ_contratto == "TI" ~ "prize3"))


## Cleaning variables on disability and occupations ---------
disab$ide_data <- paste0(disab$ide_soggetto, "_", disab$axl_dat_attribuzione)

## Occupations
colnames(ISCOCPitalia) <- c("ISCO08_code", "ISCO08_name", "CP2011_code", "CP2011_name")
ISCOCPitalia <- ISCOCPitalia[-c(1), ]
ISCOCPitalia$ISCO08_1digit <- str_sub(ISCOCPitalia$ISCO08_code, 1, 1)
ISCOCPitalia$CP2011_1digit <- str_sub(ISCOCPitalia$CP2011_code, 1, 1)



## Dataset at individual level ------------

### Removing non activated PPA
freq(data$ppa_attivata)
axlstarted <- data %>% filter(ppa_attivata == "attivata")

### Without removing those individuals treated more than once (but considering only its first treatment)
indi <- axlstarted %>%
  group_by(ide_soggetto) %>%
  filter(axl_data_attribuzione_d == min(axl_data_attribuzione_d))


### Removing suspicious cases --------------
## 1. Suspicious A = actual hours of treatment > maximum hours treatment by design
indi <- indi %>% mutate(suspicious = case_when(
  axl_fascia == "BASSA" & attiv_indiv_ore_eff > 7 ~ "Yes",
  axl_fascia == "MEDIA" & attiv_indiv_ore_eff > 13 ~ "Yes",
  axl_fascia == "ALTA" & attiv_indiv_ore_eff > 27 ~ "Yes",
  TRUE ~ "No"
)) # The number is the maximum hours of job search measures + 1 (according to the design)
freq(indi$suspicious)

## 2. Suspicious B = AXL duration < 0
indi <- indi %>% mutate(negativeaxlduration = case_when(
  axlduration < 0 ~ "Yes",
  axlduration == 0 ~ "No",
  TRUE ~ "No"
))
freq(indi$negativeaxlduration)

# problematicos_axl_levl <- data %>% filter(ide_soggetto %in% problematicos$ide_soggetto) # 21 observations vs. 19 individuals. Only 2 individuals have been treated more than once: 3124300 and 3451008 (this with AXL "retired" the second time)

summary <- indi %>%
  filter(suspicious == "No" | negativeaxlduration == "No") %>%
  group_by(axl_fascia) %>%
  summarise(
    jobsearch_h_m = mean(attiv_indiv_ore_eff),
    jobsearch_h_sd = sd(attiv_indiv_ore_eff),
    training_h_m = mean(attiv_form_ore_prev),
    training_h_sd = sd(attiv_form_ore_prev),
    jobsearch_med = median(attiv_indiv_ore_eff),
    training_med = median(attiv_form_ore_prev),
    employed_during = 1 - mean(ente_ris_occ_contratto == "no risultato occupaz")
  )

summary <- adorn_rounding(summary, digits = 2)

# openxlsx::write.xlsx(summary, 'C:/Users/alvar/UAB/OneDrive - Universitat Autònoma de Barcelona/PhD thesis/00A_thesis/3_mix_Veneto/Intermediate_outputs/summary_actualtr2201.xlsx')

# 3. Filtering
indi_ns <- indi %>%
  filter(suspicious == "No") %>%
  filter(negativeaxlduration == "No")
nrow(indi) - nrow(indi_ns) # 39 cases were removed

indi_ns$ide_data <- paste0(indi_ns$ide_soggetto, "_", indi_ns$axl_data_attribuzione)

### Pretreatment covariates ---------------
### Maximum educational attainment
indi_ns <- indi_ns %>% mutate(studio2_grouped = case_when(
  studio2 == "Licenza elementare" ~ "1_primary",
  studio2 == "Licenza media" ~ "2_lowersec",
  studio2 == "Diploma di qualifica professionale" ~ "353_uppersecP_3y",
  studio2 == "Diploma Tecnico" ~ "353_uppersecP_4y_tec",
  studio2 == "Diploma Professionale" ~ "353_uppersecP_4y_prof",
  studio2 == "Diploma Liceale" ~ "344_uppersecG",
  studio2 == "Laurea I livello (triennale)" | studio3 == "Diploma universitario" ~ "660_terBach_lev1", # Art. 17 legge 30 dicembre 2010, n. 240.
  studio2 %in% c("Laurea - vecchio o nuovo ordinamento", "Post laurea") ~ "760+_terMastery+_lev2",
  studio3 %in% c(
    "Altri titoli di istruzione superiore", "Diploma Conservatorio musicale",
    "Maestro d'arte", "Scuola magistrale (triennale)", "Diploma di istruzione artistica", "Diploma interprete, traduttore, archivista"
  ) ~ "6or7_artshhss",
  TRUE ~ "Nondetermined"
)) # the first number is the ISCED-A

indi_ns$studio2_grouped2 <- car::recode(indi_ns$studio2_grouped,
  'c("660_terBach_lev1", "6or7_artshhss") = "660_terBachy+_lev1"; "Nondetermined" = NA',
  as.factor = T
)

### Lagged outcome variables (number of working days)
indi_ns$prewd1_6 <- indi_ns$m_meno1_gg_lav + indi_ns$m_meno2_gg_lav + indi_ns$m_meno3_gg_lav + indi_ns$m_meno4_gg_lav +
  indi_ns$m_meno5_gg_lav + indi_ns$m_meno6_gg_lav

indi_ns$prewd7_12 <- indi_ns$m_meno7_gg_lav + indi_ns$m_meno8_gg_lav + indi_ns$m_meno9_gg_lav +
  indi_ns$m_meno10_gg_lav + indi_ns$m_meno11_gg_lav + indi_ns$m_meno12_gg_lav

indi_ns$prewd13_18 <- indi_ns$m_meno13_gg_lav + indi_ns$m_meno14_gg_lav + indi_ns$m_meno15_gg_lav + indi_ns$m_meno16_gg_lav +
  indi_ns$m_meno17_gg_lav + indi_ns$m_meno18_gg_lav

indi_ns$prewd19_24 <- indi_ns$m_meno19_gg_lav + indi_ns$m_meno20_gg_lav + indi_ns$m_meno21_gg_lav +
  indi_ns$m_meno22_gg_lav + indi_ns$m_meno23_gg_lav + indi_ns$m_meno24_gg_lav

indi_ns$prewd1_12 <- indi_ns$m_meno1_gg_lav + indi_ns$m_meno2_gg_lav + indi_ns$m_meno3_gg_lav + indi_ns$m_meno4_gg_lav +
  indi_ns$m_meno5_gg_lav + indi_ns$m_meno6_gg_lav + indi_ns$m_meno7_gg_lav + indi_ns$m_meno8_gg_lav + indi_ns$m_meno9_gg_lav +
  indi_ns$m_meno10_gg_lav + indi_ns$m_meno11_gg_lav + indi_ns$m_meno12_gg_lav

indi_ns$prewd13_24 <- indi_ns$m_meno13_gg_lav + indi_ns$m_meno14_gg_lav + indi_ns$m_meno15_gg_lav + indi_ns$m_meno16_gg_lav +
  indi_ns$m_meno17_gg_lav + indi_ns$m_meno18_gg_lav + indi_ns$m_meno19_gg_lav + indi_ns$m_meno20_gg_lav + indi_ns$m_meno21_gg_lav +
  indi_ns$m_meno22_gg_lav + indi_ns$m_meno23_gg_lav + indi_ns$m_meno24_gg_lav

indi_ns$prewd25_36 <- indi_ns$m_meno25_gg_lav + indi_ns$m_meno26_gg_lav + indi_ns$m_meno27_gg_lav + indi_ns$m_meno28_gg_lav +
  indi_ns$m_meno29_gg_lav + indi_ns$m_meno30_gg_lav + indi_ns$m_meno31_gg_lav + indi_ns$m_meno32_gg_lav + indi_ns$m_meno33_gg_lav +
  indi_ns$m_meno34_gg_lav + indi_ns$m_meno35_gg_lav + indi_ns$m_meno36_gg_lav

### Duration of the last unemployment spell (in number of days)
indi_ns$pre_oss_ultimo_data_cessazione_d <- ymd(indi_ns$pre_oss_ultimo_data_cessazione)
indi_ns$pre_no_oss_ultimo_data_cessazione_d <- ymd(indi_ns$pre_no_oss_ultimo_data_cessazione)

indi_ns$last_nempl_spell_oss <- indi_ns$axl_data_attribuzione_d - indi_ns$pre_oss_ultimo_data_cessazione_d
indi_ns$last_nempl_spell_no_oss <- indi_ns$axl_data_attribuzione_d - indi_ns$pre_no_oss_ultimo_data_cessazione_d # 143 cases with incoherent values (negatives)

proble_ultimoepi <- subset(indi_ns, indi_ns$last_nempl_spell_no_oss < 0)
proble_ultimoepi <- proble_ultimoepi[, c("axl_data_attribuzione_d", "pre_no_oss_ultimo_data_cessazione_d")] # observing these incoherent values, better to consider them as NAs
rm(proble_ultimoepi)

indi_ns$last_nempl_spell_no_oss_ok <- ifelse(indi_ns$last_nempl_spell_no_oss < 0, NA, indi_ns$last_nempl_spell_no_oss) # class numeric
summary(indi_ns$last_nempl_spell_no_oss_ok)

indi_ns$last_nempl_spell_oss_ok <- as.numeric(indi_ns$last_nempl_spell_oss)

indi_ns$last_nempl_spell <- pmin(indi_ns$last_nempl_spell_oss_ok, indi_ns$last_nempl_spell_no_oss_ok, na.rm = T)

### Type of the last contract
indi_ns$group_last_contract <- ifelse(indi_ns$last_nempl_spell == indi_ns$last_nempl_spell_oss_ok, "standard", "non-standard")

indi_ns$type_last_contract <- ifelse(indi_ns$group_last_contract == "standard", indi_ns$pre_oss_ultimo_contratto2_ass, indi_ns$pre_no_oss_ultimo_contratto2_ass)

### Sector of the last contract
indi_ns$sector_last_contract <- ifelse(indi_ns$group_last_contract == "standard", indi_ns$pre_oss_ultimo_cod_ateco, indi_ns$pre_no_oss_ultimo_cod_ateco)

indi_ns$sector2_last_contract <- str_sub(indi_ns$sector_last_contract, start = 1, end = 2)
indi_ns$sector3_last_contract <- paste0(indi_ns$sector2_last_contract, str_sub(indi_ns$sector_last_contract, start = 4, end = 4))
indi_ns$sector4_last_contract <- paste0(indi_ns$sector2_last_contract, str_sub(indi_ns$sector_last_contract, start = 4, end = 5))
head(indi_ns$sector4_last_contract)

indi_ns$sector5_last_contract <- paste0(indi_ns$sector4_last_contract, str_sub(indi_ns$sector_last_contract, start = 7, end = 7))
head(indi_ns$sector5_last_contract)

indi_ns <- indi_ns %>% mutate(sectorVL = case_when(
  sector2_last_contract %in% c("10", "11", "12", "13", "14", "16", "31") |
    sector3_last_contract %in% c("151", "152", "231", "234", "237", "321", "322", "323", "324") |
    sector5_last_contract %in% c("32504", "32505") ~ "MadeInItaly",
  sector2_last_contract %in% c("17", "18", "19", "20", "22", "21", "35", "36", "37", "38", "39", "05", "06", "07", "08", "09") |
    sector5_last_contract %in% c("32501", "32502", "32503") |
    sector3_last_contract %in% c("232", "233", "235", "236", "239", "329") ~ "AltreIndustrie", # desde el 35 pertenece a utilities; también añadidas las industria extractivas ("05", "06", "07", "08", "09")
  sector2_last_contract %in% c("41", "42", "43") ~ "Costruzione",
  sector2_last_contract %in% c("47", "55", "56", "79", "90", "91", "92", "93") |
    sector3_last_contract %in% c("823") | sector4_last_contract %in% c("4532", "9604") ~ "CommercioTempoLib",
  sector2_last_contract %in% c("84", "85", "75", "86", "87", "88", "97", "95", "94", "96", "98", "99") |
    sector3_last_contract %in% c("452", "772") | sector5_last_contract %in% c("45403") ~ "ServiziPersona",
  sector2_last_contract %in% c("46", "49", "50", "51", "52", "53") | sector3_last_contract %in% c("451", "453", "454") ~ "IngrossoLogistica",
  sector2_last_contract %in% c("58", "59", "60", "61", "62", "63", "69", "70", "71", "73", "74", "72") |
    sector3_last_contract %in% c("774", "781") | sector4_last_contract %in% c("6391") ~ "TerziarioAvanzato",
  sector2_last_contract %in% c("80", "82", "77", "68", "64", "65", "66") |
    sector3_last_contract %in% c("812", "811", "813", "782", "783") ~ "AltriServizi", # ultimas dos de agenzie somministrazione, también añadidos los de servicios financieros ("64", "65", "66")
  sector2_last_contract %in% c("24", "25", "28", "33", "26", "27", "29", "30") ~ "MetalMeccanico",
  is.na(sector2_last_contract) ~ "Unobserved",
  sector2_last_contract %in% c("01", "02", "03") ~ "Agricoltura"
)) # not included in the profiling model


### Occupation of the last contract
indi_ns$edu_last_contract <- ifelse(indi_ns$group_last_contract == "standard", indi_ns$pre_oss_ultimo_cod_qualifica, indi_ns$pre_no_oss_ultimo_cod_qualifica)

indi_ns$edu1_last_contract <- str_sub(indi_ns$edu_last_contract, 1, 1)
freq(indi_ns$edu1_last_contract)

indi_ns <- indi_ns %>% mutate(edulastcontract = case_when(
  edu1_last_contract %in% c("1", "2") ~ "Intellettuali",
  edu1_last_contract == "3" ~ "Tecniche",
  edu1_last_contract == "5" ~ "Qualiservizi",
  edu1_last_contract == "6" ~ "OpeSpecializzati",
  edu1_last_contract == "7" ~ "Semispecializzati",
  edu1_last_contract == "8" ~ "Nonqualificate",
  is.na(edu1_last_contract) ~ "Nonobserved",
  edu1_last_contract == "4" ~ "Impiegati",
  TRUE ~ "NA"
))
indi_ns$edu_lastcont <- car::recode(indi_ns$edulastcontract, '"NA" = NA', as.factor = T)

# From CP2011 (Italian classification) to ISCO-08
indi_ns <- left_join(indi_ns, ISCOCPitalia[, c("CP2011_code", "ISCO08_code", "ISCO08_1digit")],
  by = c("edu_last_contract" = "CP2011_code")
)


### Duration of the last employment spell (in number of days)
indi_ns$duration_last_es <- ifelse(indi_ns$group_last_contract == "standard", indi_ns$pre_oss_ultimo_durata, indi_ns$pre_no_oss_ultimo_durata)


### Reviewing all candidate control variables
indi_ns$last_nempl_spell_ok <- ifelse(is.na(indi_ns$last_nempl_spell), 1095, indi_ns$last_nempl_spell)

indi_ns$type_last_contract_ok <- ifelse(is.na(indi_ns$type_last_contract), "Unobserved", indi_ns$type_last_contract)

indi_ns$sector2_last_contract_ok <- ifelse(is.na(indi_ns$sector2_last_contract), "Unobserved", indi_ns$sector2_last_contract)

indi_ns$sector2_last_contract_num <- as.numeric(indi_ns$sector2_last_contract_ok)
indi_ns <- indi_ns %>% mutate(sector1_last_contract_ok = case_when(
  sector2_last_contract_num %in% c(1:3) ~ "A",
  sector2_last_contract_num %in% c(5:9) ~ "B",
  sector2_last_contract_num %in% c(10:33) ~ "C",
  sector2_last_contract_num %in% c(35) ~ "D",
  sector2_last_contract_num %in% c(36:39) ~ "E",
  sector2_last_contract_num %in% c(41:43) ~ "F",
  sector2_last_contract_num %in% c(45:47) ~ "G",
  sector2_last_contract_num %in% c(49:53) ~ "H",
  sector2_last_contract_num %in% c(55:56) ~ "I",
  sector2_last_contract_num %in% c(58:63) ~ "J",
  sector2_last_contract_num %in% c(64:66) ~ "K",
  sector2_last_contract_num %in% c(68) ~ "L",
  sector2_last_contract_num %in% c(69:75) ~ "M",
  sector2_last_contract_num %in% c(77:82) ~ "N",
  sector2_last_contract_num %in% c(84) ~ "O",
  sector2_last_contract_num %in% c(85) ~ "P",
  sector2_last_contract_num %in% c(86:88) ~ "Q",
  sector2_last_contract_num %in% c(90:93) ~ "R",
  sector2_last_contract_num %in% c(94:96) ~ "S",
  sector2_last_contract_num %in% c(97:98) ~ "T",
  sector2_last_contract_num %in% c(99) ~ "U",
  TRUE ~ "Unobserved"
))


indi_ns$duration_last_es36_ok <- ifelse(is.na(indi_ns$duration_last_es), 0, indi_ns$duration_last_es)

## Exporting datasets to use them in Stata --------------
indi_ns_stata <- indi_ns %>%
  select(c(
    "genere", "eta", "studio3", "studio2", "flg_italiano", "prewd1_12", "prewd13_24", "prewd25_36", "last_nempl_spell_ok", "sector1_last_contract_ok",
    "type_last_contract_ok", "sector2_last_contract_ok", "edu_lastcont", "duration_last_es36_ok", "pre_oss_max_durata", "pre_no_oss_max_durata", "ISCO08_code",
    "ide_soggetto", "axl_data_attribuzione", "scoringD1_0", "scoringD2_0", "post_interval6", "post_interval12", "post_interval18", "post_interval24", "prizetype"))



indi_ns_stata6 <- indi_ns[, c(
  "genere", "eta", "studio2", "studio3", "flg_italiano", "prewd1_12", "prewd13_24", "prewd25_36", "last_nempl_spell_ok",
  "prewd1_6", "prewd7_12", "prewd13_18", "prewd19_24", "edu_lastcont", "sector1_last_contract_ok", "ISCO08_code",
  "type_last_contract_ok", "sector2_last_contract_ok", "duration_last_es36_ok", "pre_oss_max_durata", "pre_no_oss_max_durata",
  "ide_soggetto", "axl_data_attribuzione", "scoringD1_0", "scoringD2_0", "post_interval6", "post_interval12", "post_interval18", "post_interval24",
  "post_interval712", "post_interval1318", "post_interval1924", "jshours", "attiv_form_ore_prev", "prizetype",
  "indice_previsto", "axl_scoring", "ppa_data_avvio_d",
  "axl_profil_punteggio", "axl_profil_fascia", "spent_incentivejs", "spent_incentive", "spent_js",
  "studio2_grouped", "studio2_grouped2", "sectorVL", "data_fine_ppa_d",
  "m_meno1_gg_lav", "m_meno2_gg_lav", "m_meno3_gg_lav", "m_meno4_gg_lav", "m_meno5_gg_lav", "m_meno6_gg_lav",
  "axl_ente", "axl_sportello", "cpi_titolare"
)] # needed for heterogeneity effects

indi_ns_stataD1 <- subset(indi_ns_stata6, indi_ns_stata6$scoringD2_0 < 0) # This includes group A and group B
indi_ns_stataD2 <- subset(indi_ns_stata6, indi_ns_stata6$scoringD1_0 > 0) # This includes group B and group C

indi_ns_stataD1$treatedD1 <- ifelse(indi_ns_stataD1$scoringD1_0 > 0, "Treated", "Nontreated")

#haven::write_dta(indi_ns_stataD1, "intermediate/script01/indi_ns_stata_D1.dta")
#haven::write_dta(indi_ns_stataD2, "intermediate/script01/indi_ns_stata_D2.dta")


## Dataset on training and occupations -----------------
sum(duplicated(formadomi$ide_soggetto))

longest$post_oss_max_durata_n <- as.numeric(longest$post_oss_max_durata)

occupations <- longest %>%
  mutate(ide_data = paste(ide_soggetto, axl_data_attribuzione, sep = "_")) %>%
  mutate(occ_longeste = case_when(
    post_oss_max_durata_n > post_no_oss_max_durata_n ~ post_oss_max_cod_qualifica,
    post_oss_max_durata_n < post_no_oss_max_durata_n ~ post_no_oss_max_cod_qualifica,
    post_oss_max_durata_n == post_no_oss_max_durata_n ~ post_oss_max_cod_qualifica,
    TRUE ~ post_no_oss_max_cod_qualifica)) %>%
  select(ide_data, occ_longeste)

# 3. Summary statistics -----------------------
### Table of summary statistics by group and total (Table 1)
indi_ns$trgroup <- car::recode(indi_ns$axl_profil_fascia,
  '"ALTA" = "C"; "MEDIA" = "B"; "BASSA" = "A"',
  as.factor = T,
  levels = c("A", "B", "C")
)

# level of study
indi_ns$studio2_grouped3 <- car::recode(indi_ns$studio2_grouped2,
  'c("353_uppersecP_4y_prof", "353_uppersecP_4y_tec") = "353_uppersecP_4y"',
  as.factor = T
)

indi_ns <- indi_ns %>%
  mutate(study_agg = case_when(studio2_grouped3 == "1_primary" ~ "Primary",
                               studio2_grouped3 == "2_lowersec" ~ "Lower secondary",
                               studio2_grouped3 %in% c("344_uppersecG") ~ "Upper secondary general",
                               studio2_grouped3 %in% c("353_uppersecP_3y", "353_uppersecP_4y") ~ "Upper secondary vocational",
                               studio2_grouped3 %in% c("660_terBachy+_lev1", "760+_terMastery+_lev2") ~ "University or equivalent (arts included)",
                               is.na(studio2_grouped3) ~ "Unobserved"))


f_study <- indi_ns %>%
  tabyl(study_agg, trgroup, show_na = F) %>%
  adorn_percentages(denominator = "col", na.rm = T) %>%
  adorn_pct_formatting(affix_sign = F) %>%
  slice(2, 1, 5, 6, 3, 4)

f_study_all <- indi_ns %>%
  tabyl(study_agg) %>%
  adorn_percentages(denominator = "col") %>%
  adorn_pct_formatting(affix_sign = F) %>%
  slice(2, 1, 5, 6, 3, 4)

freqstudy <- cbind(f_study, all = f_study_all[, 3])
colnames(freqstudy)[5] <- "all"

# occupation
indi_ns$ISCO08_code1 <- str_sub(indi_ns$ISCO08_code, 1, 1)

indi_ns$ISCO08_code1r <- car::recode(indi_ns$ISCO08_code1,
  '"0" = "0 Armed forces"; "1" = "1 Managers";
                                     "2" = "2 Professionals"; "3" = "3 Technicians and associate professionals";
                                     "4" = "4 Clerical support workers"; "5" = "5 Services and sales workers";
                                     "6" = "6 Skilled agricultural, forestry and fishery workers";
                                     "7" = "7 Craft and related trades workers";
                                     "8" = "8 Plant and machine operators, and assemblers";
                                     "9" = "9 Elementary occupations";
                                     NA = "X Unobserved"',
  as.factor = T
)

indi_ns <- indi_ns %>%
  mutate(ISCO08_agg = case_when(ISCO08_code1 %in% c("1", "2", "3") ~ "High skilled white collar (ISCO 1, 2, and 3)",
                                ISCO08_code1 %in% c("4", "5") ~ "Low skilled white collar (ISCO 4 and 5)",
                                ISCO08_code1 %in% c("6", "7") ~ "High skilled blue collar (ISCO 6 and 7)",
                                ISCO08_code1 %in% c("8", "9") ~ "Low skilled blue collar (ISCO 8 and 9)",
                                is.na(ISCO08_code1) ~ "Unobserved"))

flastoccu <- indi_ns %>%
  tabyl(ISCO08_agg, trgroup) %>%
  adorn_percentages(denominator = "col") %>%
  adorn_pct_formatting(affix_sign = F) %>%
  slice(2, 4, 1, 3, 5) # custom order

flastoccu_all <- indi_ns %>%
  tabyl(ISCO08_agg) %>%
  adorn_percentages(denominator = "col") %>%
  adorn_pct_formatting(affix_sign = F) %>%
  slice(2, 4, 1, 3, 5) # custom order

freqlo <- cbind(flastoccu, all = flastoccu_all[, 3])
colnames(freqlo)[5] <- "all"

# rest of variables
indi_ns$disabil <- ifelse(indi_ns$ide_data %in% disab$ide_data, "Yes", "No")

f_asfd <- indi_ns %>%
  group_by(trgroup) %>%
  summarise(
    age = mean(eta),
    sex_w = mean(genere == "F") * 100,
    foreign = mean(flg_italiano == "no") * 100,
    disability = mean(disabil == "Yes") * 100
  ) %>%
  pivot_longer(-1) %>%
  pivot_wider(names_from = "trgroup", values_from = "value") %>%
  adorn_rounding(1) %>%
  slice(1,4,3,2) # setting custom order

f_asfd$all <- c(
  mean(indi_ns$eta), mean(indi_ns$genere == "F") * 100,
  mean(indi_ns$flg_italiano == "no") * 100,
  mean(indi_ns$disabil == "Yes") * 100
)

f_asfdx <- adorn_rounding(f_asfd, 1)

colnames(f_asfdx)[1] <- "variable"
colnames(freqstudy)[1] <- "variable"
colnames(freqlo)[1] <- "variable"

samplesize <- indi_ns %>%
  tabyl(trgroup) %>%
  adorn_totals()
ssize <- as.data.frame(cbind("n", (t(samplesize$n))))
colnames(ssize) <- c("variable", "A", "B", "C", "all")
summaryf <- rbind(f_asfdx, freqstudy, freqlo, ssize)

#openxlsx::write.xlsx(summaryf, 'intermediate/script01/T2_summaryftable.xlsx')


## Subsamples according to treatment received --------------------
indi_ns_ss1 <- subset(indi_ns, indi_ns$scoringD1_0 < 0.2) # subsample for treatment 1
indi_ns_ss2 <- subset(indi_ns, indi_ns$scoringD1_0 == 0 | indi_ns$scoringD1_0 > 0) # subsample for treatment 2

indi_ns_ss1$ide_data <- paste0(indi_ns_ss1$ide_soggetto, "_", indi_ns_ss1$axl_data_attribuzione)
indi_ns_ss1$disabil <- ifelse(indi_ns_ss1$ide_data %in% disab$ide_data, "Yes", "No")
indi_ns_ss1$disability <- as.factor(indi_ns_ss1$disabil)

indi_ns_ss2$ide_data <- paste0(indi_ns_ss2$ide_soggetto, "_", indi_ns_ss2$axl_data_attribuzione)
indi_ns_ss2$disabil <- ifelse(indi_ns_ss2$ide_data %in% disab$ide_data, "Yes", "No")
indi_ns_ss2$disability <- as.factor(indi_ns_ss2$disabil)


# Saving -----------
saveRDS(indi_ns_ss1, "intermediate/script01/indi_ns_ss1_190225.RDS")
saveRDS(indi_ns_ss2, "intermediate/script01/indi_ns_ss2_190225.RDS")
saveRDS(longest, "intermediate/script01/longest_190225.RDS")

#saveRDS(indi_ns_ss1, "indi_ns_ss1_saved.RDS")
#saveRDS(indi_ns_ss2, "indi_ns_ss2_saved.RDS")
#saveRDS(longest, "longest_saved.RDS")

