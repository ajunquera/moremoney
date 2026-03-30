# ...............................................................................
# ASSEGNO PER IL LAVORO - 02-02 Pick best estimates according to AMSE estimates
# Author: Álvaro F. Junquera
# ...............................................................................

library(tidyverse)
library(readxl)

outs <- data.frame(outcome = c("Y1", "Y1", "Y2", "Y2", "Y3", "Y3", "Y4", "Y4",
                               "Y_HS", "Y_HS", "Y_HT", "Y_HT"))

# AMSE estimates ---------
## For D1 ----------------
amse_allouts_d1 <- read_excel("intermediate/script02/stata/resu_amse_1.xlsx")

amse_allouts_d1_l <- amse_allouts_d1 %>%
  pivot_longer(cols = c1:c12,
               names_to = "case",
               values_to = "estimate") %>%
  bind_cols(., outs) %>%
  mutate(degree = rep(1:2, 6)) %>%
  select(-case) %>%
  relocate(outcome) %>%
  relocate(degree, .after = outcome) %>%
  as.data.frame()

amse_allouts_d1_l_best <- amse_allouts_d1_l %>%
  group_by(outcome) %>%
  slice_min(estimate, n = 1)

amse_allouts_d1_l_best$degree # all the best estimates come from a linear degree regression

## For D2 ------------
amse_allouts_d2 <- read_excel("intermediate/script02/stata/resu_amse_2.xlsx")

amse_allouts_d2_l <- amse_allouts_d2 %>%
  pivot_longer(cols = c1:c12,
               names_to = "case",
               values_to = "estimate") %>%
  bind_cols(., outs) %>%
  mutate(degree = rep(1:2, 6)) %>%
  select(-case) %>%
  relocate(outcome) %>%
  relocate(degree, .after = outcome) %>%
  as.data.frame()

amse_allouts_d2_l_best <- amse_allouts_d2_l %>%
  group_by(outcome) %>%
  slice_min(estimate, n = 1)

amse_allouts_d2_l_best$degree # all the best estimates come from a linear degree regression


# Assistance reception - table for LATE estimates -------------
## For D1 -----------
est_jshours_d1 <- read_delim("intermediate/script02/E11_cont_D1_jshours_nocl.csv")

table_js_d1 <- data.frame(col1 = c("D_1", "", "Bandwidth", "n_h"),
                        Y_HS = est_jshours_d1$`(1)`[1:4])

## For D2 ---------
est_jshours_d2 <- read_delim("intermediate/script02/E13_cont_D2_jshours_nocl.csv")

table_js_d2 <- data.frame(col1 = c("D_2", "", "Bandwidth", "n_h"),
                          Y_HS = est_jshours_d2$`(1)`[1:4])

## Join both and save -----------
table_jsa <- bind_rows(table_js_d1, table_js_d2)

table_jsa_compact <- table_jsa %>%
  mutate(
    across(Y_HS, ~ case_when(
      row_number() == 1 ~ str_c(.x[row_number() %in% c(1,2)], collapse = "<br>"),
      row_number() == 5 ~ str_c(.x[row_number() %in% c(5,6)], collapse = "<br>"),
      TRUE ~ .x
    ))
  ) %>%
  slice(-c(2,6))

write.csv(table_jsa_compact, "intermediate/script02/table_jsa_compact.csv")


# Employment - table for LATE estimates -------------
## For D1 --------
## Read
est_y1_d1 <- read_delim("intermediate/script02/E3_cont_D1_post6_nocl.csv")
est_y2_d1 <- read_delim("intermediate/script02/E4_cont_D1_post12_nocl.csv")
est_y3_d1 <- read_delim("intermediate/script02/E5_cont_D1_post18_nocl.csv")
est_y4_d1 <- read_delim("intermediate/script02/E6_cont_D1_post24_nocl.csv")

## Pick estimates
table3_d1 <- data.frame(col1 = c("D_1", "", "Bandwidth", "n_h"),
                        Y_1 = est_y1_d1$`(1)`[1:4],
                        Y_2 = est_y2_d1$`(1)`[1:4],
                        Y_3 = est_y3_d1$`(1)`[1:4],
                        Y_4 = est_y4_d1$`(1)`[1:4])

## For D2 --------
## Read

est_y1_d2 <- read_delim("intermediate/script02/E7_cont_D2_post6_nocl.csv")
est_y2_d2 <- read_delim("intermediate/script02/E8_cont_D2_post712_nocl.csv")
est_y3_d2 <- read_delim("intermediate/script02/E9_cont_D2_post1318_nocl.csv")
est_y4_d2 <- read_delim("intermediate/script02/E10_cont_D2_post1924_nocl.csv")

## Pick estimates
table3_d2 <- data.frame(col1 = c("D_2", "", "Bandwidth", "n_h"),
                        Y_1 = est_y1_d2$`(1)`[1:4],
                        Y_2 = est_y2_d2$`(1)`[1:4],
                        Y_3 = est_y3_d2$`(1)`[1:4],
                        Y_4 = est_y4_d2$`(1)`[1:4])

## Join both and save ----------
table3 <- bind_rows(table3_d1, table3_d2)

table3_compact <- table3 %>%
  mutate(
    across(Y_1:Y_4, ~ case_when(
      row_number() == 1 ~ str_c(.x[row_number() %in% c(1,2)], collapse = "<br>"),
      row_number() == 5 ~ str_c(.x[row_number() %in% c(5,6)], collapse = "<br>"),
      TRUE ~ .x
    ))
  ) %>%
  slice(-c(2,6)) %>%
  select(-`...1`)

write.csv(table3, "intermediate/script02/table3.csv")
write.csv(table3_compact, "intermediate/script02/table3_compact.csv")


# Supplementary ----------
## Assistance reception - table for LATE estimates -------------
### For D1 -----------
est_training_d1 <- read_delim("intermediate/script02/E12_cont_D1_trhours_nocl.csv")

table_tr_d1 <- data.frame(col1 = c("D_1", "", "Bandwidth", "n_h"),
                          Y_HT = est_training_d1$`(1)`[1:4])

### For D2 ---------
est_training_d2 <- read_delim("intermediate/script02/E14_cont_D2_trhours_nocl.csv")

table_tr_d2 <- data.frame(col1 = c("D_2", "", "Bandwidth", "n_h"),
                          Y_HT = est_training_d2$`(1)`[1:4])

### Join both and save -----------
table_tr <- bind_rows(table_tr_d1, table_tr_d2)

table_tr_compact <- table_tr %>%
  mutate(
    across(Y_HT, ~ case_when(
      row_number() == 1 ~ str_c(.x[row_number() %in% c(1,2)], collapse = "<br>"),
      row_number() == 5 ~ str_c(.x[row_number() %in% c(5,6)], collapse = "<br>"),
      TRUE ~ .x
    ))
  ) %>%
  slice(-c(2,6))

write.csv(table_tr_compact, "intermediate/script02/table_tr_compact.csv")
