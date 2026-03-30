# ...............................................................................
# ASSEGNO PER IL LAVORO - 05 Heterogeneity analyses with forests
# Author: Álvaro F. Junquera
# ...............................................................................

library(rdrobust)
library(rdhte)

library(grf)
library(fastpolicytree)
library(policytree)

library(tidyverse)
library(readxl)
library(sandwich)

library(tidymodels)
library(rpart.plot)

checknas <- function(data_frame) {
  na_counts <- sapply(data_frame, function(x) sum(is.na(x)))
  return(na_counts)
}

source("source/utils.R")

# 0. Reading -------------
indi_ns_ss1_c <- readRDS("intermediate/script02/indi_ns_ss1_c.RDS")
indi_ns_ss2_c <- readRDS("intermediate/script02/indi_ns_ss2_c.RDS")

provcomuni <- read_excel("inputdata/provcomuniVEN.xlsx") %>%
  drop_na(municipality) %>%
  mutate(munici = case_when(municipality == "San Donà di Piave" ~ "San Dona' di Piave",
                            municipality == "Schio" ~ "Schio-Thiene",
                            TRUE ~ municipality)) %>%
  select(-municipality)

# 1. Preparing variables -------------
## Outcomes --------
i1_Y1 <- indi_ns_ss1_c$post_interval6
i1_Y2 <- indi_ns_ss1_c$post_interval712
i1_Y3 <- indi_ns_ss1_c$post_interval1318
i1_Y4 <- indi_ns_ss1_c$post_interval1924
i1_Ycumu <- indi_ns_ss1_c$post_interval6 + indi_ns_ss1_c$post_interval712 +
  indi_ns_ss1_c$post_interval1318 + indi_ns_ss1_c$post_interval1924

i2_Y1 <- indi_ns_ss2_c$post_interval6
i2_Y2 <- indi_ns_ss2_c$post_interval712
i2_Y3 <- indi_ns_ss2_c$post_interval1318
i2_Y4 <- indi_ns_ss2_c$post_interval1924
i2_Ycumu <- indi_ns_ss2_c$post_interval6 + indi_ns_ss2_c$post_interval712 +
  indi_ns_ss2_c$post_interval1318 + indi_ns_ss2_c$post_interval1924

## Covariates --------
# Possible covariates:
# - Participant covariates:
#   - Categorical: province of residence, sex, level of education, disability,
#     last occupation, last sector, last type of contract
#   - Continuous: age, duration of the last contract (days)
# - Provider covariate: provider

### D1
indi_ns_ss1_c <- indi_ns_ss1_c %>%
  left_join(., provcomuni, by = c("cpi_titolare" = "munici")) %>%
  mutate(province = if_else(is.na(province), "NA", province))

indi_ns_ss1_c$leveledu <- if_else(is.na(indi_ns_ss1_c$studio2_grouped2), "Unobserved", indi_ns_ss1_c$studio2_grouped2)

provs_atleast300 <- indi_ns_ss1_c %>%
  janitor::tabyl(axl_ente) %>%
  filter(n > 299) %>%
  arrange(desc(n)) %>%
  pull(axl_ente)

indi_ns_ss1_c <- indi_ns_ss1_c %>%
  mutate(provider = if_else(axl_ente %in% provs_atleast300, axl_ente, "SmallProvider"))

indi_ns_ss1_c$lastocc <- str_sub(indi_ns_ss1_c$edu_last_contract, 1, 1) # 1 digit

indi_ns_ss1_c <- indi_ns_ss1_c %>%
  mutate(lastocc = str_sub(edu_last_contract, 1, 1)) %>%
  mutate(lastoccu = case_when(lastocc == "9" | is.na(lastocc) ~ "NA",
                              TRUE ~ lastocc)) %>%
  select(-lastocc)

# indi_ns_ss1_c$occu <- paste0(str_sub(indi_ns_ss1$edu_last_contract, 1, 1), str_sub(indi_ns_ss1$edu_last_contract, 3, 3)) # 2 digits

indi_ns_ss1_c <- indi_ns_ss1_c %>%
  mutate(lasttype = case_when(pre_oss_ultimo_contratto2_ass == "a- Cti" ~ "CTI",
                              pre_oss_ultimo_contratto2_ass == "b- Cap" ~ "AP",
                              pre_oss_ultimo_contratto2_ass == "c- Ctd" ~ "TD",
                              pre_oss_ultimo_contratto2_ass == "d- Som" ~ "SOM",
                              pre_no_oss_ultimo_contratto2_ass == "e- Int" ~ "INT",
                              pre_no_oss_ultimo_contratto2_ass == "f- Dom" ~ "DOM",
                              pre_no_oss_ultimo_contratto2_ass == "g- Par" ~ "PAR",
                              TRUE ~ "Unobserved"))

i1_Xprev <- indi_ns_ss1_c %>%
  ungroup() %>%
  select(province, genere, eta, leveledu, disabil,
         lastoccu, sectorVL, lasttype, pre_oss_ultimo_durata, pre_no_oss_ultimo_durata) # provider removed!

## Columns in matrix form for {grf}
i1_X <- model.matrix(~ 0 + ., i1_Xprev)

### D2
indi_ns_ss2_c <- indi_ns_ss2_c %>%
  left_join(., provcomuni, by = c("cpi_titolare" = "munici")) %>%
  mutate(province = if_else(is.na(province), "NA", province))

indi_ns_ss2_c$leveledu <- if_else(is.na(indi_ns_ss2_c$studio2_grouped2), "Unobserved", indi_ns_ss2_c$studio2_grouped2)

provs_atleast300_s2 <- indi_ns_ss2_c %>%
  janitor::tabyl(axl_ente) %>%
  filter(n > 299) %>%
  arrange(desc(n)) %>%
  pull(axl_ente)

indi_ns_ss2_c <- indi_ns_ss2_c %>%
  mutate(provider = if_else(axl_ente %in% provs_atleast300_s2, axl_ente, "SmallProvider"))

indi_ns_ss2_c$lastocc <- str_sub(indi_ns_ss2_c$edu_last_contract, 1, 1) # 1 digit

indi_ns_ss2_c <- indi_ns_ss2_c %>%
  mutate(lastocc = str_sub(edu_last_contract, 1, 1)) %>%
  mutate(lastoccu = case_when(lastocc == "9" | is.na(lastocc) ~ "NA",
                              TRUE ~ lastocc)) %>%
  select(-lastocc)

indi_ns_ss2_c <- indi_ns_ss2_c %>%
  mutate(lasttype = case_when(pre_oss_ultimo_contratto2_ass == "a- Cti" ~ "CTI",
                              pre_oss_ultimo_contratto2_ass == "b- Cap" ~ "AP",
                              pre_oss_ultimo_contratto2_ass == "c- Ctd" ~ "TD",
                              pre_oss_ultimo_contratto2_ass == "d- Som" ~ "SOM",
                              pre_no_oss_ultimo_contratto2_ass == "e- Int" ~ "INT",
                              pre_no_oss_ultimo_contratto2_ass == "f- Dom" ~ "DOM",
                              pre_no_oss_ultimo_contratto2_ass == "g- Par" ~ "PAR",
                              TRUE ~ "Unobserved"))

i2_Xprev <- indi_ns_ss2_c %>%
  ungroup() %>%
  select(province, genere, eta, leveledu, disabil,
         lastoccu, sectorVL, lasttype, pre_oss_ultimo_durata, pre_no_oss_ultimo_durata) # provider removed!

## Columns in matrix form for {grf}
i2_X <- model.matrix(~ 0 + ., i2_Xprev)

## Treatment variables --------
cutoff <- 0

i1_Z <- indi_ns_ss1_c$scoringD1_0
i1_W <- as.numeric(i1_Z >= cutoff)

i2_Z <- indi_ns_ss2_c$scoringD2_0
i2_W <- as.numeric(i2_Z >= cutoff)


# 2. Running the forest --------
### (Kernel) linear model forest

### Forests for D1 and D2 -----------
indi_ns_ss1_c$post_Ycumu <- indi_ns_ss1_c$post_interval6 + indi_ns_ss1_c$post_interval712 +
  indi_ns_ss1_c$post_interval1318 + indi_ns_ss1_c$post_interval1924

indi_ns_ss2_c$post_Ycumu <- indi_ns_ss2_c$post_interval6 + indi_ns_ss2_c$post_interval712 +
  indi_ns_ss2_c$post_interval1318 + indi_ns_ss2_c$post_interval1924

d1yC_3ft <- forest_tree(outcome = "post_Ycumu",
                        running = "scoringD1_0",
                        cutoff = 0,
                        df = indi_ns_ss1_c,
                        vector_y = i1_Ycumu,
                        matrix_treatrunning = cbind(i1_W, i1_Z),
                        matrix_x = i1_X,
                        treedepth = 3)

View(d1yC_3ft[[2]])

indi_ns_ss1_c_BIND <- cbind(post_Ycumu = indi_ns_ss1_c$post_Ycumu,
                            scoringD1_0 = indi_ns_ss1_c$scoringD1_0,
                            i1_W, i1_Z,
                            i1_X) |>
  as.data.frame()

indi_ns_ss1_c_BIND <- indi_ns_ss1_c_BIND |>
  mutate(sub1 = pre_oss_ultimo_durata <= 110,
         sub2 = pre_oss_ultimo_durata <= 110 & leveledu344_uppersecG == 0,
         sub3 = pre_oss_ultimo_durata > 110 & pre_oss_ultimo_durata <= 670,
         sub4 = pre_oss_ultimo_durata <= 110 & leveledu344_uppersecG <= 0 & leveledu353_uppersecP_3y == 0,
         sub5 = pre_oss_ultimo_durata <= 110 & leveledu344_uppersecG > 0 & lasttypeTD == 0,
         sub6 = pre_oss_ultimo_durata > 110 & pre_oss_ultimo_durata <= 670 & sectorVLTerziarioAvanzato == 0,
         sub7 = pre_oss_ultimo_durata > 110 & pre_oss_ultimo_durata > 670 & lastoccu4 == 0) |>
  mutate(across(.cols = sub1:sub7, .fns = as.factor))

subshte <- c(indi_ns_ss1_c_BIND$sub1,
             indi_ns_ss1_c_BIND$sub2,
             indi_ns_ss1_c_BIND$sub3,
             indi_ns_ss1_c_BIND$sub4,
             indi_ns_ss1_c_BIND$sub5,
             indi_ns_ss1_c_BIND$sub6,
             indi_ns_ss1_c_BIND$sub7)

rdmodel_F1 <- rdhte(y = indi_ns_ss1_c_BIND$post_Ycumu,
                    x = indi_ns_ss1_c_BIND$scoringD1_0,
                    covs.hte = subshte)

summary(rdmodel_F1)


d2yC_3ft <- forest_tree(outcome = "post_Ycumu",
                        running = "scoringD2_0",
                        cutoff = 0,
                        df = indi_ns_ss2_c,
                        vector_y = i2_Ycumu,
                        matrix_treatrunning = cbind(i2_W, i2_Z),
                        matrix_x = i2_X,
                        treedepth = 3)

View(d2yC_3ft[[2]])

#### Multiple hypothesis testing corrections
p.adjust(p = d2yC_3ft[[2]]$diff_pvalue,
         method = "holm")


#### Policy value
policyvalues <- data.frame(contrast = c(rep("D1", 4), rep("D2", 4)),
                           outcome = c(rep(c("Y1", "Y2", "Y3", "Y4"), 2)),
                           value = c(d1y1_3[[1]], d1y2_3[[1]], d1y3_3[[1]], d1y4_3[[1]],
                                     d2y1_3[[1]], d2y2_3[[1]], d2y3_3[[1]], d2y4_3[[1]]))

policyvalues <- policyvalues |>
  mutate(value_r = round(value, 2))

rf_sg_d1 <- rbind(cbind(outcome = rep("Y1", 4), d1y1_3[[2]]),
                  cbind(outcome = rep("Y2", 4), d1y2_3[[2]]),
                  cbind(outcome = rep("Y3", 4), d1y3_3[[2]]),
                  cbind(outcome = rep("Y4", 4), d1y4_3[[2]]))

rf_sg_d1 <- rf_sg_d1 |>
  mutate(var1 = case_when(var1 == "pre_oss_ultimo_durata" ~ "duration_last_good_emp",
                          var1 == "pre_no_oss_ultimo_durata" ~ "duration_last_bad_emp",
                          TRUE ~ var1),
         var2 = case_when(var2 == "pre_oss_ultimo_durata" ~ "duration_last_good_emp",
                          var2 == "pre_no_oss_ultimo_durata" ~ "duration_last_bad_emp",
                          TRUE ~ var2),
         var3 = case_when(var3 == "pre_oss_ultimo_durata" ~ "duration_last_good_emp",
                          var3 == "pre_no_oss_ultimo_durata" ~ "duration_last_bad_emp",
                          TRUE ~ var3))

rf_sg_d2 <- rbind(cbind(outcome = rep("Y1", 4), d2y1_3[[2]]),
                  cbind(outcome = rep("Y2", 4), d2y2_3[[2]]),
                  cbind(outcome = rep("Y3", 4), d2y3_3[[2]]),
                  cbind(outcome = rep("Y4", 4), d2y4_3[[2]]))

rf_sg_d2 <- rf_sg_d2 |>
  mutate(var1 = case_when(var1 == "pre_oss_ultimo_durata" ~ "duration_last_good_emp",
                          var1 == "pre_no_oss_ultimo_durata" ~ "duration_last_bad_emp",
                          TRUE ~ var1),
         var2 = case_when(var2 == "pre_oss_ultimo_durata" ~ "duration_last_good_emp",
                          var2 == "pre_no_oss_ultimo_durata" ~ "duration_last_bad_emp",
                          TRUE ~ var2),
         var3 = case_when(var3 == "pre_oss_ultimo_durata" ~ "duration_last_good_emp",
                          var3 == "pre_no_oss_ultimo_durata" ~ "duration_last_bad_emp",
                          TRUE ~ var3))

write.csv(policyvalues, "intermediate/script05/policyvalues.csv")
write.csv(rf_sg_d1, "intermediate/script05/rf_sg_d1.csv")
write.csv(rf_sg_d2, "intermediate/script05/rf_sg_d2.csv")

### Multiple hypotheses corrections -----------
rf_sg_d1 <- read_csv("intermediate/script05/rf_sg_d1.csv")
rf_sg_d2 <- read_csv("intermediate/script05/rf_sg_d2.csv")


rf_sg_d1$pholm <- p.adjust(rf_sg_d1$diff_pvalue, method = "holm")
rf_sg_d1$phochberg <- p.adjust(rf_sg_d1$diff_pvalue, method = "hochberg")
rf_sg_d1$phommel <- p.adjust(rf_sg_d1$diff_pvalue, method = "hommel")

sum(rf_sg_d1$pholm < 0.1)
sum(rf_sg_d1$phochberg < 0.1)
sum(rf_sg_d1$phommel < 0.1)



rf_sg_d2$pholm <- p.adjust(rf_sg_d2$diff_pvalue, method = "holm")
rf_sg_d2$phochberg <- p.adjust(rf_sg_d2$diff_pvalue, method = "hochberg")
rf_sg_d2$phommel <- p.adjust(rf_sg_d2$diff_pvalue, method = "hommel")

sum(rf_sg_d2$pholm < 0.1, na.rm = T)
sum(rf_sg_d2$phochberg < 0.1, na.rm = T)
sum(rf_sg_d2$phommel < 0.1, na.rm = T)



### Other interpretability measures of the forest ---------

#### ATE by quintile ----------
bandwidth <- rdrobust::rdbwselect(y = i1_Y4, x = i1_Z,
                                  c = cutoff)$bws[1,1]

# Compute kernel weights for a triangular kernel.
dist <- abs((indi_ns_ss1_c$scoringD1_0 - cutoff) / bandwidth)
sample.weights <- (1 - dist) * (dist <= 1) / bandwidth

# Estimate a local linear regression with the running variable Z conditional on covariates X = x:
# Y = c(x) + tau(x) W + b(x) Z.
# Specify gradient.weights = c(1, 0) to target heterogeneity in the RDD coefficient tau(x).
# Also, fit forest on subset with non-zero weights for faster estimation.

subset <- sample.weights > 0

num.rankings <- 5
num.folds <- 10
folds <- sort(seq(nrow(i1_X[subset, ])) %% num.folds) + 1


lmf_Y4 <- lm_forest(X = i1_X[subset, ],
                    Y = i1_Y4[subset],
                    W = cbind(i1_W, i1_Z)[subset, ], # W (K=1) is treatment indicator, Z (K=2) is running variable
                    sample.weights = sample.weights[subset],
                    gradient.weights = c(1, 0),
                    clusters = folds)

# Retrieve out-of-bag predictions.
tau.hat <- predict(lmf_Y4)$predictions[, 1, ] # coefficient of W (K=1)
y.hat.lmf <- lmf_Y4$Y.hat

# Rank observations *within each fold* into quintiles according to their CATE predictions.
ranking <- rep(NA, nrow(i1_X[subset, ]))
for (fold in seq(num.folds)) {
  tau.hat.quantiles <- quantile(tau.hat[folds == fold], probs = seq(0, 1, by=1/num.rankings))
  ranking[folds == fold] <- cut(tau.hat[folds == fold], tau.hat.quantiles, include.lowest=TRUE,labels=seq(num.rankings))
}


# Formula y ~ 0 + ranking + ranking:w
dfrdd1 <- cbind(tau.hat, cbind(i1_W, i1_Z)[subset, ], i1_Xprev[subset, ], rank = factor(ranking))

ols.ate <- lm(tau.hat ~ 0 + rank,
              data = dfrdd1,
              weights = sample.weights[subset]) # W is treatment indicator, Z is score
summary(ols.ate)

ols.ate.test <- lmtest::coeftest(ols.ate, vcov = sandwich::vcovHC(ols.ate, type='HC2'))
ols.ate.test

ols.ate.ok <- data.frame("ols", paste0("Q", seq(num.rankings)), ols.ate.test[1:5, 1:2])
ols.ate.ok

rownames(ols.ate.ok) <- NULL # just for display
colnames(ols.ate.ok) <- c("method", "ranking", "estimate", "std.err")
ols.ate.ok

ggplot(ols.ate.ok) +
  aes(x = ranking, y = estimate) +
  geom_point(position=position_dodge(0.2)) +
  geom_errorbar(aes(ymin=estimate-2*std.err, ymax=estimate+2*std.err), width=.2, position=position_dodge(0.2)) +
  ylab("") + xlab("") +
  theme_minimal() +
  theme(legend.position="bottom", legend.title = element_blank())

### Testing considering multiple hypothesis testing problem
# Auxiliary function to computes adjusted p-values
# following the Romano-Wolf method.
# For a reference, see http://ftp.iza.org/dp12845.pdf page 8
#  t.orig: vector of t-statistics from original model
#  t.boot: matrix of t-statistics from bootstrapped models
romano_wolf_correction <- function(t.orig, t.boot) {
  abs.t.orig <- abs(t.orig)
  abs.t.boot <- abs(t.boot)
  abs.t.sorted <- sort(abs.t.orig, decreasing = TRUE)

  max.order <- order(abs.t.orig, decreasing = TRUE)
  rev.order <- order(max.order)

  M <- nrow(t.boot)
  S <- ncol(t.boot)

  p.adj <- rep(0, S)
  p.adj[1] <- mean(apply(abs.t.boot, 1, max) > abs.t.sorted[1])
  for (s in seq(2, S)) {
    cur.index <- max.order[s:S]
    p.init <- mean(apply(abs.t.boot[, cur.index, drop=FALSE], 1, max) > abs.t.sorted[s])
    p.adj[s] <- max(p.init, p.adj[s-1])
  }
  p.adj[rev.order]
}

summary_rw_lm <- function(model, indices=NULL, cov.type="HC2", num.boot=10000) {

  if (is.null(indices)) {
    indices <- 1:nrow(coef(summary(model)))
  }
  # Grab the original t values.
  summary <- coef(summary(model))[indices,,drop=FALSE]
  t.orig <- summary[, "t value"]

  # Null resampling.
  # This is a trick to speed up bootstrapping linear models.
  # Here, we don't really need to re-fit linear regressions, which would be a bit slow.
  # We know that betahat ~ N(beta, Sigma), and we have an estimate Sigmahat.
  # So we can approximate "null t-values" by
  #  - Draw beta.boot ~ N(0, Sigma-hat) --- note the 0 here, this is what makes it a *null* t-value.
  #  - Compute t.boot = beta.boot / sqrt(diag(Sigma.hat))
  Sigma.hat <- sandwich::vcovHC(model, type=cov.type)[indices, indices]
  se.orig <- sqrt(diag(Sigma.hat))
  num.coef <- length(se.orig)
  beta.boot <- MASS::mvrnorm(n=num.boot, mu=rep(0, num.coef), Sigma=Sigma.hat)
  t.boot <- sweep(beta.boot, 2, se.orig, "/")
  p.adj <- romano_wolf_correction(t.orig, t.boot)

  result <- cbind(summary[,c(1,2,4),drop=F], p.adj)
  colnames(result) <- c('Estimate', 'Std. Error', 'Orig. p-value', 'Adj. p-value')
  result
}

res <- summary_rw_lm(ols.ate, indices=2:num.rankings)
rownames(res) <- paste("Rank", 2:num.rankings, "- Rank 1") # just for display
res

##### Average values of covariates by quantile --------------
#tau.hat ~ 0 + rank,
data <- cbind(i1_X[subset, ], rank = factor(ranking))

covariates <- colnames(data)
covariates <- gsub("\\+", "more", covariates)
colnames(data) <- covariates

covariates <- covariates[covariates != "rank"]


df <- mapply(function(covariate) {
  # Looping over covariate names
  # Compute average covariate value per ranking (with correct standard errors)
  fmla <- formula(paste0(covariate, "~ 0 + ranking"))
  ols <- lm(fmla, data=transform(data, ranking=factor(ranking)))
  ols.res <- lmtest::coeftest(ols, vcov=vcovHC(ols, "HC2"))

  # Retrieve results
  avg <- ols.res[,1]
  stderr <- ols.res[,2]

  # Tally up results
  data.frame(covariate, avg, stderr, ranking=paste0("Q", seq(num.rankings)),
             # Used for coloring
             scaling=pnorm((avg - mean(avg))/sd(avg)),

             # We will order based on how much variation is 'explain' by the averages
             # relative to the total variation of the covariate in the data
             variation=sd(avg) / sd(data[,covariate]),
             # String to print in each cell in heatmap below
             labels=paste0(signif(avg, 3), "\n", "(", signif(stderr, 3), ")"))
}, covariates, SIMPLIFY = FALSE)

df <- do.call(rbind, df)

# a small optional trick to ensure heatmap will be in decreasing order of 'variation'
#df$covariate <- reorder(df$covariate, order(df$variation))

# plot heatmap
ggplot(df |> filter(str_starts(covariate, c("^province|genere")))) +
  aes(ranking, covariate) +
  geom_tile(aes(fill = scaling)) +
  geom_text(aes(label = labels)) +
  scale_fill_gradient(low = "#E1BE6A", high = "#40B0A6") +
  theme_minimal() +
  ylab("") + xlab("CATE estimate ranking") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text=element_text(size=11))

ggplot(df |> filter(str_starts(covariate, c("^leveledu|disabil")))) +
  aes(ranking, covariate) +
  geom_tile(aes(fill = scaling)) +
  geom_text(aes(label = labels)) +
  scale_fill_gradient(low = "#E1BE6A", high = "#40B0A6") +
  theme_minimal() +
  ylab("") + xlab("CATE estimate ranking") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text=element_text(size=11))

ggplot(df |> filter(str_starts(covariate, c("^lastoccu|sectorVL")))) +
  aes(ranking, covariate) +
  geom_tile(aes(fill = scaling)) +
  geom_text(aes(label = labels)) +
  scale_fill_gradient(low = "#E1BE6A", high = "#40B0A6") +
  theme_minimal() +
  ylab("") + xlab("CATE estimate ranking") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text=element_text(size=11))

ggplot(df |> filter(str_starts(covariate, c("^lasttype|eta|pre_oss_ultimo_durata|pre_no_oss_ultimo_durata")))) +
  aes(ranking, covariate) +
  geom_tile(aes(fill = scaling)) +
  geom_text(aes(label = labels)) +
  scale_fill_gradient(low = "#E1BE6A", high = "#40B0A6") +
  theme_minimal() +
  ylab("") + xlab("CATE estimate ranking") +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text=element_text(size=11))

#### Variable Importance --------
var_imp <- variable_importance(lmf_Y4)
names(var_imp) <- names(as.data.frame(i1_X))
sorted_var_imp <- sort(var_imp, decreasing = TRUE)
sorted_var_imp[1:5]  # showing only first few


#### Best Linear Projection -----------
### Manual interpretation: Only main effects
blp_lmfY4 <- lm(tau.hat ~ i1_X[subset, ],
                weights = sample.weights[subset])

blp_lmfY4_co <- summary(blp_lmfY4)

blp_lmfY4_results <- data.frame(coeffi = blp_lmfY4_co$coefficients[, 2],
                                pvalue = blp_lmfY4_co$coefficients[, 4]) |>
  filter(pvalue < 0.05) |>
  arrange(desc(coeffi))


##### Manual interpretation: Up to 2-terms interactions ------
#### unpenalized LM
i1_X2 <- model.matrix(~ 0 + (.)^2, i1_Xprev)

blp_lmfY4_pairw <- lm(tau.hat ~ i1_X2[subset, ],
                weights = sample.weights[subset])

blp_lmfY4_pairw_co <- summary(blp_lmfY4_pairw)

blp_lmfY4_pw_results <- data.frame(coeffi = blp_lmfY4_pairw_co$coefficients[, 2],
                                pvalue = blp_lmfY4_pairw_co$coefficients[, 4]) |>
  filter(pvalue < 0.05) |>
  arrange(desc(coeffi)) # this under-consider the coefficients of continuous covariates!!

MSE_pairw <- mean((fitted(blp_lmfY4_pairw) - residuals(blp_lmfY4_pairw))^2)
MSE_pairw

#### Validation
indi_ns_ss1_2w_sub1 <- indi_ns_ss1_c |>
  filter(lastoccu == "8")

val_2w_sub1 <- rdrobust(
  y = indi_ns_ss1_2w_sub1$post_interval1924,
  x = indi_ns_ss1_2w_sub1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL, all = T)

summary(val_2w_sub1)

#### penalized LM
i1_X2_sp <- sparse.model.matrix(~ 0 + (.)^2, i1_Xprev)

blp_lmfY4_2w_sp <- glmnet(y = tau.hat, x = i1_X2_sp[subset, ],
                       family = "gaussian",
                       weights = sample.weights[subset])

cvfit_2w_sp <- cv.glmnet(y = tau.hat, x = i1_X2_sp[subset, ],
                   family = "gaussian",
                   weights = sample.weights[subset])

MSE_train_glmnet_2w <- mean((tau.hat - predict(cvfit_2w_sp, newx = i1_X2_sp[subset, ]))^2)
MSE_train_glmnet_2w

coefs_cvfit_2w_sp <- as.matrix(coef(cvfit_2w_sp, s = "lambda.min")) |>
  as.data.frame() |>
  rename(estimate = lambda.min) |>
  arrange(desc(estimate))

#### Validation
indi_ns_ss1_sub1_2w <- indi_ns_ss1_c |>
  filter(leveledu == "2_lowersec" &
           lastoccu == "2")

val_sub1_2w <- rdrobust(
  y = indi_ns_ss1_sub1_2w$post_interval1924,
  x = indi_ns_ss1_sub1_2w$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL, all = T)

summary(val_sub1_2w)

rdplot(y = indi_ns_ss1_sub1_2w$post_interval1924,
       x = indi_ns_ss1_sub1_2w$scoringD1_0)

##### Manual interpretation: Up to 3-terms interactions -------------
i1_X3 <- sparse.model.matrix(~ 0 + (.)^3, i1_Xprev)

blp_lmfY4_3w <- glmnet(y = tau.hat, x = i1_X3[subset, ],
                       family = "gaussian",
                       weights = sample.weights[subset])

cvfit <- cv.glmnet(y = tau.hat, x = i1_X3[subset, ],
                   family = "gaussian",
                   weights = sample.weights[subset])

MSE_train_glmnet <- mean((tau.hat - predict(cvfit, newx = i1_X3[subset, ]))^2)
MSE_train_glmnet

coefs_cvfit <- as.matrix(coef(cvfit, s = "lambda.min")) |>
  as.data.frame() |>
  rename(estimate = lambda.min) |>
  arrange(desc(estimate))

##### Validation ------------
indi_ns_ss1_sub1 <- indi_ns_ss1_c |>
  filter(province == "Vicenza" &
           leveledu == "353_uppersecP_3y" &
           lastoccu == "NA")

val_sub1 <- rdrobust(
  y = indi_ns_ss1_sub1$post_interval1924,
  x = indi_ns_ss1_sub1$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL, all = T)

summary(val_sub1)


indi_ns_ss1_sub2 <- indi_ns_ss1_c |>
  filter(province == "Padova" &
           leveledu == "344_uppersecG" &
           lasttype == "Unobserved")

val_sub2 <- rdrobust(
  y = indi_ns_ss1_sub2$post_interval1924,
  x = indi_ns_ss1_sub2$scoringD1_0,
  kernel = "triangular",
  c = 0, p = 1, bwselect = "mserd", cluster = NULL, all = T)

summary(val_sub2)

rdplot(y = indi_ns_ss1_sub2$post_interval1924,
       x = indi_ns_ss1_sub2$scoringD1_0,
       nbins = c(43, 57))


#### Surrogate tree -----------
df_used <- cbind(tau.hat, i1_Xprev[subset, ])

### 0. Data splitting
set.seed(1988)
rfpreds_split <- initial_split(as.data.frame(df_used), prop = 7/10)

rfpreds_train <- training(rfpreds_split)
rfpreds_test <- testing(rfpreds_split)

### 1. Formula & preprocessing {recipe}
stree_rec <-
  recipe(tau.hat ~ ., data = rfpreds_train) %>%
  step_zv(all_predictors()) # zv = zero variance. Removes columns that have only 1 value (rare categories)

### 2. Specify model {parsnip}
stree_mod <-
  decision_tree(cost_complexity = 0.005) %>%
  set_mode("regression") %>%
  set_engine("rpart")

### 3. Join {parsnip} and {recipe}
stree_wflow <-
  workflow() %>%
  add_model(stree_mod) %>%
  add_recipe(stree_rec)

### 4. Estimation of coefficients
set.seed(2024)

stree_fit <-
  stree_wflow %>%
  fit(data = rfpreds_train)

### 5. Performance
### ...in training subset
stree_preds_train <-
  augment(stree_fit, rfpreds_train)

stree_preds_train %>%
  rsq(truth = tau.hat, estimate = .pred)

### ...in test subset
stree_preds_test <-
  augment(stree_fit, rfpreds_test)

stree_preds_test %>%
  rsq(truth = tau.hat, estimate = .pred)

### 6. Plotting
tree_fit_rpart <- extract_fit_engine(stree_fit)

rpart.plot(tree_fit_rpart)