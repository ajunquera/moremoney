# ...............................................................................
# ASSEGNO PER IL LAVORO - 04 Plots, heterogeneity, cream skimming
# Author: Álvaro F. Junquera (UAB)
# ...............................................................................

library(tidyverse)
library(rdrobust)
library(rddensity)
library(modelsummary)
library(openxlsx)

# Reading data -------------
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


# Some necessary functions -----------
starizer <- function(x) {
  if (x$pv[3] < 0.01) {
    x$stars <- "***"
  } else if (x$pv[3] < 0.05) {
    x$stars <- "**"
  } else if (x$pv[3] < 0.1) {
    x$stars <- "*"
  } else {
    x$stars <- ""
  }

  return(x)
}


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


# 1. Plots  --------------------
## Treatment 1 (D1) ----------------
### Outcome 0a: number of hours of job search measures received ---------
summary(indi_ns_ss1$jshours)

plot_jsh_qs_mv <- rdplot(
  y = indi_ns_ss1$jshours, x = indi_ns_ss1$scoringD1_0, binselect = "qsmv",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 13), shade = F
) # quantile-spaced
summary(plot_jsh_qs_mv)


### Generate rdplot using ggplot2: preparing inputs
# plot1 <- plot24qs_imse
plot1 <- plot_jsh_qs_mv


# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss1$scoringD1_0
y <- indi_ns_ss1$jshours
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y["HS"], ": hours of JS measures received"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[1], " (", S[1], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y["HS"]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 13)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 0.75) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 0.75) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/late_d1_yhs.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

### Outcome 0b: number of hours of training measures received ---------
summary(indi_ns_ss1$attiv_form_ore_prev)

plot_trh_qs_mv <- rdplot(
  y = indi_ns_ss1$attiv_form_ore_prev, x = indi_ns_ss1$scoringD1_0, binselect = "qsmv",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 136), shade = F
) # quantile-spaced
summary(plot_trh_qs_mv)


### Generate rdplot using ggplot2: preparing inputs
# plot1 <- plot24qs_imse
plot1 <- plot_trh_qs_mv


# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss1$scoringD1_0
y <- indi_ns_ss1$attiv_form_ore_prev
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y["HT"], ": hours of training received"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[1], " (", S[1], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y["HT"]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 136)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 0.75) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 0.75) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # Saved at 395x300


ggsave(filename = "intermediate/script04/plots/late_d1_yht.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


### Outcome 1: Number of worked days in post [1, 6] months ------
summary(indi_ns_ss1$post_interval6)
summary(indi_ns_ss1$scoringD1_0)

summary(indi_ns_ss1$scoringD1_0)
summary(indi_ns_ss1$post_interval6)

plot6es_mv <- rdplot(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0, binselect = "esmv",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 185), shade = F
) # evenly-spaced
summary(plot6es_mv)

plot6qs_mv <- rdplot(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0, binselect = "qsmv",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 185), shade = F
) # quantile-spaced
summary(plot6qs_mv)

plot6qs_imse <- rdplot(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0, binselect = "qs",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 185), ci = 95, hide = TRUE
) # quantile-spaced
summary(plot6qs_imse)

### Option A. rdplot
summary(indi_ns_ss1$post_interval6)

plot_s1d1_qs_mv <- rdplot(
  y = indi_ns_ss1$post_interval6, x = indi_ns_ss1$scoringD1_0, binselect = "qsmv",
  x.lim = c(-0.37, 0.2), y.lim = c(0, 185), shade = F,
  title = "Y1: worked days during months [1, 6]", x.label = "Score for D1 (S1)", y.label = "LSM of Y1",
  col.dots = "dimgrey"
) # quantile-spaced, saved at  500 x 379
summary(plot_s1d1_qs_mv)



### Option B. Generate rdplot using ggplot2: preparing inputs

# plot1 <- plot6qs_imse
plot1 <- plot6qs_mv

# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss1$scoringD1_0
y <- indi_ns_ss1$post_interval6
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y[1], ": worked days during months [+1, +6]"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[1], " (", S[1], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y[1]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 185)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 0.75) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 0.75) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # 395 x 300

ggsave(filename = "intermediate/script04/plots/late_d1_y1.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


### Outcome 2: Number of worked days in post [7, 12] months ------
summary(indi_ns_ss1$post_interval712)
summary(indi_ns_ss1$scoringD1_0)

summary(indi_ns_ss1$scoringD1_0)
summary(indi_ns_ss1$post_interval12)

plot12es_mv <- rdplot(
  y = indi_ns_ss1$post_interval712, x = indi_ns_ss1$scoringD1_0, binselect = "esmv",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 185), shade = F
) # evenly-spaced
summary(plot12es_mv)

plot12qs_mv <- rdplot(
  y = indi_ns_ss1$post_interval712, x = indi_ns_ss1$scoringD1_0, binselect = "qsmv",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 185), shade = F
) # quantile-spaced
summary(plot12qs_mv)

plot12qs_imse <- rdplot(
  y = indi_ns_ss1$post_interval12, x = indi_ns_ss1$scoringD1_0, binselect = "qs",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 375), ci = 95, hide = TRUE
) # quantile-spaced
summary(plot12qs_imse)


### Generate rdplot using ggplot2: preparing inputs
# plot1 <- plot12qs_imse
plot1 <- plot12qs_mv

# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss1$scoringD1_0
y <- indi_ns_ss1$post_interval712
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y[2], ": worked days during months [+7, +12]"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[1], " (", S[1], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y[2]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 185)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 0.75) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 0.75) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/late_d1_y2.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


### Outcome 3: Number of worked days in post [13, 18] months ------
summary(indi_ns_ss1$post_interval1318)
summary(indi_ns_ss1$scoringD1_0)

summary(indi_ns_ss1$scoringD1_0)
summary(indi_ns_ss1$post_interval18)

plot18es_mv <- rdplot(
  y = indi_ns_ss1$post_interval1318, x = indi_ns_ss1$scoringD1_0, binselect = "esmv",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 185), shade = F
) # evenly-spaced
summary(plot18es_mv)

plot18qs_mv <- rdplot(
  y = indi_ns_ss1$post_interval1318, x = indi_ns_ss1$scoringD1_0, binselect = "qsmv",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 185), shade = F
) # quantile-spaced
summary(plot18qs_mv)



### Generate rdplot using ggplot2: preparing inputs
# plot1 <- plot18qs_imse
plot1 <- plot18qs_mv

# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss1$scoringD1_0
y <- indi_ns_ss1$post_interval1318
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y[3], ": worked days during months [+13, +18]"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[1], " (", S[1], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y[3]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 185)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 0.75) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 0.75) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/late_d1_y3.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


### Outcome 4: Number of worked days in post [19, 24] months ------
summary(indi_ns_ss1$post_interval1924)
summary(indi_ns_ss1$scoringD1_0)

summary(indi_ns_ss1$scoringD1_0)
summary(indi_ns_ss1$post_interval24)

plot24es_mv <- rdplot(
  y = indi_ns_ss1$post_interval1924, x = indi_ns_ss1$scoringD1_0, binselect = "esmv",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 750), shade = F
) # evenly-spaced
summary(plot24es_mv)

plot24qs_mv <- rdplot(
  y = indi_ns_ss1$post_interval1924, x = indi_ns_ss1$scoringD1_0, binselect = "qsmv",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 185), shade = F
) # quantile-spaced
summary(plot24qs_mv)

plot24qs_imse <- rdplot(
  y = indi_ns_ss1$post_interval1924, x = indi_ns_ss1$scoringD1_0, binselect = "qs",
  x.lim = c(-0.40, 0.25), y.lim = c(0, 750), ci = 95, hide = TRUE
) # quantile-spaced
summary(plot24qs_imse)


### Generate rdplot using ggplot2: preparing inputs
# plot1 <- plot24qs_imse
plot1 <- plot24qs_mv


# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss1$scoringD1_0
y <- indi_ns_ss1$post_interval1924
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y[4], ": worked days during months [+19, +24]"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[1], " (", S[1], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y[4]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 185)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 0.75) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 0.75) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/late_d1_y4.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


## Treatment 2 (D2) ------------------------
### Outcome 0a: number of hours of job search measures received ---------
summary(indi_ns_ss2$jshours)

plot_jsh_qs_mv <- rdplot(
  y = indi_ns_ss2$jshours, x = indi_ns_ss2$scoringD2_0, binselect = "qsmv",
  x.lim = c(-0.25, 0.40), y.lim = c(0, 27), shade = F
) # quantile-spaced
summary(plot_jsh_qs_mv)


### Generate rdplot using ggplot2: preparing inputs
# plot1 <- plot24qs_imse
plot1 <- plot_jsh_qs_mv


# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss2$scoringD2_0
y <- indi_ns_ss2$jshours
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y[HS], ": hours of JS measures received"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[2], " (", S[2], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y["HS"]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 27)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 1) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/late_d2_yhs.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


### Outcome 0b: number of hours of training measures received ---------
summary(indi_ns_ss2$attiv_form_ore_prev)

plot_trh_qs_mv <- rdplot(
  y = indi_ns_ss2$attiv_form_ore_prev, x = indi_ns_ss2$scoringD2_0, binselect = "qsmv",
  x.lim = c(-0.25, 0.40), y.lim = c(0, 146), shade = F
) # quantile-spaced
summary(plot_trh_qs_mv)


### Generate rdplot using ggplot2: preparing inputs
# plot1 <- plot24qs_imse
plot1 <- plot_trh_qs_mv


# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss2$scoringD2_0
y <- indi_ns_ss2$attiv_form_ore_prev
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y[HT], ": hours of training received"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[2], " (", S[2], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y["HT"]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 146)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 1) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/late_d2_yht.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


### Outcome 1: post [1, 6] months ---------------------
summary(indi_ns_ss2$scoringD2_0)
summary(indi_ns_ss2$post_interval6)

plotD2_6_qs_imse <- rdplot(
  y = indi_ns_ss2$post_interval6, x = indi_ns_ss2$scoringD2_0, binselect = "qs",
  x.lim = c(-0.25, 0.40), y.lim = c(0, 200), ci = 95, hide = TRUE
) # quantile-spaced
summary(plotD2_6_qs_imse)

plotD2_6_qs_mv <- rdplot(
  y = indi_ns_ss2$post_interval6, x = indi_ns_ss2$scoringD2_0, binselect = "qsmv",
  x.lim = c(-0.25, 0.40), y.lim = c(0, 185), shade = F
) # quantile-spaced mv
summary(plotD2_6_qs_mv)

### Generate rdplot using ggplot2: preparing inputs
plot1 <- plotD2_6_qs_mv

# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss2$scoringD2_0
y <- indi_ns_ss2$post_interval6
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y[1], ": worked days during months [+1, +6]"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[2], " (", S[2], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y[1]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 185)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 0.75) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 0.75) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = c(-0.25, 0.4), ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/late_d2_y1.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


### Outcome 2: post [7, 12] months ---------------------
summary(indi_ns_ss2$scoringD2_0)
summary(indi_ns_ss2$post_interval12)


plotD2_12_qs_mv <- rdplot(
  y = indi_ns_ss2$post_interval712, x = indi_ns_ss2$scoringD2_0, binselect = "qsmv",
  x.lim = c(-0.25, 0.40), y.lim = c(0, 185), shade = F
) # quantile-spaced mv
summary(plotD2_12_qs_mv)


### Generate rdplot using ggplot2: preparing inputs
plot1 <- plotD2_12_qs_mv

# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss2$scoringD2_0
y <- indi_ns_ss2$post_interval712
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y[2], ": worked days during months [+7, +12]"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[2], " (", S[2], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y[2]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 185)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 0.75) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 0.75) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # Saved at 395x300

# Shade (alpha controls transparency)
temp_plot +
  geom_ribbon(aes(x = rdplot_mean_bin, ymin = rdplot_cil_bin, ymax = rdplot_cir_bin), alpha = 0.2) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/late_d2_y2.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


### Outcome 3: post 18 months ---------------------
summary(indi_ns_ss2$scoringD2_0)
summary(indi_ns_ss2$post_interval18)

plotD2_18_qs_imse <- rdplot(
  y = indi_ns_ss2$post_interval18, x = indi_ns_ss2$scoringD2_0, binselect = "qs",
  x.lim = c(-0.25, 0.40), y.lim = c(0, 575), ci = 95, hide = TRUE
) # quantile-spaced
summary(plotD2_18_qs_imse)

plotD2_18_qs_mv <- rdplot(
  y = indi_ns_ss2$post_interval1318, x = indi_ns_ss2$scoringD2_0, binselect = "qsmv",
  x.lim = c(-0.25, 0.40), y.lim = c(0, 185), shade = F
) # quantile-spaced mv
summary(plotD2_18_qs_mv)


### Generate rdplot using ggplot2: preparing inputs
plot1 <- plotD2_18_qs_mv

# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss2$scoringD2_0
y <- indi_ns_ss2$post_interval1318
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y[3], ": worked days during months [+13, +18]"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[2], " (", S[2], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y[3]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 185)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 0.75) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 0.75) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot # Saved at 395x300

# Shade (alpha controls transparency)
temp_plot +
  geom_ribbon(aes(x = rdplot_mean_bin, ymin = rdplot_cil_bin, ymax = rdplot_cir_bin), alpha = 0.2) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/late_d2_y3.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


### Outcome 4: post 24 months ---------------------
summary(indi_ns_ss2$scoringD2_0)
summary(indi_ns_ss2$post_interval24)

plotD2_24_qs_imse <- rdplot(
  y = indi_ns_ss2$post_interval24, x = indi_ns_ss2$scoringD2_0, binselect = "qs",
  x.lim = c(-0.25, 0.40), y.lim = c(0, 750), ci = 95, hide = TRUE
) # quantile-spaced
summary(plotD2_24_qs_imse)

plotD2_24_qs_mv <- rdplot(
  y = indi_ns_ss2$post_interval1924, x = indi_ns_ss2$scoringD2_0, binselect = "qsmv",
  x.lim = c(-0.25, 0.40), y.lim = c(0, 185), shade = F
) # quantile-spaced mv
summary(plotD2_24_qs_mv)


### Generate rdplot using ggplot2: preparing inputs
plot1 <- plotD2_24_qs_mv

# plot1 = rdplot(y,x, ci=95, hide=TRUE)
c <- 0
x <- indi_ns_ss2$scoringD2_0
y <- indi_ns_ss2$post_interval1924
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
rdplot_mean_y <- plot1$vars_bins[, "rdplot_mean_y"]
y_hat <- plot1$vars_poly[, "rdplot_y"]
x_plot <- plot1$vars_poly[, "rdplot_x"]
rdplot_cil_bin <- plot1$vars_bins[, "rdplot_ci_l"]
rdplot_cir_bin <- plot1$vars_bins[, "rdplot_ci_r"]
rdplot_mean_bin <- plot1$vars_bins[, "rdplot_mean_bin"]
y_hat_r <- y_hat[x_plot >= c]
y_hat_l <- y_hat[x_plot < c]
x_plot_r <- x_plot[x_plot >= c]
x_plot_l <- x_plot[x_plot < c]

col.lines <- "blue"
col.dots <- 1
type.dots <- 20
# title="Y1: worked days during months [+1, +6]"
title <- expression(paste(Y[4], ": worked days during months [+19, +24]"))
# x.label="Score for D1 (S1)"
x.label <- expression(paste("Score for ", D[2], " (", S[2], ")"))
# y.label="LSM of Y1"
y.label <- expression(paste("LSM of ", Y[4]))
x.lim <- c(min(x, na.rm = T), max(x, na.rm = T))
# y.lim=c(min(y, na.rm=T), max(y, na.rm=T))
y.lim <- c(0, 185)

### Generate rdplot using ggplot2: executing ggplot2

temp_plot <- ggplot() +
  theme_bw() +
  geom_point(aes(x = rdplot_mean_bin, y = rdplot_mean_y), col = "dimgrey", na.rm = TRUE, size = 1) +
  geom_line(aes(x = x_plot_l, y = y_hat_l), col = "red", na.rm = TRUE, size = 0.75) +
  geom_line(aes(x = x_plot_r, y = y_hat_r), col = "blue", na.rm = TRUE, size = 0.75) +
  labs(x = x.label, y = y.label) +
  ggtitle(title) +
  labs(title = title, y = y.label, x = x.label) +
  coord_cartesian(xlim = x.lim, ylim = y.lim) +
  theme(legend.position = "None") +
  geom_vline(xintercept = c, size = 0.5) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  )
temp_plot

# Shade (alpha controls transparency)
temp_plot +
  geom_ribbon(aes(x = rdplot_mean_bin, ymin = rdplot_cil_bin, ymax = rdplot_cir_bin), alpha = 0.2) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/late_d2_y4.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")


## Plots for logit models (employment quality) ---------
indi_ns_ss1_c$treatedD1 <- ifelse(indi_ns_ss1_c$scoringD1_0 > 0, TRUE, FALSE)
indi_ns_ss2_c$treatedD2 <- ifelse(indi_ns_ss2_c$scoringD2_0 > 0, TRUE, FALSE)

### D1---------
# D1 OEC
logit_oec <- binsglm(
  y = indi_ns_ss1_c$lcris_OEC,
  x = indi_ns_ss1_c$scoringD1_0,
  by = indi_ns_ss1_c$treatedD1,
  randcut = 1, family = binomial(link = "logit"),
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss1_c$scoringD1_0[indi_ns_ss1_c$treatedD1 == F]
outcome1d1_left <- indi_ns_ss1_c$lcris_OEC[indi_ns_ss1_c$treatedD1 == F]

running1_right <- indi_ns_ss1_c$scoringD1_0[indi_ns_ss1_c$treatedD1 == T]
outcome1d1_right <- indi_ns_ss1_c$lcris_OEC[indi_ns_ss1_c$treatedD1 == T]

xlabeld1 <- expression(paste("Score for ", D[1], " (", S[1], ")"))
ylabelo1 <- expression(paste("Prop. of ", Y[OEC], " = YES"))
titleo1 <- expression(paste(Y[OEC], " = YES in [+1, +24]"))

xlimite <- c(min(indi_ns_ss1_c$scoringD1_0, na.rm = T), max(indi_ns_ss1_c$scoringD1_0, na.rm = T))
ylimite <- c(0, 1)

dots_logit_oec <- rbind(logit_oec[["data.plot"]][["Group FALSE"]][["data.dots"]],
                        logit_oec[["data.plot"]][["Group TRUE"]][["data.dots"]])

temp_plot <- ggplot() +
  geom_point(data=dots_logit_oec, aes(x=x, y=fit), color="dimgrey", size=2) +
  geom_smooth(aes(running1_left, outcome1d1_left),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "red", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_smooth(aes(running1_right, outcome1d1_right),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "blue", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo1, y = ylabelo1, x = xlabeld1) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) +
  theme_bw() # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/latej_d1_yoec.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# D1 FTCm12
logit_FTCm12 <- binsglm(
  y = indi_ns_ss1_c$lcris_FTCm12,
  x = indi_ns_ss1_c$scoringD1_0,
  by = indi_ns_ss1_c$treatedD1,
  randcut = 1, family = binomial(link = "logit"),
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss1_c$scoringD1_0[indi_ns_ss1_c$treatedD1 == F]
outcome1d1_left <- indi_ns_ss1_c$lcris_FTCm12[indi_ns_ss1_c$treatedD1 == F]

running1_right <- indi_ns_ss1_c$scoringD1_0[indi_ns_ss1_c$treatedD1 == T]
outcome1d1_right <- indi_ns_ss1_c$lcris_FTCm12[indi_ns_ss1_c$treatedD1 == T]

xlabeld1 <- expression(paste("Score for ", D[1], " (", S[1], ")"))
ylabelo1 <- expression(paste("Prop. of ", Y[FTCm12], " = YES"))
titleo1 <- expression(paste(Y[FTCm12], " = YES in [+1, +24]"))

xlimite <- c(min(indi_ns_ss1_c$scoringD1_0, na.rm = T), max(indi_ns_ss1_c$scoringD1_0, na.rm = T))
ylimite <- c(0, 1)

dots_logit_FTCm12 <- rbind(logit_FTCm12[["data.plot"]][["Group FALSE"]][["data.dots"]],
                           logit_FTCm12[["data.plot"]][["Group TRUE"]][["data.dots"]])

temp_plot <- ggplot() +
  geom_point(data=dots_logit_FTCm12, aes(x=x, y=fit), color="dimgrey", size=2) +
  geom_smooth(aes(running1_left, outcome1d1_left),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "red", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_smooth(aes(running1_right, outcome1d1_right),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "blue", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo1, y = ylabelo1, x = xlabeld1) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) +
  theme_bw() # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/latej_d1_yftcm12.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# D1 FTC612
logit_FTC612 <- binsglm(
  y = indi_ns_ss1_c$lcris_FTC612,
  x = indi_ns_ss1_c$scoringD1_0,
  by = indi_ns_ss1_c$treatedD1,
  randcut = 1, family = binomial(link = "logit"),
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss1_c$scoringD1_0[indi_ns_ss1_c$treatedD1 == F]
outcome1d1_left <- indi_ns_ss1_c$lcris_FTC612[indi_ns_ss1_c$treatedD1 == F]

running1_right <- indi_ns_ss1_c$scoringD1_0[indi_ns_ss1_c$treatedD1 == T]
outcome1d1_right <- indi_ns_ss1_c$lcris_FTC612[indi_ns_ss1_c$treatedD1 == T]

xlabeld1 <- expression(paste("Score for ", D[1], " (", S[1], ")"))
ylabelo1 <- expression(paste("Prop. of ", Y[FTC612], " = YES"))
titleo1 <- expression(paste(Y[FTC612], " = YES in [+1, +24]"))

xlimite <- c(min(indi_ns_ss1_c$scoringD1_0, na.rm = T), max(indi_ns_ss1_c$scoringD1_0, na.rm = T))
ylimite <- c(0, 1)

dots_logit_FTC612 <- rbind(logit_FTC612[["data.plot"]][["Group FALSE"]][["data.dots"]],
                           logit_FTC612[["data.plot"]][["Group TRUE"]][["data.dots"]])

temp_plot <- ggplot() +
  geom_point(data=dots_logit_FTC612, aes(x=x, y=fit), color="dimgrey", size=2) +
  geom_smooth(aes(running1_left, outcome1d1_left),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "red", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_smooth(aes(running1_right, outcome1d1_right),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "blue", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo1, y = ylabelo1, x = xlabeld1) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) +
  theme_bw() # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/latej_d1_yftc612.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# D1 FTC05
logit_FTC05 <- binsglm(
  y = indi_ns_ss1_c$lcris_FTCme6,
  x = indi_ns_ss1_c$scoringD1_0,
  by = indi_ns_ss1_c$treatedD1,
  randcut = 1, family = binomial(link = "logit"),
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss1_c$scoringD1_0[indi_ns_ss1_c$treatedD1 == F]
outcome1d1_left <- indi_ns_ss1_c$lcris_FTCme6[indi_ns_ss1_c$treatedD1 == F]

running1_right <- indi_ns_ss1_c$scoringD1_0[indi_ns_ss1_c$treatedD1 == T]
outcome1d1_right <- indi_ns_ss1_c$lcris_FTCme6[indi_ns_ss1_c$treatedD1 == T]

xlabeld1 <- expression(paste("Score for ", D[1], " (", S[1], ")"))
ylabelo1 <- expression(paste("Prop. of ", Y[FTC05], " = YES"))
titleo1 <- expression(paste(Y[FTC05], " = YES in [+1, +24]"))

xlimite <- c(min(indi_ns_ss1_c$scoringD1_0, na.rm = T), max(indi_ns_ss1_c$scoringD1_0, na.rm = T))
ylimite <- c(0, 1)

dots_logit_FTC05 <- rbind(logit_FTC05[["data.plot"]][["Group FALSE"]][["data.dots"]],
                          logit_FTC05[["data.plot"]][["Group TRUE"]][["data.dots"]])

temp_plot <- ggplot() +
  geom_point(data=dots_logit_FTC05, aes(x=x, y=fit), color="dimgrey", size=2) +
  geom_smooth(aes(running1_left, outcome1d1_left),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "red", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_smooth(aes(running1_right, outcome1d1_right),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "blue", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo1, y = ylabelo1, x = xlabeld1) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) +
  theme_bw() # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/latej_d1_yftc05.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

### D2---------
# D2 OEC
logit_oec <- binsglm(
  y = indi_ns_ss2_c$lcris_OEC,
  x = indi_ns_ss2_c$scoringD2_0,
  by = indi_ns_ss2_c$treatedD2,
  randcut = 1, family = binomial(link = "logit"),
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss2_c$scoringD2_0[indi_ns_ss2_c$treatedD2 == F]
outcome1D2_left <- indi_ns_ss2_c$lcris_OEC[indi_ns_ss2_c$treatedD2 == F]

running1_right <- indi_ns_ss2_c$scoringD2_0[indi_ns_ss2_c$treatedD2 == T]
outcome1D2_right <- indi_ns_ss2_c$lcris_OEC[indi_ns_ss2_c$treatedD2 == T]

xlabelD2 <- expression(paste("Score for ", D[2], " (", S[2], ")"))
ylabelo1 <- expression(paste("Prop. of ", Y[OEC], " = YES"))
titleo1 <- expression(paste(Y[OEC], " = YES in [+1, +24]"))

xlimite <- c(min(indi_ns_ss2_c$scoringD2_0, na.rm = T), max(indi_ns_ss2_c$scoringD2_0, na.rm = T))
ylimite <- c(0, 1)

dots_logit_oec <- rbind(logit_oec[["data.plot"]][["Group FALSE"]][["data.dots"]],
                        logit_oec[["data.plot"]][["Group TRUE"]][["data.dots"]])

temp_plot <- ggplot() +
  geom_point(data=dots_logit_oec, aes(x=x, y=fit), color="dimgrey", size=2) +
  geom_smooth(aes(running1_left, outcome1D2_left),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "red", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_smooth(aes(running1_right, outcome1D2_right),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "blue", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo1, y = ylabelo1, x = xlabelD2) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) +
  theme_bw() # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/latej_d2_yoec.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# D2 FTCm12
logit_FTCm12 <- binsglm(
  y = indi_ns_ss2_c$lcris_FTCm12,
  x = indi_ns_ss2_c$scoringD2_0,
  by = indi_ns_ss2_c$treatedD2,
  randcut = 1, family = binomial(link = "logit"),
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss2_c$scoringD2_0[indi_ns_ss2_c$treatedD2 == F]
outcome1D2_left <- indi_ns_ss2_c$lcris_FTCm12[indi_ns_ss2_c$treatedD2 == F]

running1_right <- indi_ns_ss2_c$scoringD2_0[indi_ns_ss2_c$treatedD2 == T]
outcome1D2_right <- indi_ns_ss2_c$lcris_FTCm12[indi_ns_ss2_c$treatedD2 == T]

xlabelD2 <- expression(paste("Score for ", D[2], " (", S[2], ")"))
ylabelo1 <- expression(paste("Prop. of ", Y[FTCm12], " = YES"))
titleo1 <- expression(paste(Y[FTCm12], " = YES in [+1, +24]"))

xlimite <- c(min(indi_ns_ss2_c$scoringD2_0, na.rm = T), max(indi_ns_ss2_c$scoringD2_0, na.rm = T))
ylimite <- c(0, 1)

dots_logit_FTCm12 <- rbind(logit_FTCm12[["data.plot"]][["Group FALSE"]][["data.dots"]],
                           logit_FTCm12[["data.plot"]][["Group TRUE"]][["data.dots"]])

temp_plot <- ggplot() +
  geom_point(data=dots_logit_FTCm12, aes(x=x, y=fit), color="dimgrey", size=2) +
  geom_smooth(aes(running1_left, outcome1D2_left),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "red", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_smooth(aes(running1_right, outcome1D2_right),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "blue", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo1, y = ylabelo1, x = xlabelD2) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) +
  theme_bw() # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/latej_d2_yftcm12.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# D2 FTC612
logit_FTC612 <- binsglm(
  y = indi_ns_ss2_c$lcris_FTC612,
  x = indi_ns_ss2_c$scoringD2_0,
  by = indi_ns_ss2_c$treatedD2,
  randcut = 1, family = binomial(link = "logit"),
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss2_c$scoringD2_0[indi_ns_ss2_c$treatedD2 == F]
outcome1D2_left <- indi_ns_ss2_c$lcris_FTC612[indi_ns_ss2_c$treatedD2 == F]

running1_right <- indi_ns_ss2_c$scoringD2_0[indi_ns_ss2_c$treatedD2 == T]
outcome1D2_right <- indi_ns_ss2_c$lcris_FTC612[indi_ns_ss2_c$treatedD2 == T]

xlabelD2 <- expression(paste("Score for ", D[2], " (", S[2], ")"))
ylabelo1 <- expression(paste("Prop. of ", Y[FTC612], " = YES"))
titleo1 <- expression(paste(Y[FTC612], " = YES in [+1, +24]"))

xlimite <- c(min(indi_ns_ss2_c$scoringD2_0, na.rm = T), max(indi_ns_ss2_c$scoringD2_0, na.rm = T))
ylimite <- c(0, 1)


dots_logit_FTC612 <- rbind(logit_FTC612[["data.plot"]][["Group FALSE"]][["data.dots"]],
                           logit_FTC612[["data.plot"]][["Group TRUE"]][["data.dots"]])

temp_plot <- ggplot() +
  geom_point(data=dots_logit_FTC612, aes(x=x, y=fit), color="dimgrey", size=2) +
  geom_smooth(aes(running1_left, outcome1D2_left),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "red", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_smooth(aes(running1_right, outcome1D2_right),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "blue", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo1, y = ylabelo1, x = xlabelD2) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) +
  theme_bw() # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/latej_d2_yftc612.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# D2 FTC05
logit_FTC05 <- binsglm(
  y = indi_ns_ss2_c$lcris_FTCme6,
  x = indi_ns_ss2_c$scoringD2_0,
  by = indi_ns_ss2_c$treatedD2,
  randcut = 1, family = binomial(link = "logit"),
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss2_c$scoringD2_0[indi_ns_ss2_c$treatedD2 == F]
outcome1D2_left <- indi_ns_ss2_c$lcris_FTCme6[indi_ns_ss2_c$treatedD2 == F]

running1_right <- indi_ns_ss2_c$scoringD2_0[indi_ns_ss2_c$treatedD2 == T]
outcome1D2_right <- indi_ns_ss2_c$lcris_FTCme6[indi_ns_ss2_c$treatedD2 == T]

xlabelD2 <- expression(paste("Score for ", D[2], " (", S[2], ")"))
ylabelo1 <- expression(paste("Prop. of ", Y[FTC05], " = YES"))
titleo1 <- expression(paste(Y[FTC05], " = YES in [+1, +24]"))

xlimite <- c(min(indi_ns_ss2_c$scoringD2_0, na.rm = T), max(indi_ns_ss2_c$scoringD2_0, na.rm = T))
ylimite <- c(0, 1)


dots_logit_FTC05 <- rbind(logit_FTC05[["data.plot"]][["Group FALSE"]][["data.dots"]],
                          logit_FTC05[["data.plot"]][["Group TRUE"]][["data.dots"]])

temp_plot <- ggplot() +
  geom_point(data=dots_logit_FTC05, aes(x=x, y=fit), color="dimgrey", size=2) +
  geom_smooth(aes(running1_left, outcome1D2_left),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "red", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_smooth(aes(running1_right, outcome1D2_right),
              method = "glm", method.args = list(family = "binomial"),
              se = F, col = "blue", linewidth = 0.75,
              formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo1, y = ylabelo1, x = xlabelD2) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) +
  theme_bw() # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/latej_d2_yftc05.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

## Plots for quantile effects -------------
indi_ns_ss1$treatedD1 <- ifelse(indi_ns_ss1$scoringD1_0 > 0, TRUE, FALSE)
indi_ns_ss2$treatedD2 <- ifelse(indi_ns_ss2$scoringD2_0 > 0, TRUE, FALSE)

### D1 -------------
# o1
estbs1_y1 <- binsqreg(
  y = indi_ns_ss1$post_interval6,
  x = indi_ns_ss1$scoringD1_0,
  by = indi_ns_ss1$treatedD1,
  quantile = 0.5, randcut = 1,
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss1$scoringD1_0[indi_ns_ss1$treatedD1 == F]
outcome1d1_left <- indi_ns_ss1$post_interval6[indi_ns_ss1$treatedD1 == F]

running1_right <- indi_ns_ss1$scoringD1_0[indi_ns_ss1$treatedD1 == T]
outcome1d1_right <- indi_ns_ss1$post_interval6[indi_ns_ss1$treatedD1 == T]

xlabeld1 <- expression(paste("Score for ", D[1], " (", S[1], ")"))
ylabelo1 <- expression(paste("LSMed of ", Y[1]))
titleo1 <- expression(paste(Y[1], ": worked days during months [+1, +6]"))

xlimite <- c(min(indi_ns_ss1$scoringD1_0, na.rm = T), max(indi_ns_ss1$scoringD1_0, na.rm = T))
ylimite <- c(0, 185)

temp_plot <- estbs1_y1$bins_plot +
  geom_quantile(aes(running1_left, outcome1d1_left),
                quantiles = 0.5, col = "red", linewidth = 0.75,
                formula = y ~ poly(x, 4)) +
  geom_quantile(aes(running1_right, outcome1d1_right),
                quantiles = 0.5, col = "blue", linewidth = 0.75,
                formula = y ~ poly(x, 4)) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo1, y = ylabelo1, x = xlabeld1) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/lqte_d1_y1.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# o2
estbs1_y2 <- binsqreg(
  y = indi_ns_ss1$post_interval712,
  x = indi_ns_ss1$scoringD1_0,
  by = indi_ns_ss1$treatedD1,
  quantile = 0.5, randcut = 1,
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss1$scoringD1_0[indi_ns_ss1$treatedD1 == F]
outcome2d1_left <- indi_ns_ss1$post_interval712[indi_ns_ss1$treatedD1 == F]

running1_right <- indi_ns_ss1$scoringD1_0[indi_ns_ss1$treatedD1 == T]
outcome2d1_right <- indi_ns_ss1$post_interval712[indi_ns_ss1$treatedD1 == T]

xlabeld1 <- expression(paste("Score for ", D[1], " (", S[1], ")"))
ylabelo2 <- expression(paste("LSMed of ", Y[2]))
titleo2 <- expression(paste(Y[2], ": worked days during months [+7, +12]"))

xlimite <- c(min(indi_ns_ss1$scoringD1_0, na.rm = T), max(indi_ns_ss1$scoringD1_0, na.rm = T))
ylimite <- c(0, 185)

temp_plot <- estbs1_y2$bins_plot +
  geom_quantile(aes(running1_left, outcome2d1_left),
                quantiles = 0.5, col = "red", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_quantile(aes(running1_right, outcome2d1_right),
                quantiles = 0.5, col = "blue", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo2, y = ylabelo2, x = xlabeld1) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/lqte_d1_y2.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# o3
estbs1_y3 <- binsqreg(
  y = indi_ns_ss1$post_interval1318,
  x = indi_ns_ss1$scoringD1_0,
  by = indi_ns_ss1$treatedD1,
  quantile = 0.5, randcut = 1,
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss1$scoringD1_0[indi_ns_ss1$treatedD1 == F]
outcome3d1_left <- indi_ns_ss1$post_interval1318[indi_ns_ss1$treatedD1 == F]

running1_right <- indi_ns_ss1$scoringD1_0[indi_ns_ss1$treatedD1 == T]
outcome3d1_right <- indi_ns_ss1$post_interval1318[indi_ns_ss1$treatedD1 == T]

xlabeld1 <- expression(paste("Score for ", D[1], " (", S[1], ")"))
ylabelo3 <- expression(paste("LSMed of ", Y[3]))
titleo3 <- expression(paste(Y[3], ": worked days during months [+13, +18]"))

xlimite <- c(min(indi_ns_ss1$scoringD1_0, na.rm = T), max(indi_ns_ss1$scoringD1_0, na.rm = T))
ylimite <- c(0, 185)

temp_plot <- estbs1_y3$bins_plot +
  geom_quantile(aes(running1_left, outcome3d1_left),
                quantiles = 0.5, col = "red", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_quantile(aes(running1_right, outcome3d1_right),
                quantiles = 0.5, col = "blue", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo3, y = ylabelo3, x = xlabeld1) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/lqte_d1_y3.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# o4
estbs1_y4 <- binsqreg(
  y = indi_ns_ss1$post_interval1924,
  x = indi_ns_ss1$scoringD1_0,
  by = indi_ns_ss1$treatedD1,
  quantile = 0.5, randcut = 1,
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running1_left <- indi_ns_ss1$scoringD1_0[indi_ns_ss1$treatedD1 == F]
outcome4d1_left <- indi_ns_ss1$post_interval1924[indi_ns_ss1$treatedD1 == F]

running1_right <- indi_ns_ss1$scoringD1_0[indi_ns_ss1$treatedD1 == T]
outcome4d1_right <- indi_ns_ss1$post_interval1924[indi_ns_ss1$treatedD1 == T]

xlabeld1 <- expression(paste("Score for ", D[1], " (", S[1], ")"))
ylabelo4 <- expression(paste("LSMed of ", Y[4]))
titleo4 <- expression(paste(Y[4], ": worked days during months [+19, +24]"))

ylim <- c(0, 185)

temp_plot <- estbs1_y4$bins_plot +
  geom_quantile(aes(running1_left, outcome4d1_left),
                quantiles = 0.5, col = "red", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_quantile(aes(running1_right, outcome4d1_right),
                quantiles = 0.5, col = "blue", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(ylim = ylim) +
  labs(title = titleo4, y = ylabelo4, x = xlabeld1) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/lqte_d1_y4.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

### D2 -------------
# o1
estbs2_y1 <- binsqreg(
  y = indi_ns_ss2$post_interval6,
  x = indi_ns_ss2$scoringD2_0,
  by = indi_ns_ss2$treatedD2,
  quantile = 0.5, randcut = 1,
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

running2_left <- indi_ns_ss2$scoringD2_0[indi_ns_ss2$treatedD2 == F]
outcome1d2_left <- indi_ns_ss2$post_interval6[indi_ns_ss2$treatedD2 == F]

running2_right <- indi_ns_ss2$scoringD2_0[indi_ns_ss2$treatedD2 == T]
outcome1d2_right <- indi_ns_ss2$post_interval6[indi_ns_ss2$treatedD2 == T]

xlabeld2 <- expression(paste("Score for ", D[2], " (", S[2], ")"))
ylabelo1 <- expression(paste("LSMed of ", Y[1]))
titleo1 <- expression(paste(Y[1], ": worked days during months [+1, +6]"))

xlimite <- c(min(indi_ns_ss2$scoringD2_0, na.rm = T), max(indi_ns_ss2$scoringD2_0, na.rm = T))
ylimite <- c(0, 185)

temp_plot <- estbs2_y1$bins_plot +
  geom_quantile(aes(running2_left, outcome1d2_left),
                quantiles = 0.5, col = "red", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_quantile(aes(running2_right, outcome1d2_right),
                quantiles = 0.5, col = "blue", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo1, y = ylabelo1, x = xlabeld2) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/lqte_d2_y1.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# o2
estbs2_y2 <- binsqreg(
  y = indi_ns_ss2$post_interval712,
  x = indi_ns_ss2$scoringD2_0,
  by = indi_ns_ss2$treatedD2,
  quantile = 0.5, randcut = 1,
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

# running2_left <- indi_ns_ss2$scoringD2_0[indi_ns_ss2$treatedD2 == F]
outcome2d2_left <- indi_ns_ss2$post_interval712[indi_ns_ss2$treatedD2 == F]

# running2_right <- indi_ns_ss2$scoringD2_0[indi_ns_ss2$treatedD2 == T]
outcome2d2_right <- indi_ns_ss2$post_interval712[indi_ns_ss2$treatedD2 == T]

# xlabeld1 <- expression(paste("Score for ", D[1]," (", S[1], ")"))
ylabelo2 <- expression(paste("LSMed of ", Y[2]))
titleo2 <- expression(paste(Y[2], ": worked days during months [+7, +12]"))

xlimite <- c(min(indi_ns_ss2$scoringD2_0, na.rm = T), max(indi_ns_ss2$scoringD2_0, na.rm = T))
ylimite <- c(0, 185)

temp_plot <- estbs2_y2$bins_plot +
  geom_quantile(aes(running2_left, outcome2d2_left),
                quantiles = 0.5, col = "red", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_quantile(aes(running2_right, outcome2d2_right),
                quantiles = 0.5, col = "blue", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo2, y = ylabelo2, x = xlabeld2) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/lqte_d2_y2.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# o3
estbs2_y3 <- binsqreg(
  y = indi_ns_ss2$post_interval1318,
  x = indi_ns_ss2$scoringD2_0,
  by = indi_ns_ss2$treatedD2,
  quantile = 0.5, randcut = 1,
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

# running1_left <- indi_ns_ss2$scoringD2_0[indi_ns_ss2$treatedD2 == F]
outcome3d2_left <- indi_ns_ss2$post_interval1318[indi_ns_ss2$treatedD2 == F]

# running1_right <- indi_ns_ss2$scoringD2_0[indi_ns_ss2$treatedD2 == T]
outcome3d2_right <- indi_ns_ss2$post_interval1318[indi_ns_ss2$treatedD2 == T]

# xlabeld1 <- expression(paste("Score for ", D[1]," (", S[1], ")"))
ylabelo3 <- expression(paste("LSMed of ", Y[3]))
titleo3 <- expression(paste(Y[3], ": worked days during months [+13, +18]"))

xlimite <- c(min(indi_ns_ss2$scoringD2_0, na.rm = T), max(indi_ns_ss2$scoringD2_0, na.rm = T))
ylimite <- c(0, 185)

temp_plot <- estbs2_y3$bins_plot +
  geom_quantile(aes(running2_left, outcome3d2_left),
                quantiles = 0.5, col = "red", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_quantile(aes(running2_right, outcome3d2_right),
                quantiles = 0.5, col = "blue", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(xlim = xlimite, ylim = ylimite) +
  labs(title = titleo3, y = ylabelo3, x = xlabeld2) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/lqte_d2_y3.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# o4
estbs2_y4 <- binsqreg(
  y = indi_ns_ss2$post_interval1924,
  x = indi_ns_ss2$scoringD2_0,
  by = indi_ns_ss2$treatedD2,
  quantile = 0.5, randcut = 1,
  deriv = 0,
  bycolors = c("dimgrey", "dimgrey"),
  legendoff = T
)

# running1_left <- indi_ns_ss2$scoringD2_0[indi_ns_ss2$treatedD2 == F]
outcome4d2_left <- indi_ns_ss2$post_interval1924[indi_ns_ss2$treatedD2 == F]

# running1_right <- indi_ns_ss2$scoringD2_0[indi_ns_ss2$treatedD2 == T]
outcome4d2_right <- indi_ns_ss2$post_interval1924[indi_ns_ss2$treatedD2 == T]

# xlabeld1 <- expression(paste("Score for ", D[1]," (", S[1], ")"))
ylabelo4 <- expression(paste("LSMed of ", Y[4]))
titleo4 <- expression(paste(Y[4], ": worked days during months [+19, +24]"))

ylim <- c(0, 185)

temp_plot <- estbs2_y4$bins_plot +
  geom_quantile(aes(running2_left, outcome4d2_left),
                quantiles = 0.5, col = "red", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_quantile(aes(running2_right, outcome4d2_right),
                quantiles = 0.5, col = "blue", linewidth = 0.75,
                formula = y ~ poly(x, 4)
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  coord_cartesian(ylim = ylim) +
  labs(title = titleo4, y = ylabelo4, x = xlabeld2) +
  theme(
    axis.text.x = element_text(size = rel(1.2)),
    axis.text.y = element_text(size = rel(1.2))
  ) # Saved at 395x300

ggsave(filename = "intermediate/script04/plots/lqte_d2_y4.svg",  plot = temp_plot,
       width = 395 / 96, height = 300 / 96, dpi = 96, # Resolución de 96 dpi (píxeles por pulgada)
       units = "in", device = "svg")

# 2. LATE heterogeneity (subgroup analyses) ----------------
#saveRDS(indi_ns_ss1_c, "indi_ns_ss1_c_24-02-23.RDS")
#saveRDS(indi_ns_ss2_c, "indi_ns_ss2_c_24-02-23.RDS")

#saveRDS(indi_ns_ss1, "indi_ns_ss1_24-02-23.RDS")
#saveRDS(indi_ns_ss2, "indi_ns_ss2_24-02-23.RDS")

#indi_ns_ss1 <- readRDS("indi_ns_ss1_24-02-23.RDS")
#indi_ns_ss2 <- readRDS("indi_ns_ss2_24-02-23.RDS")

#indi_ns_ss1_c <- readRDS("indi_ns_ss1_c_24-02-23.RDS")
#indi_ns_ss2_c <- readRDS("indi_ns_ss2_c_24-02-23.RDS")


## Depending on sex ------------
### D1 ---------------
indi_ns_ss1_F <- indi_ns_ss1_c %>% filter(genere == "F")
indi_ns_ss1_M <- indi_ns_ss1_c %>% filter(genere == "M")

### Outcome 1 for women: days worked post months [1, 6] --> Unsignificant effect (unsignificant difference)
rddD1ei6_p1kTw <- rdrobust(y = indi_ns_ss1_F$post_interval6, x = indi_ns_ss1_F$scoringD1_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd", cluster = NULL)
summary(rddD1ei6_p1kTw)


### Outcome 2 for women: days worked post months [7, 12] --> Unsignificant effect (unsignificant difference)
rddD1ei12_p1kTw <- rdrobust(y = indi_ns_ss1_F$post_interval712, x = indi_ns_ss1_F$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD1ei12_p1kTw)


### Outcome 3 for women: days worked post [13, 18] months --> Unsignificant effect (unsignificant difference)
rddD1ei18_p1kTw <- rdrobust(y = indi_ns_ss1_F$post_interval1318, x = indi_ns_ss1_F$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD1ei18_p1kTw)

### Outcome 4 for women: days worked post [19, 24] months --> Unsignificant effect (unsignificant difference)
rddD1ei24_p1kTw <- rdrobust(y = indi_ns_ss1_F$post_interval1924, x = indi_ns_ss1_F$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD1ei24_p1kTw)

### Outcome 1 for men: days worked post months [1, 6] --> Unsignificant effect (unsignificant difference)
rddD1ei6_p1kTm <- rdrobust(y = indi_ns_ss1_M$post_interval6, x = indi_ns_ss1_M$scoringD1_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd", cluster = NULL)
summary(rddD1ei6_p1kTm)


### Outcome 2 for men: days worked post months [7, 12] --> Unsignificant effect (unsignificant difference)
rddD1ei12_p1kTm <- rdrobust(y = indi_ns_ss1_M$post_interval712, x = indi_ns_ss1_M$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD1ei12_p1kTm)

### Outcome 3 for men: days worked post [13, 18] months --> Unsignificant effect (unsignificant difference)
rddD1ei18_p1kTm <- rdrobust(y = indi_ns_ss1_M$post_interval1318, x = indi_ns_ss1_M$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD1ei18_p1kTm)

### Outcome 4 for men: days worked post [19, 24] months --> Unsignificant effect (unsignificant difference)
rddD1ei24_p1kTm <- rdrobust(y = indi_ns_ss1_M$post_interval1924, x = indi_ns_ss1_M$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD1ei24_p1kTm)






### D2 ------------
indi_ns_ss2_F <- subset(indi_ns_ss2_c, indi_ns_ss2_c$genere == "F")
indi_ns_ss2_M <- subset(indi_ns_ss2_c, indi_ns_ss2_c$genere == "M")


### Outcome 1 for women: days worked post months [1, 6] --> Unsignificant effect (unsignificant difference)
rddD2ei6_p1kTw <- rdrobust(y = indi_ns_ss2_F$post_interval6, x = indi_ns_ss2_F$scoringD2_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd", cluster = NULL)
summary(rddD2ei6_p1kTw)


### Outcome 2 for women: days worked post months [7, 12] --> Unsignificant effect (unsignificant difference)
rddD2ei12_p1kTw <- rdrobust(y = indi_ns_ss2_F$post_interval712, x = indi_ns_ss2_F$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD2ei12_p1kTw)


### Outcome 3 for women: days worked post [13, 18] months --> Unsignificant effect (unsignificant difference)
rddD2ei18_p1kTw <- rdrobust(y = indi_ns_ss2_F$post_interval1318, x = indi_ns_ss2_F$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD2ei18_p1kTw)

### Outcome 4 for women: days worked post [19, 24] months --> Unsignificant effect (unsignificant difference)
rddD2ei24_p1kTw <- rdrobust(y = indi_ns_ss2_F$post_interval1924, x = indi_ns_ss2_F$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD2ei24_p1kTw)

### Outcome 1 for men: days worked post months [1, 6] --> Unsignificant effect (unsignificant difference)
rddD2ei6_p1kTm <- rdrobust(y = indi_ns_ss2_M$post_interval6, x = indi_ns_ss2_M$scoringD2_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd", cluster = NULL)
summary(rddD2ei6_p1kTm)


### Outcome 2 for men: days worked post months [7, 12] --> Unsignificant effect (unsignificant difference)
rddD2ei12_p1kTm <- rdrobust(y = indi_ns_ss2_M$post_interval712, x = indi_ns_ss2_M$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD2ei12_p1kTm)

### Outcome 3 for men: days worked post [13, 18] months --> Unsignificant effect (unsignificant difference)
rddD2ei18_p1kTm <- rdrobust(y = indi_ns_ss2_M$post_interval1318, x = indi_ns_ss2_M$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD2ei18_p1kTm)

### Outcome 4 for men: days worked post [19, 24] months --> Unsignificant effect (unsignificant difference)
rddD2ei24_p1kTm <- rdrobust(y = indi_ns_ss2_M$post_interval1924, x = indi_ns_ss2_M$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)
summary(rddD2ei24_p1kTm)


## Depending on age --------------
### D1 ---------------
indi_ns_ss1_ageup <- indi_ns_ss1_c %>% filter(eta > 50 | eta == 50)
indi_ns_ss1_agedown <- indi_ns_ss1_c %>% filter(eta < 50)


### Outcome 1 for age above 50
rddD1ei6_p1kTu <- rdrobust(y = indi_ns_ss1_ageup$post_interval6, x = indi_ns_ss1_ageup$scoringD1_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd", cluster = NULL)

### Outcome 2 for age above50
rddD1ei12_p1kTu <- rdrobust(y = indi_ns_ss1_ageup$post_interval712, x = indi_ns_ss1_ageup$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 3 for age above 50
rddD1ei18_p1kTu <- rdrobust(y = indi_ns_ss1_ageup$post_interval1318, x = indi_ns_ss1_ageup$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 4 for age above 50
rddD1ei24_p1kTu <- rdrobust(y = indi_ns_ss1_ageup$post_interval1924, x = indi_ns_ss1_ageup$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)


### Outcome 1 for age below 50
rddD1ei6_p1kTd <- rdrobust(y = indi_ns_ss1_agedown$post_interval6, x = indi_ns_ss1_agedown$scoringD1_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd", cluster = NULL)

### Outcome 2 for age below 50
rddD1ei12_p1kTd <- rdrobust(y = indi_ns_ss1_agedown$post_interval712, x = indi_ns_ss1_agedown$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 3 for age below 50
rddD1ei18_p1kTd <- rdrobust(y = indi_ns_ss1_agedown$post_interval1318, x = indi_ns_ss1_agedown$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 4 for age below 50
rddD1ei24_p1kTd <- rdrobust(y = indi_ns_ss1_agedown$post_interval1924, x = indi_ns_ss1_agedown$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)





### D2 ---------------
indi_ns_ss2_ageup <- subset(indi_ns_ss2_c, indi_ns_ss2_c$eta > 50 | indi_ns_ss2_c$eta == 50)
indi_ns_ss2_agedown <- subset(indi_ns_ss2_c, indi_ns_ss2_c$eta < 50)


### Outcome 1 for age above 50
rddD2ei6_p1kTu <- rdrobust(y = indi_ns_ss2_ageup$post_interval6, x = indi_ns_ss2_ageup$scoringD2_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd", cluster = NULL)

### Outcome 2 for age above50
rddD2ei12_p1kTu <- rdrobust(y = indi_ns_ss2_ageup$post_interval712, x = indi_ns_ss2_ageup$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 3 for age above 50
rddD2ei18_p1kTu <- rdrobust(y = indi_ns_ss2_ageup$post_interval1318, x = indi_ns_ss2_ageup$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 4 for age above 50
rddD2ei24_p1kTu <- rdrobust(y = indi_ns_ss2_ageup$post_interval1924, x = indi_ns_ss2_ageup$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)


### Outcome 1 for age below 50
rddD2ei6_p1kTd <- rdrobust(y = indi_ns_ss2_agedown$post_interval6, x = indi_ns_ss2_agedown$scoringD2_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd", cluster = NULL)

### Outcome 2 for age below 50
rddD2ei12_p1kTd <- rdrobust(y = indi_ns_ss2_agedown$post_interval712, x = indi_ns_ss2_agedown$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 3 for age below 50
rddD2ei18_p1kTd <- rdrobust(y = indi_ns_ss2_agedown$post_interval1318, x = indi_ns_ss2_agedown$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 4 for age below 50
rddD2ei24_p1kTd <- rdrobust(y = indi_ns_ss2_agedown$post_interval1924, x = indi_ns_ss2_agedown$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)


## Depending on education ------------
### D1--------------
indi_ns_ss1_1_PriLS <- indi_ns_ss1_c %>%
  filter(studio2_grouped2 == "1_primary" | studio2_grouped2 == "2_lowersec")

indi_ns_ss1_1_UppSec <- indi_ns_ss1_c %>%
  filter(studio2_grouped2 %in% c("344_uppersecG", "353_uppersecP_3y", "353_uppersecP_4y_prof",
                                 "353_uppersecP_4y_tec"))

indi_ns_ss1_1_Uni <- indi_ns_ss1_c %>%
  filter(studio2_grouped2 %in% c("660_terBachy+_lev1", "760+_terMastery+_lev2"))


## PriLS, UppSec, Uni


### Outcome 1 for PriLS: days worked post months [1, 6] --> Not significant effect (not significant difference) !!!!!
rddD1ei6_p1kT_PriLS <- rdrobust(y = indi_ns_ss1_1_PriLS$post_interval6, x = indi_ns_ss1_1_PriLS$scoringD1_0,
                                kernel = "triangular",
                                c = 0, p = 1, bwselect = "mserd",
                                cluster = NULL)

### Outcome 2 for PriLS: days worked post months [7, 12] --> Not significant effect (not ssignificant difference) !!!!!!!!!!!!!!!!
rddD1ei12_p1kT_PriLS <- rdrobust(y = indi_ns_ss1_1_PriLS$post_interval712, x = indi_ns_ss1_1_PriLS$scoringD1_0,
                                 kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                 cluster = NULL)

### Outcome 3 for PriLS: days worked post [13, 18] months
rddD1ei18_p1kT_PriLS <- rdrobust(y = indi_ns_ss1_1_PriLS$post_interval1318, x = indi_ns_ss1_1_PriLS$scoringD1_0,
                                 kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                 cluster = NULL)

### Outcome 4 for PriLS: days worked post [19, 24] months
rddD1ei24_p1kT_PriLS <- rdrobust(y = indi_ns_ss1_1_PriLS$post_interval1924, x = indi_ns_ss1_1_PriLS$scoringD1_0,
                                 kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                 cluster = NULL)


### Outcome 1 for UppSec
# days worked post months [1, 6]
rddD1ei6_p1kT_UppSec <- rdrobust(y = indi_ns_ss1_1_UppSec$post_interval6, x = indi_ns_ss1_1_UppSec$scoringD1_0,
                                 kernel = "triangular",
                                 c = 0, p = 1, bwselect = "mserd",
                                 cluster = NULL)

### Outcome 2 for UppSec: days worked post months [7, 12]
rddD1ei12_p1kT_UppSec <- rdrobust(y = indi_ns_ss1_1_UppSec$post_interval712, x = indi_ns_ss1_1_UppSec$scoringD1_0,
                                  kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                  cluster = NULL)

### Outcome 3 for UppSec: days worked post [13, 18] months
rddD1ei18_p1kT_UppSec <- rdrobust(y = indi_ns_ss1_1_UppSec$post_interval1318, x = indi_ns_ss1_1_UppSec$scoringD1_0,
                                  kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                  cluster = NULL)


### Outcome 4 for UppSec: days worked post [19, 24] months
rddD1ei24_p1kT_UppSec <- rdrobust(y = indi_ns_ss1_1_UppSec$post_interval1924, x = indi_ns_ss1_1_UppSec$scoringD1_0,
                                  kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                  cluster = NULL)

### Outcome 1 for Uni
# days worked post months [1, 6]
rddD1ei6_p1kT_Uni <- rdrobust(y = indi_ns_ss1_1_Uni$post_interval6, x = indi_ns_ss1_1_Uni$scoringD1_0,
                              kernel = "triangular",
                              c = 0, p = 1, bwselect = "mserd",
                              cluster = NULL)

### Outcome 2 for Uni: days worked post months [7, 12]
rddD1ei12_p1kT_Uni <- rdrobust(y = indi_ns_ss1_1_Uni$post_interval712, x = indi_ns_ss1_1_Uni$scoringD1_0,
                               kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                               cluster = NULL)

### Outcome 3 for Uni: days worked post [13, 18] months
rddD1ei18_p1kT_Uni <- rdrobust(y = indi_ns_ss1_1_Uni$post_interval1318, x = indi_ns_ss1_1_Uni$scoringD1_0,
                               kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                               cluster = NULL)


### Outcome 4 for Uni: days worked post [19, 24] months
rddD1ei24_p1kT_Uni <- rdrobust(y = indi_ns_ss1_1_Uni$post_interval1924, x = indi_ns_ss1_1_Uni$scoringD1_0,
                               kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                               cluster = NULL)


#### Deeping: Impact on HS and HT ----------

# HS: hours of job search services
rddD1_HS_p1kT_UppSec <- rdrobust(y = indi_ns_ss1_1_UppSec$jshours, x = indi_ns_ss1_1_UppSec$scoringD1_0,
                                 kernel = "triangular",
                                 c = 0, p = 1, bwselect = "mserd",
                                 cluster = NULL)
summary(rddD1_HS_p1kT_UppSec)

# HT: hours of adult training
rddD1_HT_p1kT_UppSec <- rdrobust(y = indi_ns_ss1_1_UppSec$attiv_form_ore_prev, x = indi_ns_ss1_1_UppSec$scoringD1_0,
                                 kernel = "triangular",
                                 c = 0, p = 1, bwselect = "mserd",
                                 cluster = NULL)
summary(rddD1_HT_p1kT_UppSec)



### D2 ----------
indi_ns_ss2_1_PriLS <- subset(indi_ns_ss2_c,
                              studio2_grouped2 == "1_primary" | studio2_grouped2 == "2_lowersec")

indi_ns_ss2_1_UppSec <- subset(indi_ns_ss2_c,
                               studio2_grouped2 %in% c("344_uppersecG", "353_uppersecP_3y",
                                                       "353_uppersecP_4y_prof", "353_uppersecP_4y_tec"))

indi_ns_ss2_1_Uni <- subset(indi_ns_ss2_c,
                            studio2_grouped2 %in% c("660_terBachy+_lev1", "760+_terMastery+_lev2"))


## PriLS, UppSec, Uni


### Outcome 1 for PriLS: days worked post months [1, 6] --> Not significant effect (not significant difference) !!!!!
rddD2ei6_p1kT_PriLS <- rdrobust(y = indi_ns_ss2_1_PriLS$post_interval6, x = indi_ns_ss2_1_PriLS$scoringD2_0,
                                kernel = "triangular",
                                c = 0, p = 1, bwselect = "mserd",
                                cluster = NULL)

### Outcome 2 for PriLS: days worked post months [7, 12] --> Not significant effect (not ssignificant difference) !!!!!!!!!!!!!!!!
rddD2ei12_p1kT_PriLS <- rdrobust(y = indi_ns_ss2_1_PriLS$post_interval712, x = indi_ns_ss2_1_PriLS$scoringD2_0,
                                 kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                 cluster = NULL)

### Outcome 3 for PriLS: days worked post [13, 18] months
rddD2ei18_p1kT_PriLS <- rdrobust(y = indi_ns_ss2_1_PriLS$post_interval1318, x = indi_ns_ss2_1_PriLS$scoringD2_0,
                                 kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                 cluster = NULL)

### Outcome 4 for PriLS: days worked post [19, 24] months
rddD2ei24_p1kT_PriLS <- rdrobust(y = indi_ns_ss2_1_PriLS$post_interval1924, x = indi_ns_ss2_1_PriLS$scoringD2_0,
                                 kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                 cluster = NULL)


### Outcome 1 for UppSec
# days worked post months [1, 6] --> Not significant effect (not significant difference) !!!!!
rddD2ei6_p1kT_UppSec <- rdrobust(y = indi_ns_ss2_1_UppSec$post_interval6, x = indi_ns_ss2_1_UppSec$scoringD2_0,
                                 kernel = "triangular",
                                 c = 0, p = 1, bwselect = "mserd",
                                 cluster = NULL)

### Outcome 2 for UppSec: days worked post months [7, 12] --> Not significant effect (not ssignificant difference) !!!!!!!!!!!!!!!!
rddD2ei12_p1kT_UppSec <- rdrobust(y = indi_ns_ss2_1_UppSec$post_interval712, x = indi_ns_ss2_1_UppSec$scoringD2_0,
                                  kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                  cluster = NULL)

### Outcome 3 for UppSec: days worked post [13, 18] months
rddD2ei18_p1kT_UppSec <- rdrobust(y = indi_ns_ss2_1_UppSec$post_interval1318, x = indi_ns_ss2_1_UppSec$scoringD2_0,
                                  kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                  cluster = NULL)


### Outcome 4 for UppSec: days worked post [19, 24] months
rddD2ei24_p1kT_UppSec <- rdrobust(y = indi_ns_ss2_1_UppSec$post_interval1924, x = indi_ns_ss2_1_UppSec$scoringD2_0,
                                  kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                                  cluster = NULL)

### Outcome 1 for Uni
# days worked post months [1, 6] --> Not significant effect (not significant difference) !!!!!
rddD2ei6_p1kT_Uni <- rdrobust(y = indi_ns_ss2_1_Uni$post_interval6, x = indi_ns_ss2_1_Uni$scoringD2_0,
                              kernel = "triangular",
                              c = 0, p = 1, bwselect = "mserd",
                              cluster = NULL)

### Outcome 2 for Uni: days worked post months [7, 12] --> Not significant effect (not ssignificant difference) !!!!!!!!!!!!!!!!
rddD2ei12_p1kT_Uni <- rdrobust(y = indi_ns_ss2_1_Uni$post_interval712, x = indi_ns_ss2_1_Uni$scoringD2_0,
                               kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                               cluster = NULL)

### Outcome 3 for Uni: days worked post [13, 18] months
rddD2ei18_p1kT_Uni <- rdrobust(y = indi_ns_ss2_1_Uni$post_interval1318, x = indi_ns_ss2_1_Uni$scoringD2_0,
                               kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                               cluster = NULL)


### Outcome 4 for Uni: days worked post [19, 24] months
rddD2ei24_p1kT_Uni <- rdrobust(y = indi_ns_ss2_1_Uni$post_interval1924, x = indi_ns_ss2_1_Uni$scoringD2_0,
                               kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                               cluster = NULL)


## Depending on date ------------
### D1--------------

indi_ns_ss1_bCovid <- indi_ns_ss1_c %>% filter(ppa_data_avvio_d < ymd("2019-01-01"))  # started treatment before 01-01-2019 (at least 1 year of observation)

indi_ns_ss1_afCovid <- indi_ns_ss1_c %>% filter(ppa_data_avvio_d > ymd("2020-01-01")) # for Y1

## before vs after Covid


### Outcome 1 for before: days worked post months [1, 6]
rddD1ei6_p1kTb <- rdrobust(y = indi_ns_ss1_bCovid$post_interval6, x = indi_ns_ss1_bCovid$scoringD1_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd",
                           cluster = NULL)


### Outcome 2 for before: days worked post months [7, 12]
rddD1ei12_p1kTb <- rdrobust(y = indi_ns_ss1_bCovid$post_interval712, x = indi_ns_ss1_bCovid$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)


### Outcome 1 for after: days worked post months [1, 6]
rddD1ei6_p1kTa <- rdrobust(y = indi_ns_ss1_afCovid$post_interval6, x = indi_ns_ss1_afCovid$scoringD1_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd",
                           cluster = NULL)

### Outcome 2 for after: days worked post months [7, 12]
rddD1ei12_p1kTa <- rdrobust(y = indi_ns_ss1_afCovid$post_interval712, x = indi_ns_ss1_afCovid$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 3 for after: days worked post [13, 18] months
rddD1ei18_p1kTa <- rdrobust(y = indi_ns_ss1_afCovid$post_interval1318, x = indi_ns_ss1_afCovid$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 4 for after: days worked post [19, 24] months
rddD1ei24_p1kTa <- rdrobust(y = indi_ns_ss1_afCovid$post_interval1924, x = indi_ns_ss1_afCovid$scoringD1_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### D2 ----------

indi_ns_ss2_bCovid <- subset(indi_ns_ss2_c, ppa_data_avvio_d < ymd("2019-01-01"))  # started treatment before 01-01-2019 (at least 1 year of observation)

indi_ns_ss2_afCovid <- subset(indi_ns_ss2_c, ppa_data_avvio_d > ymd("2020-01-01"))  # for Y1

## before vs after Covid


### Outcome 1 for before: days worked post months [1, 6]
rddD2ei6_p1kTb <- rdrobust(y = indi_ns_ss2_bCovid$post_interval6, x = indi_ns_ss2_bCovid$scoringD2_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd",
                           cluster = NULL)


### Outcome 2 for before: days worked post months [7, 12]
rddD2ei12_p1kTb <- rdrobust(y = indi_ns_ss2_bCovid$post_interval712, x = indi_ns_ss2_bCovid$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)


### Outcome 1 for after: days worked post months [1, 6]
rddD2ei6_p1kTa <- rdrobust(y = indi_ns_ss2_afCovid$post_interval6, x = indi_ns_ss2_afCovid$scoringD2_0,
                           kernel = "triangular",
                           c = 0, p = 1, bwselect = "mserd",
                           cluster = NULL)

### Outcome 2 for after: days worked post months [7, 12]
rddD2ei12_p1kTa <- rdrobust(y = indi_ns_ss2_afCovid$post_interval712, x = indi_ns_ss2_afCovid$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 3 for after: days worked post [13, 18] months
rddD2ei18_p1kTa <- rdrobust(y = indi_ns_ss2_afCovid$post_interval1318, x = indi_ns_ss2_afCovid$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

### Outcome 4 for after: days worked post [19, 24] months
rddD2ei24_p1kTa <- rdrobust(y = indi_ns_ss2_afCovid$post_interval1924, x = indi_ns_ss2_afCovid$scoringD2_0,
                            kernel = "triangular", c = 0, p = 1, bwselect = "mserd",
                            cluster = NULL)

## Exporting -----------

### D1 ---------------

subgroupd1 <- data.frame(subgroup = c("Sex = Woman", "", "Sex = Man", "",
                                      "Age >= 50", "", "Age < 50", "",
                                      "Educ. = 1 or 2", "", "Educ. = 3", "", "Educ. = 6, 7 or 8", "",
                                      "Date < COVID-19", "", "Date >= COVID-19", ""),
                         Y1 = c(rddD1ei6_p1kTw$Estimate[1], rddD1ei6_p1kTw$Estimate[3],
                                rddD1ei6_p1kTm$Estimate[1], rddD1ei6_p1kTm$Estimate[3],
                                rddD1ei6_p1kTu$Estimate[1], rddD1ei6_p1kTu$Estimate[3],
                                rddD1ei6_p1kTd$Estimate[1], rddD1ei6_p1kTd$Estimate[3],
                                rddD1ei6_p1kT_PriLS$Estimate[1], rddD1ei6_p1kT_PriLS$Estimate[3],
                                rddD1ei6_p1kT_UppSec$Estimate[1], rddD1ei6_p1kT_UppSec$Estimate[3],
                                rddD1ei6_p1kT_Uni$Estimate[1], rddD1ei6_p1kT_Uni$Estimate[3],
                                rddD1ei6_p1kTb$Estimate[1], rddD1ei6_p1kTb$Estimate[3],
                                rddD1ei6_p1kTa$Estimate[1], rddD1ei6_p1kTa$Estimate[3]),
                         Y2 = c(rddD1ei12_p1kTw$Estimate[1], rddD1ei12_p1kTw$Estimate[3],
                                rddD1ei12_p1kTm$Estimate[1], rddD1ei12_p1kTm$Estimate[3],
                                rddD1ei12_p1kTu$Estimate[1], rddD1ei12_p1kTu$Estimate[3],
                                rddD1ei12_p1kTd$Estimate[1], rddD1ei12_p1kTd$Estimate[3],
                                rddD1ei12_p1kT_PriLS$Estimate[1], rddD1ei12_p1kT_PriLS$Estimate[3],
                                rddD1ei12_p1kT_UppSec$Estimate[1], rddD1ei12_p1kT_UppSec$Estimate[3],
                                rddD1ei12_p1kT_Uni$Estimate[1], rddD1ei12_p1kT_Uni$Estimate[3],
                                rddD1ei12_p1kTb$Estimate[1], rddD1ei12_p1kTb$Estimate[3],
                                rddD1ei12_p1kTa$Estimate[1], rddD1ei12_p1kTa$Estimate[3]),
                         Y3 = c(rddD1ei18_p1kTw$Estimate[1], rddD1ei18_p1kTw$Estimate[3],
                                rddD1ei18_p1kTm$Estimate[1], rddD1ei18_p1kTm$Estimate[3],
                                rddD1ei18_p1kTu$Estimate[1], rddD1ei18_p1kTu$Estimate[3],
                                rddD1ei18_p1kTd$Estimate[1], rddD1ei18_p1kTd$Estimate[3],
                                rddD1ei18_p1kT_PriLS$Estimate[1], rddD1ei18_p1kT_PriLS$Estimate[3],
                                rddD1ei18_p1kT_UppSec$Estimate[1], rddD1ei18_p1kT_UppSec$Estimate[3],
                                rddD1ei18_p1kT_Uni$Estimate[1], rddD1ei18_p1kT_Uni$Estimate[3],
                                NA, NA,
                                rddD1ei18_p1kTa$Estimate[1], rddD1ei18_p1kTa$Estimate[3]),
                         Y4 = c(rddD1ei24_p1kTw$Estimate[1], rddD1ei24_p1kTw$Estimate[3],
                                rddD1ei24_p1kTm$Estimate[1], rddD1ei24_p1kTm$Estimate[3],
                                rddD1ei24_p1kTu$Estimate[1], rddD1ei24_p1kTu$Estimate[3],
                                rddD1ei24_p1kTd$Estimate[1], rddD1ei24_p1kTd$Estimate[3],
                                rddD1ei24_p1kT_PriLS$Estimate[1], rddD1ei24_p1kT_PriLS$Estimate[3],
                                rddD1ei24_p1kT_UppSec$Estimate[1], rddD1ei24_p1kT_UppSec$Estimate[3],
                                rddD1ei24_p1kT_Uni$Estimate[1], rddD1ei24_p1kT_Uni$Estimate[3],
                                NA, NA,
                                rddD1ei24_p1kTa$Estimate[1], rddD1ei24_p1kTa$Estimate[3])
)

subgroupd1$Y1 <- round(subgroupd1$Y1, 3)
subgroupd1$Y2 <- round(subgroupd1$Y2, 3)
subgroupd1$Y3 <- round(subgroupd1$Y3, 3)
subgroupd1$Y4 <- round(subgroupd1$Y4, 3)


rddD1ei6_p1kTw <- starizer(rddD1ei6_p1kTw)
rddD1ei6_p1kTm <- starizer(rddD1ei6_p1kTm)
rddD1ei6_p1kTu <- starizer(rddD1ei6_p1kTu)
rddD1ei6_p1kTd <- starizer(rddD1ei6_p1kTd)
rddD1ei6_p1kT_PriLS <- starizer(rddD1ei6_p1kT_PriLS)
rddD1ei6_p1kT_UppSec <- starizer(rddD1ei6_p1kT_UppSec)
rddD1ei6_p1kT_Uni <- starizer(rddD1ei6_p1kT_Uni)
rddD1ei6_p1kTb <- starizer(rddD1ei6_p1kTb)
rddD1ei6_p1kTa <- starizer(rddD1ei6_p1kTa)

rddD1ei12_p1kTw <- starizer(rddD1ei12_p1kTw)
rddD1ei12_p1kTm <- starizer(rddD1ei12_p1kTm)
rddD1ei12_p1kTu <- starizer(rddD1ei12_p1kTu)
rddD1ei12_p1kTd <- starizer(rddD1ei12_p1kTd)
rddD1ei12_p1kT_PriLS <- starizer(rddD1ei12_p1kT_PriLS)
rddD1ei12_p1kT_UppSec <- starizer(rddD1ei12_p1kT_UppSec)
rddD1ei12_p1kT_Uni <- starizer(rddD1ei12_p1kT_Uni)
rddD1ei12_p1kTb <- starizer(rddD1ei12_p1kTb)
rddD1ei12_p1kTa <- starizer(rddD1ei12_p1kTa)

rddD1ei18_p1kTw <- starizer(rddD1ei18_p1kTw)
rddD1ei18_p1kTm <- starizer(rddD1ei18_p1kTm)
rddD1ei18_p1kTu <- starizer(rddD1ei18_p1kTu)
rddD1ei18_p1kTd <- starizer(rddD1ei18_p1kTd)
rddD1ei18_p1kT_PriLS <- starizer(rddD1ei18_p1kT_PriLS)
rddD1ei18_p1kT_UppSec <- starizer(rddD1ei18_p1kT_UppSec)
rddD1ei18_p1kT_Uni <- starizer(rddD1ei18_p1kT_Uni)
rddD1ei18_p1kTa <- starizer(rddD1ei18_p1kTa)

rddD1ei24_p1kTw <- starizer(rddD1ei24_p1kTw)
rddD1ei24_p1kTm <- starizer(rddD1ei24_p1kTm)
rddD1ei24_p1kTu <- starizer(rddD1ei24_p1kTu)
rddD1ei24_p1kTd <- starizer(rddD1ei24_p1kTd)
rddD1ei24_p1kT_PriLS <- starizer(rddD1ei24_p1kT_PriLS)
rddD1ei24_p1kT_UppSec <- starizer(rddD1ei24_p1kT_UppSec)
rddD1ei24_p1kT_Uni <- starizer(rddD1ei24_p1kT_Uni)
rddD1ei24_p1kTa <- starizer(rddD1ei24_p1kTa)


# Lista de nombres de objetos por cada fila impar (en orden)
objects_Y1 <- c("rddD1ei6_p1kTw", "rddD1ei6_p1kTm", "rddD1ei6_p1kTu", "rddD1ei6_p1kTd",
                "rddD1ei6_p1kT_PriLS", "rddD1ei6_p1kT_UppSec", "rddD1ei6_p1kT_Uni",
                "rddD1ei6_p1kTb", "rddD1ei6_p1kTa")
objects_Y2 <- gsub("ei6", "ei12", objects_Y1)
objects_Y3 <- gsub("ei6", "ei18", objects_Y1)
objects_Y3 <- objects_Y3[objects_Y3 != "rddD1ei18_p1kTb"]

objects_Y4 <- gsub("ei6", "ei24", objects_Y1)
objects_Y4 <- objects_Y4[objects_Y4 != "rddD1ei24_p1kTb"]


# Filas impares, donde se encuentran las estimaciones puntuales
row_indices <- seq(1, 17, by = 2)

# Adding stars
for (i in seq_along(row_indices)) {
  row <- row_indices[i]

  # Solo pegamos stars si el valor en subgroupd1 no es NA
  if (!is.na(subgroupd1$Y1[row])) {
    stars_Y1 <- get(objects_Y1[i])$stars
    subgroupd1$Y1[row] <- paste0(subgroupd1$Y1[row], stars_Y1)
  }

  if (!is.na(subgroupd1$Y2[row])) {
    stars_Y2 <- get(objects_Y2[i])$stars
    subgroupd1$Y2[row] <- paste0(subgroupd1$Y2[row], stars_Y2)
  }
  if (!is.na(subgroupd1$Y3[row]) & i <= length(objects_Y3)) {
    stars_Y3 <- get(objects_Y3[i])$stars
    subgroupd1$Y3[row] <- paste0(subgroupd1$Y3[row], stars_Y3)
  }

  if (!is.na(subgroupd1$Y4[row]) & i <= length(objects_Y4)) {
    stars_Y4 <- get(objects_Y4[i])$stars
    subgroupd1$Y4[row] <- paste0(subgroupd1$Y4[row], stars_Y4)
  }
}

# Adding parentheses
row_pares <- seq(2, nrow(subgroupd1), by = 2)

# Cols to modify
cols <- c("Y1", "Y2", "Y3", "Y4")

# Añadir paréntesis a los valores no NA
for (col in cols) {
  for (row in row_pares) {
    if (!is.na(subgroupd1[[col]][row])) {
      subgroupd1[[col]][row] <- paste0("(", subgroupd1[[col]][row], ")")
    }
  }
}


### D2 ------------------

subgroupD2 <- data.frame(subgroup = c("Sex = Woman", "", "Sex = Man", "",
                                      "Age >= 50", "", "Age < 50", "",
                                      "Educ. = 1 or 2", "", "Educ. = 3", "", "Educ. = 6, 7 or 8", "",
                                      "Date < COVID-19", "", "Date >= COVID-19", ""),
                         Y1 = c(rddD2ei6_p1kTw$Estimate[1], rddD2ei6_p1kTw$Estimate[3],
                                rddD2ei6_p1kTm$Estimate[1], rddD2ei6_p1kTm$Estimate[3],
                                rddD2ei6_p1kTu$Estimate[1], rddD2ei6_p1kTu$Estimate[3],
                                rddD2ei6_p1kTd$Estimate[1], rddD2ei6_p1kTd$Estimate[3],
                                rddD2ei6_p1kT_PriLS$Estimate[1], rddD2ei6_p1kT_PriLS$Estimate[3],
                                rddD2ei6_p1kT_UppSec$Estimate[1], rddD2ei6_p1kT_UppSec$Estimate[3],
                                rddD2ei6_p1kT_Uni$Estimate[1], rddD2ei6_p1kT_Uni$Estimate[3],
                                rddD2ei6_p1kTb$Estimate[1], rddD2ei6_p1kTb$Estimate[3],
                                rddD2ei6_p1kTa$Estimate[1], rddD2ei6_p1kTa$Estimate[3]),
                         Y2 = c(rddD2ei12_p1kTw$Estimate[1], rddD2ei12_p1kTw$Estimate[3],
                                rddD2ei12_p1kTm$Estimate[1], rddD2ei12_p1kTm$Estimate[3],
                                rddD2ei12_p1kTu$Estimate[1], rddD2ei12_p1kTu$Estimate[3],
                                rddD2ei12_p1kTd$Estimate[1], rddD2ei12_p1kTd$Estimate[3],
                                rddD2ei12_p1kT_PriLS$Estimate[1], rddD2ei12_p1kT_PriLS$Estimate[3],
                                rddD2ei12_p1kT_UppSec$Estimate[1], rddD2ei12_p1kT_UppSec$Estimate[3],
                                rddD2ei12_p1kT_Uni$Estimate[1], rddD2ei12_p1kT_Uni$Estimate[3],
                                rddD2ei12_p1kTb$Estimate[1], rddD2ei12_p1kTb$Estimate[3],
                                rddD2ei12_p1kTa$Estimate[1], rddD2ei12_p1kTa$Estimate[3]),
                         Y3 = c(rddD2ei18_p1kTw$Estimate[1], rddD2ei18_p1kTw$Estimate[3],
                                rddD2ei18_p1kTm$Estimate[1], rddD2ei18_p1kTm$Estimate[3],
                                rddD2ei18_p1kTu$Estimate[1], rddD2ei18_p1kTu$Estimate[3],
                                rddD2ei18_p1kTd$Estimate[1], rddD2ei18_p1kTd$Estimate[3],
                                rddD2ei18_p1kT_PriLS$Estimate[1], rddD2ei18_p1kT_PriLS$Estimate[3],
                                rddD2ei18_p1kT_UppSec$Estimate[1], rddD2ei18_p1kT_UppSec$Estimate[3],
                                rddD2ei18_p1kT_Uni$Estimate[1], rddD2ei18_p1kT_Uni$Estimate[3],
                                NA, NA,
                                rddD2ei18_p1kTa$Estimate[1], rddD2ei18_p1kTa$Estimate[3]),
                         Y4 = c(rddD2ei24_p1kTw$Estimate[1], rddD2ei24_p1kTw$Estimate[3],
                                rddD2ei24_p1kTm$Estimate[1], rddD2ei24_p1kTm$Estimate[3],
                                rddD2ei24_p1kTu$Estimate[1], rddD2ei24_p1kTu$Estimate[3],
                                rddD2ei24_p1kTd$Estimate[1], rddD2ei24_p1kTd$Estimate[3],
                                rddD2ei24_p1kT_PriLS$Estimate[1], rddD2ei24_p1kT_PriLS$Estimate[3],
                                rddD2ei24_p1kT_UppSec$Estimate[1], rddD2ei24_p1kT_UppSec$Estimate[3],
                                rddD2ei24_p1kT_Uni$Estimate[1], rddD2ei24_p1kT_Uni$Estimate[3],
                                NA, NA,
                                rddD2ei24_p1kTa$Estimate[1], rddD2ei24_p1kTa$Estimate[3])
)

subgroupD2$Y1 <- round(subgroupD2$Y1, 3)
subgroupD2$Y2 <- round(subgroupD2$Y2, 3)
subgroupD2$Y3 <- round(subgroupD2$Y3, 3)
subgroupD2$Y4 <- round(subgroupD2$Y4, 3)


rddD2ei6_p1kTw <- starizer(rddD2ei6_p1kTw)
rddD2ei6_p1kTm <- starizer(rddD2ei6_p1kTm)
rddD2ei6_p1kTu <- starizer(rddD2ei6_p1kTu)
rddD2ei6_p1kTd <- starizer(rddD2ei6_p1kTd)
rddD2ei6_p1kT_PriLS <- starizer(rddD2ei6_p1kT_PriLS)
rddD2ei6_p1kT_UppSec <- starizer(rddD2ei6_p1kT_UppSec)
rddD2ei6_p1kT_Uni <- starizer(rddD2ei6_p1kT_Uni)
rddD2ei6_p1kTb <- starizer(rddD2ei6_p1kTb)
rddD2ei6_p1kTa <- starizer(rddD2ei6_p1kTa)

rddD2ei12_p1kTw <- starizer(rddD2ei12_p1kTw)
rddD2ei12_p1kTm <- starizer(rddD2ei12_p1kTm)
rddD2ei12_p1kTu <- starizer(rddD2ei12_p1kTu)
rddD2ei12_p1kTd <- starizer(rddD2ei12_p1kTd)
rddD2ei12_p1kT_PriLS <- starizer(rddD2ei12_p1kT_PriLS)
rddD2ei12_p1kT_UppSec <- starizer(rddD2ei12_p1kT_UppSec)
rddD2ei12_p1kT_Uni <- starizer(rddD2ei12_p1kT_Uni)
rddD2ei12_p1kTb <- starizer(rddD2ei12_p1kTb)
rddD2ei12_p1kTa <- starizer(rddD2ei12_p1kTa)

rddD2ei18_p1kTw <- starizer(rddD2ei18_p1kTw)
rddD2ei18_p1kTm <- starizer(rddD2ei18_p1kTm)
rddD2ei18_p1kTu <- starizer(rddD2ei18_p1kTu)
rddD2ei18_p1kTd <- starizer(rddD2ei18_p1kTd)
rddD2ei18_p1kT_PriLS <- starizer(rddD2ei18_p1kT_PriLS)
rddD2ei18_p1kT_UppSec <- starizer(rddD2ei18_p1kT_UppSec)
rddD2ei18_p1kT_Uni <- starizer(rddD2ei18_p1kT_Uni)
rddD2ei18_p1kTa <- starizer(rddD2ei18_p1kTa)

rddD2ei24_p1kTw <- starizer(rddD2ei24_p1kTw)
rddD2ei24_p1kTm <- starizer(rddD2ei24_p1kTm)
rddD2ei24_p1kTu <- starizer(rddD2ei24_p1kTu)
rddD2ei24_p1kTd <- starizer(rddD2ei24_p1kTd)
rddD2ei24_p1kT_PriLS <- starizer(rddD2ei24_p1kT_PriLS)
rddD2ei24_p1kT_UppSec <- starizer(rddD2ei24_p1kT_UppSec)
rddD2ei24_p1kT_Uni <- starizer(rddD2ei24_p1kT_Uni)
rddD2ei24_p1kTa <- starizer(rddD2ei24_p1kTa)


# Lista de nombres de objetos por cada fila impar (en orden)
objects_Y1 <- c("rddD2ei6_p1kTw", "rddD2ei6_p1kTm", "rddD2ei6_p1kTu", "rddD2ei6_p1kTd",
                "rddD2ei6_p1kT_PriLS", "rddD2ei6_p1kT_UppSec", "rddD2ei6_p1kT_Uni",
                "rddD2ei6_p1kTb", "rddD2ei6_p1kTa")
objects_Y2 <- gsub("ei6", "ei12", objects_Y1)
objects_Y3 <- gsub("ei6", "ei18", objects_Y1)
objects_Y3 <- objects_Y3[objects_Y3 != "rddD2ei18_p1kTb"]

objects_Y4 <- gsub("ei6", "ei24", objects_Y1)
objects_Y4 <- objects_Y4[objects_Y4 != "rddD2ei24_p1kTb"]


# Filas impares, donde se encuentran las estimaciones puntuales
row_indices <- seq(1, 17, by = 2)

# Adding stars
for (i in seq_along(row_indices)) {
  row <- row_indices[i]

  # Solo pegamos stars si el valor en subgroupD2 no es NA
  if (!is.na(subgroupD2$Y1[row])) {
    stars_Y1 <- get(objects_Y1[i])$stars
    subgroupD2$Y1[row] <- paste0(subgroupD2$Y1[row], stars_Y1)
  }

  if (!is.na(subgroupD2$Y2[row])) {
    stars_Y2 <- get(objects_Y2[i])$stars
    subgroupD2$Y2[row] <- paste0(subgroupD2$Y2[row], stars_Y2)
  }
  if (!is.na(subgroupD2$Y3[row]) & i <= length(objects_Y3)) {
    stars_Y3 <- get(objects_Y3[i])$stars
    subgroupD2$Y3[row] <- paste0(subgroupD2$Y3[row], stars_Y3)
  }

  if (!is.na(subgroupD2$Y4[row]) & i <= length(objects_Y4)) {
    stars_Y4 <- get(objects_Y4[i])$stars
    subgroupD2$Y4[row] <- paste0(subgroupD2$Y4[row], stars_Y4)
  }
}

# Adding parentheses
row_pares <- seq(2, nrow(subgroupD2), by = 2)

## Cols to modify
cols <- c("Y1", "Y2", "Y3", "Y4")

## Añadir paréntesis a los valores no NA
for (col in cols) {
  for (row in row_pares) {
    if (!is.na(subgroupD2[[col]][row])) {
      subgroupD2[[col]][row] <- paste0("(", subgroupD2[[col]][row], ")")
    }
  }
}

# 3. LATE_j heterogeneity (subgroup analyses) [multinomial and binomial logit models] ----------------

## Depending on sex ------------
### D1--------------

rddD1_type_p1kT_F <- rdcate_multinom(score = indi_ns_ss1_F$scoringD1_0,
                                     outcome = indi_ns_ss1_F$longestcontract_standardo)

rddD1_type_p1kT_M <- rdcate_multinom(score = indi_ns_ss1_M$scoringD1_0,
                                     outcome = indi_ns_ss1_M$longestcontract_standardo)

### Outcome lcris_OEC
rddD1_OEC_p1kT_F <- rdcate(outcome = indi_ns_ss1_F$lcris_OEC,
                           score = indi_ns_ss1_F$scoringD1_0)

### Outcome lcris_FTCm12
rddD1_FTCm12_p1kT_F <- rdcate(outcome = indi_ns_ss1_F$lcris_FTCm12,
                              score = indi_ns_ss1_F$scoringD1_0)

### Outcome lcris_FTC612
rddD1_FTC612_p1kT_F <- rdcate(outcome = indi_ns_ss1_F$lcris_FTC612,
                              score = indi_ns_ss1_F$scoringD1_0)


### Outcome lcris_FTC612
rddD1_FTC05_p1kT_F <- rdcate(outcome = indi_ns_ss1_F$lcris_FTCme6,
                             score = indi_ns_ss1_F$scoringD1_0)


### Outcome lcris_OEC
rddD1_OEC_p1kT_M <- rdcate(outcome = indi_ns_ss1_M$lcris_OEC,
                           score = indi_ns_ss1_M$scoringD1_0)

### Outcome lcris_FTCm12
rddD1_FTCm12_p1kT_M <- rdcate(outcome = indi_ns_ss1_M$lcris_FTCm12,
                              score = indi_ns_ss1_M$scoringD1_0)

### Outcome lcris_FTC612
rddD1_FTC612_p1kT_M <- rdcate(outcome = indi_ns_ss1_M$lcris_FTC612,
                              score = indi_ns_ss1_M$scoringD1_0)

### Outcome lcris_FTC05
rddD1_FTC05_p1kT_M <- rdcate(outcome = indi_ns_ss1_M$lcris_FTCme6,
                             score = indi_ns_ss1_M$scoringD1_0)


quality_sex_D1_hetero <- data.frame(
  treatment = c(rep("D1", 2)),
  prOEC_F = c(
    paste0(rddD1_OEC_p1kT_F$pointest, rddD1_OEC_p1kT_F$stars),
    rddD1_OEC_p1kT_F$ci95
  ),
  prFTCm12_F = c(
    paste0(rddD1_FTCm12_p1kT_F$pointest, rddD1_FTCm12_p1kT_F$stars),
    rddD1_FTCm12_p1kT_F$ci95
  ),
  prFTC612_F = c(
    paste0(rddD1_FTC612_p1kT_F$pointest, rddD1_FTC612_p1kT_F$stars),
    rddD1_FTC612_p1kT_F$ci95
  ),
  prFTC05_F = c(
    paste0(rddD1_FTC05_p1kT_F$pointest, rddD1_FTC05_p1kT_F$stars),
    rddD1_FTC05_p1kT_F$ci95
  ),
  prOEC_M = c(
    paste0(rddD1_OEC_p1kT_M$pointest, rddD1_OEC_p1kT_M$stars),
    rddD1_OEC_p1kT_M$ci95
  ),
  prFTCm12_M = c(
    paste0(rddD1_FTCm12_p1kT_M$pointest, rddD1_FTCm12_p1kT_M$stars),
    rddD1_FTCm12_p1kT_M$ci95
  ),
  prFTC612_M = c(
    paste0(rddD1_FTC612_p1kT_M$pointest, rddD1_FTC612_p1kT_M$stars),
    rddD1_FTC612_p1kT_M$ci95
  ),
  prFTC05_M = c(
    paste0(rddD1_FTC05_p1kT_M$pointest, rddD1_FTC05_p1kT_M$stars),
    rddD1_FTC05_p1kT_M$ci95
  )
)


#write.xlsx(quality_sex_D1_hetero, 'quality_sex_D1_hetero.xlsx')

### D2 ------------------

rddD2_type_p1kT_F <- rdcate_multinom(score = indi_ns_ss2_F$scoringD2_0,
                                     outcome = indi_ns_ss2_F$longestcontract_standardo)

rddD2_type_p1kT_M <- rdcate_multinom(score = indi_ns_ss2_M$scoringD2_0,
                                     outcome = indi_ns_ss2_M$longestcontract_standardo)

### Outcome lcris_OEC
rddD2_OEC_p1kT_F <- rdcate(outcome = indi_ns_ss2_F$lcris_OEC,
                           score = indi_ns_ss2_F$scoringD2_0)

### Outcome lcris_FTCm12
rddD2_FTCm12_p1kT_F <- rdcate(outcome = indi_ns_ss2_F$lcris_FTCm12,
                              score = indi_ns_ss2_F$scoringD2_0)

### Outcome lcris_FTC612
rddD2_FTC612_p1kT_F <- rdcate(outcome = indi_ns_ss2_F$lcris_FTC612,
                              score = indi_ns_ss2_F$scoringD2_0)


### Outcome lcris_FTC612
rddD2_FTC05_p1kT_F <- rdcate(outcome = indi_ns_ss2_F$lcris_FTCme6,
                             score = indi_ns_ss2_F$scoringD2_0)


### Outcome lcris_OEC
rddD2_OEC_p1kT_M <- rdcate(outcome = indi_ns_ss2_M$lcris_OEC,
                           score = indi_ns_ss2_M$scoringD2_0)

### Outcome lcris_FTCm12
rddD2_FTCm12_p1kT_M <- rdcate(outcome = indi_ns_ss2_M$lcris_FTCm12,
                              score = indi_ns_ss2_M$scoringD2_0)

### Outcome lcris_FTC612
rddD2_FTC612_p1kT_M <- rdcate(outcome = indi_ns_ss2_M$lcris_FTC612,
                              score = indi_ns_ss2_M$scoringD2_0)

### Outcome lcris_FTC05
rddD2_FTC05_p1kT_M <- rdcate(outcome = indi_ns_ss2_M$lcris_FTCme6,
                             score = indi_ns_ss2_M$scoringD2_0)


quality_sex_D2_hetero <- data.frame(
  treatment = c(rep("D2", 2)),
  prOEC_F = c(
    paste0(rddD2_OEC_p1kT_F$pointest, rddD2_OEC_p1kT_F$stars),
    rddD2_OEC_p1kT_F$ci95
  ),
  prFTCm12_F = c(
    paste0(rddD2_FTCm12_p1kT_F$pointest, rddD2_FTCm12_p1kT_F$stars),
    rddD2_FTCm12_p1kT_F$ci95
  ),
  prFTC612_F = c(
    paste0(rddD2_FTC612_p1kT_F$pointest, rddD2_FTC612_p1kT_F$stars),
    rddD2_FTC612_p1kT_F$ci95
  ),
  prFTC05_F = c(
    paste0(rddD2_FTC05_p1kT_F$pointest, rddD2_FTC05_p1kT_F$stars),
    rddD2_FTC05_p1kT_F$ci95
  ),
  prOEC_M = c(
    paste0(rddD2_OEC_p1kT_M$pointest, rddD2_OEC_p1kT_M$stars),
    rddD2_OEC_p1kT_M$ci95
  ),
  prFTCm12_M = c(
    paste0(rddD2_FTCm12_p1kT_M$pointest, rddD2_FTCm12_p1kT_M$stars),
    rddD2_FTCm12_p1kT_M$ci95
  ),
  prFTC612_M = c(
    paste0(rddD2_FTC612_p1kT_M$pointest, rddD2_FTC612_p1kT_M$stars),
    rddD2_FTC612_p1kT_M$ci95
  ),
  prFTC05_M = c(
    paste0(rddD2_FTC05_p1kT_M$pointest, rddD2_FTC05_p1kT_M$stars),
    rddD2_FTC05_p1kT_M$ci95
  )
)


#write.xlsx(quality_sex_D2_hetero, 'quality_sex_D2_hetero.xlsx')

## Depending on age ------------
### D1--------------


rddD1_type_p1kT_ageup <- rdcate_multinom(score = indi_ns_ss1_ageup$scoringD1_0,
                                         outcome = indi_ns_ss1_ageup$longestcontract_standardo)

rddD1_type_p1kT_agedown <- rdcate_multinom(score = indi_ns_ss1_agedown$scoringD1_0,
                                           outcome = indi_ns_ss1_agedown$longestcontract_standardo)

rddD1_type_p1kT_ageup
rddD1_type_p1kT_agedown

### Outcome lcris_OEC
rddD1_OEC_p1kT_ageup <- rdcate(outcome = indi_ns_ss1_ageup$lcris_OEC,
                               score = indi_ns_ss1_ageup$scoringD1_0)

### Outcome lcris_FTCm12
rddD1_ageupTCm12_p1kT_ageup <- rdcate(outcome = indi_ns_ss1_ageup$lcris_FTCm12,
                                      score = indi_ns_ss1_ageup$scoringD1_0)

### Outcome lcris_FTC612
rddD1_ageupTC612_p1kT_ageup <- rdcate(outcome = indi_ns_ss1_ageup$lcris_FTC612,
                                      score = indi_ns_ss1_ageup$scoringD1_0)

### Outcome lcris_FTC05
rddD1_ageupTC05_p1kT_ageup <- rdcate(outcome = indi_ns_ss1_ageup$lcris_FTCme6,
                                     score = indi_ns_ss1_ageup$scoringD1_0)


### Outcome lcris_OEC
rddD1_OEC_p1kT_agedown <- rdcate(outcome = indi_ns_ss1_agedown$lcris_OEC,
                                 score = indi_ns_ss1_agedown$scoringD1_0)

### Outcome lcris_FTCm12
rddD1_ageupTCm12_p1kT_agedown <- rdcate(outcome = indi_ns_ss1_agedown$lcris_FTCm12,
                                        score = indi_ns_ss1_agedown$scoringD1_0)

### Outcome lcris_FTC612
rddD1_ageupTC612_p1kT_agedown <- rdcate(outcome = indi_ns_ss1_agedown$lcris_FTC612,
                                        score = indi_ns_ss1_agedown$scoringD1_0)

### Outcome lcris_FTC05
rddD1_ageupTC05_p1kT_agedown <- rdcate(outcome = indi_ns_ss1_agedown$lcris_FTCme6,
                                       score = indi_ns_ss1_agedown$scoringD1_0)



quality_age_D1_hetero <- data.frame(
  treatment = c(rep("D1", 2)),
  prOEC_ageup = c(
    paste0(rddD1_OEC_p1kT_ageup$pointest, rddD1_OEC_p1kT_ageup$stars),
    rddD1_OEC_p1kT_ageup$ci95
  ),
  prFTCm12_ageup = c(
    paste0(rddD1_ageupTCm12_p1kT_ageup$pointest, rddD1_ageupTCm12_p1kT_ageup$stars),
    rddD1_ageupTCm12_p1kT_ageup$ci95
  ),
  prFTC612_ageup = c(
    paste0(rddD1_ageupTC612_p1kT_ageup$pointest, rddD1_ageupTC612_p1kT_ageup$stars),
    rddD1_ageupTC612_p1kT_ageup$ci95
  ),
  prFTC05_ageup = c(
    paste0(rddD1_ageupTC05_p1kT_ageup$pointest, rddD1_ageupTC05_p1kT_ageup$stars),
    rddD1_ageupTC05_p1kT_ageup$ci95
  ),
  prOEC_agedown = c(
    paste0(rddD1_OEC_p1kT_agedown$pointest, rddD1_OEC_p1kT_agedown$stars),
    rddD1_OEC_p1kT_agedown$ci95
  ),
  prFTCm12_agedown = c(
    paste0(rddD1_ageupTCm12_p1kT_agedown$pointest, rddD1_ageupTCm12_p1kT_agedown$stars),
    rddD1_ageupTCm12_p1kT_agedown$ci95
  ),
  prFTC612_agedown = c(
    paste0(rddD1_ageupTC612_p1kT_agedown$pointest, rddD1_ageupTC612_p1kT_agedown$stars),
    rddD1_ageupTC612_p1kT_agedown$ci95
  ),
  prFTC05_agedown = c(
    paste0(rddD1_ageupTC05_p1kT_agedown$pointest, rddD1_ageupTC05_p1kT_agedown$stars),
    rddD1_ageupTC05_p1kT_agedown$ci95
  )
)


write.xlsx(quality_age_D1_hetero, 'quality_age_D1_hetero.xlsx')


### D2--------------


rddD2_type_p1kT_ageup <- rdcate_multinom(score = indi_ns_ss2_ageup$scoringD2_0,
                                         outcome = indi_ns_ss2_ageup$longestcontract_standardo)

rddD2_type_p1kT_agedown <- rdcate_multinom(score = indi_ns_ss2_agedown$scoringD2_0,
                                           outcome = indi_ns_ss2_agedown$longestcontract_standardo)

rddD2_type_p1kT_ageup
rddD2_type_p1kT_agedown

### Outcome lcris_OEC
rddD2_OEC_p1kT_ageup <- rdcate(outcome = indi_ns_ss2_ageup$lcris_OEC,
                               score = indi_ns_ss2_ageup$scoringD2_0)

### Outcome lcris_FTCm12
rddD2_ageupTCm12_p1kT_ageup <- rdcate(outcome = indi_ns_ss2_ageup$lcris_FTCm12,
                                      score = indi_ns_ss2_ageup$scoringD2_0)

### Outcome lcris_FTC612
rddD2_ageupTC612_p1kT_ageup <- rdcate(outcome = indi_ns_ss2_ageup$lcris_FTC612,
                                      score = indi_ns_ss2_ageup$scoringD2_0)

### Outcome lcris_FTC05
rddD2_ageupTC05_p1kT_ageup <- rdcate(outcome = indi_ns_ss2_ageup$lcris_FTCme6,
                                     score = indi_ns_ss2_ageup$scoringD2_0)


### Outcome lcris_OEC
rddD2_OEC_p1kT_agedown <- rdcate(outcome = indi_ns_ss2_agedown$lcris_OEC,
                                 score = indi_ns_ss2_agedown$scoringD2_0)

### Outcome lcris_FTCm12
rddD2_ageupTCm12_p1kT_agedown <- rdcate(outcome = indi_ns_ss2_agedown$lcris_FTCm12,
                                        score = indi_ns_ss2_agedown$scoringD2_0)

### Outcome lcris_FTC612
rddD2_ageupTC612_p1kT_agedown <- rdcate(outcome = indi_ns_ss2_agedown$lcris_FTC612,
                                        score = indi_ns_ss2_agedown$scoringD2_0)

### Outcome lcris_FTC05
rddD2_ageupTC05_p1kT_agedown <- rdcate(outcome = indi_ns_ss2_agedown$lcris_FTCme6,
                                       score = indi_ns_ss2_agedown$scoringD2_0)



quality_age_D2_hetero <- data.frame(
  treatment = c(rep("D2", 2)),
  prOEC_ageup = c(
    paste0(rddD2_OEC_p1kT_ageup$pointest, rddD2_OEC_p1kT_ageup$stars),
    rddD2_OEC_p1kT_ageup$ci95
  ),
  prFTCm12_ageup = c(
    paste0(rddD2_ageupTCm12_p1kT_ageup$pointest, rddD2_ageupTCm12_p1kT_ageup$stars),
    rddD2_ageupTCm12_p1kT_ageup$ci95
  ),
  prFTC612_ageup = c(
    paste0(rddD2_ageupTC612_p1kT_ageup$pointest, rddD2_ageupTC612_p1kT_ageup$stars),
    rddD2_ageupTC612_p1kT_ageup$ci95
  ),
  prFTC05_ageup = c(
    paste0(rddD2_ageupTC05_p1kT_ageup$pointest, rddD2_ageupTC05_p1kT_ageup$stars),
    rddD2_ageupTC05_p1kT_ageup$ci95
  ),
  prOEC_agedown = c(
    paste0(rddD2_OEC_p1kT_agedown$pointest, rddD2_OEC_p1kT_agedown$stars),
    rddD2_OEC_p1kT_agedown$ci95
  ),
  prFTCm12_agedown = c(
    paste0(rddD2_ageupTCm12_p1kT_agedown$pointest, rddD2_ageupTCm12_p1kT_agedown$stars),
    rddD2_ageupTCm12_p1kT_agedown$ci95
  ),
  prFTC612_agedown = c(
    paste0(rddD2_ageupTC612_p1kT_agedown$pointest, rddD2_ageupTC612_p1kT_agedown$stars),
    rddD2_ageupTC612_p1kT_agedown$ci95
  ),
  prFTC05_agedown = c(
    paste0(rddD2_ageupTC05_p1kT_agedown$pointest, rddD2_ageupTC05_p1kT_agedown$stars),
    rddD2_ageupTC05_p1kT_agedown$ci95
  )
)


write.xlsx(quality_age_D2_hetero, 'quality_age_D2_hetero.xlsx')



## Depending on education ------------
### D1--------------

rddD1_type_p1kT_PriLS <- rdcate_multinom(score = indi_ns_ss1_1_PriLS$scoringD1_0,
                                         outcome = indi_ns_ss1_1_PriLS$longestcontract_standardo)

rddD1_type_p1kT_UppSec <- rdcate_multinom(score = indi_ns_ss1_1_UppSec$scoringD1_0,
                                          outcome = indi_ns_ss1_1_UppSec$longestcontract_standardo)

rddD1_type_p1kT_Uni <- rdcate_multinom(score = indi_ns_ss1_1_Uni$scoringD1_0,
                                       outcome = indi_ns_ss1_1_Uni$longestcontract_standardo)

rddD1_type_p1kT_PriLS
rddD1_type_p1kT_UppSec
rddD1_type_p1kT_Uni


### Outcome lcris_OEC
rddD1_OEC_p1kT_PriLS <- rdcate(outcome = indi_ns_ss1_1_PriLS$lcris_OEC,
                               score = indi_ns_ss1_1_PriLS$scoringD1_0)

### Outcome lcris_FTCm12
rddD1_FTCm12_p1kT_PriLS <- rdcate(outcome = indi_ns_ss1_1_PriLS$lcris_FTCm12,
                                  score = indi_ns_ss1_1_PriLS$scoringD1_0)

### Outcome lcris_FTC612
rddD1_FTC612_p1kT_PriLS <- rdcate(outcome = indi_ns_ss1_1_PriLS$lcris_FTC612,
                                  score = indi_ns_ss1_1_PriLS$scoringD1_0)

### Outcome lcris_FTC05
rddD1_FTC05_p1kT_PriLS <- rdcate(outcome = indi_ns_ss1_1_PriLS$lcris_FTCme6,
                                 score = indi_ns_ss1_1_PriLS$scoringD1_0)

### Outcome lcris_OEC UppSec
rddD1_OEC_p1kT_UppSec <- rdcate(outcome = indi_ns_ss1_1_UppSec$lcris_OEC,
                                score = indi_ns_ss1_1_UppSec$scoringD1_0)

### Outcome lcris_FTCm12 UppSec
rddD1_FTCm12_p1kT_UppSec <- rdcate(outcome = indi_ns_ss1_1_UppSec$lcris_FTCm12,
                                   score = indi_ns_ss1_1_UppSec$scoringD1_0)

### Outcome lcris_FTC612 UppSec
rddD1_FTC612_p1kT_UppSec <- rdcate(outcome = indi_ns_ss1_1_UppSec$lcris_FTC612,
                                   score = indi_ns_ss1_1_UppSec$scoringD1_0)

### Outcome lcris_FTC05 UppSec
rddD1_FTC05_p1kT_UppSec <- rdcate(outcome = indi_ns_ss1_1_UppSec$lcris_FTCme6,
                                  score = indi_ns_ss1_1_UppSec$scoringD1_0)


### Outcome lcris_OEC Uni
rddD1_OEC_p1kT_Uni <- rdcate(outcome = indi_ns_ss1_1_Uni$lcris_OEC,
                             score = indi_ns_ss1_1_Uni$scoringD1_0)

### Outcome lcris_FTCm12 Uni
rddD1_FTCm12_p1kT_Uni <- rdcate(outcome = indi_ns_ss1_1_Uni$lcris_FTCm12,
                                score = indi_ns_ss1_1_Uni$scoringD1_0)

### Outcome lcris_FTC612 Uni
rddD1_FTC612_p1kT_Uni <- rdcate(outcome = indi_ns_ss1_1_Uni$lcris_FTC612,
                                score = indi_ns_ss1_1_Uni$scoringD1_0)

### Outcome lcris_FTC05 Uni
rddD1_FTC05_p1kT_Uni <- rdcate(outcome = indi_ns_ss1_1_Uni$lcris_FTCme6,
                               score = indi_ns_ss1_1_Uni$scoringD1_0)

quality_education_D1_hetero <- data.frame(
  treatment = c(rep("D1", 2)),
  prOEC_PriLS = c(
    paste0(rddD1_OEC_p1kT_PriLS$pointest, rddD1_OEC_p1kT_PriLS$stars),
    rddD1_OEC_p1kT_PriLS$ci95
  ),
  prFTCm12_PriLS = c(
    paste0(rddD1_FTCm12_p1kT_PriLS$pointest, rddD1_FTCm12_p1kT_PriLS$stars),
    rddD1_FTCm12_p1kT_PriLS$ci95
  ),
  prFTC612_PriLS = c(
    paste0(rddD1_FTC612_p1kT_PriLS$pointest, rddD1_FTC612_p1kT_PriLS$stars),
    rddD1_FTC612_p1kT_PriLS$ci95
  ),
  prFTC05_PriLS = c(
    paste0(rddD1_FTC05_p1kT_PriLS$pointest, rddD1_FTC05_p1kT_PriLS$stars),
    rddD1_FTC05_p1kT_PriLS$ci95
  ),
  prOEC_UppSec = c(
    paste0(rddD1_OEC_p1kT_UppSec$pointest, rddD1_OEC_p1kT_UppSec$stars),
    rddD1_OEC_p1kT_UppSec$ci95
  ),
  prFTCm12_UppSec = c(
    paste0(rddD1_FTCm12_p1kT_UppSec$pointest, rddD1_FTCm12_p1kT_UppSec$stars),
    rddD1_FTCm12_p1kT_UppSec$ci95
  ),
  prFTC612_UppSec = c(
    paste0(rddD1_FTC612_p1kT_UppSec$pointest, rddD1_FTC612_p1kT_UppSec$stars),
    rddD1_FTC612_p1kT_UppSec$ci95
  ),
  prFTC05_UppSec = c(
    paste0(rddD1_FTC05_p1kT_UppSec$pointest, rddD1_FTC05_p1kT_UppSec$stars),
    rddD1_FTC05_p1kT_UppSec$ci95
  ),
  prOEC_Uni = c(
    paste0(rddD1_OEC_p1kT_Uni$pointest, rddD1_OEC_p1kT_Uni$stars),
    rddD1_OEC_p1kT_Uni$ci95
  ),
  prFTCm12_Uni = c(
    paste0(rddD1_FTCm12_p1kT_Uni$pointest, rddD1_FTCm12_p1kT_Uni$stars),
    rddD1_FTCm12_p1kT_Uni$ci95
  ),
  prFTC612_Uni = c(
    paste0(rddD1_FTC612_p1kT_Uni$pointest, rddD1_FTC612_p1kT_Uni$stars),
    rddD1_FTC612_p1kT_Uni$ci95
  ),
  prFTC05_Uni = c(
    paste0(rddD1_FTC05_p1kT_Uni$pointest, rddD1_FTC05_p1kT_Uni$stars),
    rddD1_FTC05_p1kT_Uni$ci95
  )
)


write.xlsx(quality_education_D1_hetero, 'quality_education_D1_hetero.xlsx')

### D2--------------

rddD2_type_p1kT_PriLS <- rdcate_multinom(score = indi_ns_ss2_1_PriLS$scoringD2_0,
                                         outcome = indi_ns_ss2_1_PriLS$longestcontract_standardo)

rddD2_type_p1kT_UppSec <- rdcate_multinom(score = indi_ns_ss2_1_UppSec$scoringD2_0,
                                          outcome = indi_ns_ss2_1_UppSec$longestcontract_standardo)

rddD2_type_p1kT_Uni <- rdcate_multinom(score = indi_ns_ss2_1_Uni$scoringD2_0,
                                       outcome = indi_ns_ss2_1_Uni$longestcontract_standardo)

rddD2_type_p1kT_PriLS
rddD2_type_p1kT_UppSec
rddD2_type_p1kT_Uni


### Outcome lcris_OEC
rddD2_OEC_p1kT_PriLS <- rdcate(outcome = indi_ns_ss2_1_PriLS$lcris_OEC,
                               score = indi_ns_ss2_1_PriLS$scoringD2_0)

### Outcome lcris_FTCm12
rddD2_FTCm12_p1kT_PriLS <- rdcate(outcome = indi_ns_ss2_1_PriLS$lcris_FTCm12,
                                  score = indi_ns_ss2_1_PriLS$scoringD2_0)

### Outcome lcris_FTC612
rddD2_FTC612_p1kT_PriLS <- rdcate(outcome = indi_ns_ss2_1_PriLS$lcris_FTC612,
                                  score = indi_ns_ss2_1_PriLS$scoringD2_0)

### Outcome lcris_FTC05
rddD2_FTC05_p1kT_PriLS <- rdcate(outcome = indi_ns_ss2_1_PriLS$lcris_FTCme6,
                                 score = indi_ns_ss2_1_PriLS$scoringD2_0)

### Outcome lcris_OEC UppSec
rddD2_OEC_p1kT_UppSec <- rdcate(outcome = indi_ns_ss2_1_UppSec$lcris_OEC,
                                score = indi_ns_ss2_1_UppSec$scoringD2_0)

### Outcome lcris_FTCm12 UppSec
rddD2_FTCm12_p1kT_UppSec <- rdcate(outcome = indi_ns_ss2_1_UppSec$lcris_FTCm12,
                                   score = indi_ns_ss2_1_UppSec$scoringD2_0)

### Outcome lcris_FTC612 UppSec
rddD2_FTC612_p1kT_UppSec <- rdcate(outcome = indi_ns_ss2_1_UppSec$lcris_FTC612,
                                   score = indi_ns_ss2_1_UppSec$scoringD2_0)

### Outcome lcris_FTC05 UppSec
rddD2_FTC05_p1kT_UppSec <- rdcate(outcome = indi_ns_ss2_1_UppSec$lcris_FTCme6,
                                  score = indi_ns_ss2_1_UppSec$scoringD2_0)


### Outcome lcris_OEC Uni
rddD2_OEC_p1kT_Uni <- rdcate(outcome = indi_ns_ss2_1_Uni$lcris_OEC,
                             score = indi_ns_ss2_1_Uni$scoringD2_0)

### Outcome lcris_FTCm12 Uni
rddD2_FTCm12_p1kT_Uni <- rdcate(outcome = indi_ns_ss2_1_Uni$lcris_FTCm12,
                                score = indi_ns_ss2_1_Uni$scoringD2_0)

### Outcome lcris_FTC612 Uni
rddD2_FTC612_p1kT_Uni <- rdcate(outcome = indi_ns_ss2_1_Uni$lcris_FTC612,
                                score = indi_ns_ss2_1_Uni$scoringD2_0)

### Outcome lcris_FTC05 Uni
rddD2_FTC05_p1kT_Uni <- rdcate(outcome = indi_ns_ss2_1_Uni$lcris_FTCme6,
                               score = indi_ns_ss2_1_Uni$scoringD2_0)

quality_education_D2_hetero <- data.frame(
  treatment = c(rep("D2", 2)),
  prOEC_PriLS = c(
    paste0(rddD2_OEC_p1kT_PriLS$pointest, rddD2_OEC_p1kT_PriLS$stars),
    rddD2_OEC_p1kT_PriLS$ci95
  ),
  prFTCm12_PriLS = c(
    paste0(rddD2_FTCm12_p1kT_PriLS$pointest, rddD2_FTCm12_p1kT_PriLS$stars),
    rddD2_FTCm12_p1kT_PriLS$ci95
  ),
  prFTC612_PriLS = c(
    paste0(rddD2_FTC612_p1kT_PriLS$pointest, rddD2_FTC612_p1kT_PriLS$stars),
    rddD2_FTC612_p1kT_PriLS$ci95
  ),
  prFTC05_PriLS = c(
    paste0(rddD2_FTC05_p1kT_PriLS$pointest, rddD2_FTC05_p1kT_PriLS$stars),
    rddD2_FTC05_p1kT_PriLS$ci95
  ),
  prOEC_UppSec = c(
    paste0(rddD2_OEC_p1kT_UppSec$pointest, rddD2_OEC_p1kT_UppSec$stars),
    rddD2_OEC_p1kT_UppSec$ci95
  ),
  prFTCm12_UppSec = c(
    paste0(rddD2_FTCm12_p1kT_UppSec$pointest, rddD2_FTCm12_p1kT_UppSec$stars),
    rddD2_FTCm12_p1kT_UppSec$ci95
  ),
  prFTC612_UppSec = c(
    paste0(rddD2_FTC612_p1kT_UppSec$pointest, rddD2_FTC612_p1kT_UppSec$stars),
    rddD2_FTC612_p1kT_UppSec$ci95
  ),
  prFTC05_UppSec = c(
    paste0(rddD2_FTC05_p1kT_UppSec$pointest, rddD2_FTC05_p1kT_UppSec$stars),
    rddD2_FTC05_p1kT_UppSec$ci95
  ),
  prOEC_Uni = c(
    paste0(rddD2_OEC_p1kT_Uni$pointest, rddD2_OEC_p1kT_Uni$stars),
    rddD2_OEC_p1kT_Uni$ci95
  ),
  prFTCm12_Uni = c(
    paste0(rddD2_FTCm12_p1kT_Uni$pointest, rddD2_FTCm12_p1kT_Uni$stars),
    rddD2_FTCm12_p1kT_Uni$ci95
  ),
  prFTC612_Uni = c(
    paste0(rddD2_FTC612_p1kT_Uni$pointest, rddD2_FTC612_p1kT_Uni$stars),
    rddD2_FTC612_p1kT_Uni$ci95
  ),
  prFTC05_Uni = c(
    paste0(rddD2_FTC05_p1kT_Uni$pointest, rddD2_FTC05_p1kT_Uni$stars),
    rddD2_FTC05_p1kT_Uni$ci95
  )
)


write.xlsx(quality_education_D2_hetero, 'quality_education_D2_hetero.xlsx')

## Depending on date -------------
### D1 --------------


rddD1_type_p1kT_bCovid <- rdcate_multinom(score = indi_ns_ss1_bCovid$scoringD1_0,
                                          outcome = indi_ns_ss1_bCovid$longestcontract_standardo)

rddD1_type_p1kT_afCovid <- rdcate_multinom(score = indi_ns_ss1_afCovid$scoringD1_0,
                                           outcome = indi_ns_ss1_afCovid$longestcontract_standardo)

rddD1_type_p1kT_bCovid
rddD1_type_p1kT_afCovid


### Outcome lcris_OEC
rddD1_OEC_p1kT_after <- rdcate(outcome = indi_ns_ss1_afCovid$lcris_OEC,
                               score = indi_ns_ss1_afCovid$scoringD1_0)

### Outcome lcris_FTCm12
rddD1_FTCm12_p1kT_after <- rdcate(outcome = indi_ns_ss1_afCovid$lcris_FTCm12,
                                  score = indi_ns_ss1_afCovid$scoringD1_0)

### Outcome lcris_FTC612
rddD1_FTC612_p1kT_after <- rdcate(outcome = indi_ns_ss1_afCovid$lcris_FTC612,
                                  score = indi_ns_ss1_afCovid$scoringD1_0)

### Outcome lcris_FTC05
rddD1_FTC05_p1kT_after <- rdcate(outcome = indi_ns_ss1_afCovid$lcris_FTCme6,
                                 score = indi_ns_ss1_afCovid$scoringD1_0)

quality_date_D1_hetero <- data.frame(
  treatment = c(rep("D1", 2)),
  prOEC_after = c(
    paste0(rddD1_OEC_p1kT_after$pointest, rddD1_OEC_p1kT_after$stars),
    rddD1_OEC_p1kT_after$ci95
  ),
  prFTCm12_after = c(
    paste0(rddD1_FTCm12_p1kT_after$pointest, rddD1_FTCm12_p1kT_after$stars),
    rddD1_FTCm12_p1kT_after$ci95
  ),
  prFTC612_after = c(
    paste0(rddD1_FTC612_p1kT_after$pointest, rddD1_FTC612_p1kT_after$stars),
    rddD1_FTC612_p1kT_after$ci95
  ),
  prFTC05_after = c(
    paste0(rddD1_FTC05_p1kT_after$pointest, rddD1_FTC05_p1kT_after$stars),
    rddD1_FTC05_p1kT_after$ci95
  )
)


write.xlsx(quality_date_D1_hetero, 'quality_date_D1_hetero.xlsx')


### D2 --------------

rddD2_type_p1kT_bCovid <- rdcate_multinom(score = indi_ns_ss2_bCovid$scoringD2_0,
                                          outcome = indi_ns_ss2_bCovid$longestcontract_standardo)

rddD2_type_p1kT_afCovid <- rdcate_multinom(score = indi_ns_ss2_afCovid$scoringD2_0,
                                           outcome = indi_ns_ss2_afCovid$longestcontract_standardo)

rddD2_type_p1kT_bCovid
rddD2_type_p1kT_afCovid

### Outcome lcris_OEC
rddD2_OEC_p1kT_after <- rdcate(outcome = indi_ns_ss2_afCovid$lcris_OEC,
                               score = indi_ns_ss2_afCovid$scoringD2_0)

### Outcome lcris_FTCm12
rddD2_FTCm12_p1kT_after <- rdcate(outcome = indi_ns_ss2_afCovid$lcris_FTCm12,
                                  score = indi_ns_ss2_afCovid$scoringD2_0)

### Outcome lcris_FTC612
rddD2_FTC612_p1kT_after <- rdcate(outcome = indi_ns_ss2_afCovid$lcris_FTC612,
                                  score = indi_ns_ss2_afCovid$scoringD2_0)

### Outcome lcris_FTC05
rddD2_FTC05_p1kT_after <- rdcate(outcome = indi_ns_ss2_afCovid$lcris_FTCme6,
                                 score = indi_ns_ss2_afCovid$scoringD2_0)

quality_date_D2_hetero <- data.frame(
  treatment = c(rep("D2", 2)),
  prOEC_after = c(
    paste0(rddD2_OEC_p1kT_after$pointest, rddD2_OEC_p1kT_after$stars),
    rddD2_OEC_p1kT_after$ci95
  ),
  prFTCm12_after = c(
    paste0(rddD2_FTCm12_p1kT_after$pointest, rddD2_FTCm12_p1kT_after$stars),
    rddD2_FTCm12_p1kT_after$ci95
  ),
  prFTC612_after = c(
    paste0(rddD2_FTC612_p1kT_after$pointest, rddD2_FTC612_p1kT_after$stars),
    rddD2_FTC612_p1kT_after$ci95
  ),
  prFTC05_after = c(
    paste0(rddD2_FTC05_p1kT_after$pointest, rddD2_FTC05_p1kT_after$stars),
    rddD2_FTC05_p1kT_after$ci95
  )
)


#write.xlsx(quality_date_D2_hetero, 'quality_date_D2_hetero.xlsx')

## Exporting -------------
### D1 ----------
subgroupD1quali <- bind_rows(list(
  setNames(quality_sex_D1_hetero[1, c(2:3)], c("YOEC", "YFTCm12")),
  setNames(quality_sex_D1_hetero[1, c(6:7)], c("YOEC", "YFTCm12")),
  setNames(quality_age_D1_hetero[1, c(2:3)], c("YOEC", "YFTCm12")),
  setNames(quality_age_D1_hetero[1, c(2:3)], c("YOEC", "YFTCm12")),
  setNames(quality_education_D1_hetero[1, c(2:3)], c("YOEC", "YFTCm12")),
  setNames(quality_education_D1_hetero[1, c(6:7)], c("YOEC", "YFTCm12")),
  setNames(quality_education_D1_hetero[1, c(10:11)], c("YOEC", "YFTCm12")),
  setNames(quality_date_D1_hetero[1, c(2:3)], c("YOEC", "YFTCm12")) ))

subgroups <- subgroupd1$subgroup[nzchar(subgroupd1$subgroup)]

subgroupD1quali$subgroup <- setdiff(subgroups, c("Date < COVID-19"))

subgroupD1all <- left_join(subgroupd1, subgroupD1quali, by = "subgroup")


if (!file.exists("intermediate/script04/B2_pseudotreatments_D1.xlsx")) {

  openxlsx::write.xlsx(subgroupD1all, 'intermediate/script04/B2_pseudotreatments_D1.xlsx')

}



### D2 ---------
subgroupD2quali <- bind_rows(list(
  setNames(quality_sex_D2_hetero[1, c(2:3)], c("YOEC", "YFTCm12")),
  setNames(quality_sex_D2_hetero[1, c(6:7)], c("YOEC", "YFTCm12")),
  setNames(quality_age_D2_hetero[1, c(2:3)], c("YOEC", "YFTCm12")),
  setNames(quality_age_D2_hetero[1, c(2:3)], c("YOEC", "YFTCm12")),
  setNames(quality_education_D2_hetero[1, c(2:3)], c("YOEC", "YFTCm12")),
  setNames(quality_education_D2_hetero[1, c(6:7)], c("YOEC", "YFTCm12")),
  setNames(quality_education_D2_hetero[1, c(10:11)], c("YOEC", "YFTCm12")),
  setNames(quality_date_D2_hetero[1, c(2:3)], c("YOEC", "YFTCm12")) ))

subgroups <- subgroupD2$subgroup[nzchar(subgroupD2$subgroup)]

subgroupD2quali$subgroup <- setdiff(subgroups, c("Date < COVID-19"))

subgroupD2all <- left_join(subgroupD2, subgroupD2quali, by = "subgroup")


if (!file.exists("intermediate/script04/B3_pseudotreatments_D2.xlsx")) {

  openxlsx::write.xlsx(subgroupD2all, 'intermediate/script04/B3_pseudotreatments_D2.xlsx')

}


# 4. Cream skimming (access bias) ------------------

## D1 ---------------

providers1 <- indi_ns_ss1 %>%
  tabyl(axl_ente) %>%
  mutate(vp100 = valid_percent * 100) %>%
  arrange(desc(vp100)) %>%
  mutate(cum_vp100 = cumsum(vp100))

# Binarizing the provider variable
indi_ns_ss1_cPROV <- indi_ns_ss1_c %>% filter(!is.na(axl_ente)) # only 70 NAs out of 22,274 observations

indi_ns_ss1_cPROV$prov1_GI <- ifelse(indi_ns_ss1_cPROV$axl_ente == "GI GROUP  SPA", 1, 0)
indi_ns_ss1_cPROV$prov1_UMANA <- ifelse(indi_ns_ss1_cPROV$axl_ente == "UMANA", 1, 0)
indi_ns_ss1_cPROV$prov1_ADECCO <- ifelse(indi_ns_ss1_cPROV$axl_ente == "ADECCO ITALIA S.P.A.", 1, 0)
indi_ns_ss1_cPROV$prov1_MANPOWER <- ifelse(indi_ns_ss1_cPROV$axl_ente == "Manpower Srl", 1, 0)
indi_ns_ss1_cPROV$prov1_ENAIP <- ifelse(indi_ns_ss1_cPROV$axl_ente == "ENAIP VENETO", 1, 0)
indi_ns_ss1_cPROV$prov1_RANDSTAD <- ifelse(indi_ns_ss1_cPROV$axl_ente == "RANDSTAD ITALIA S.P.A.", 1, 0)
indi_ns_ss1_cPROV$prov1_SYNERGIE <- ifelse(indi_ns_ss1_cPROV$axl_ente == "Synergie Italia Agenzia per il Lavoro S.p.A.", 1, 0)
indi_ns_ss1_cPROV$prov1_STAFF <- ifelse(indi_ns_ss1_cPROV$axl_ente == "Staff S.p.A.", 1, 0)
indi_ns_ss1_cPROV$prov1_ENAC <- ifelse(indi_ns_ss1_cPROV$axl_ente == "Fondazione ENAC Veneto C.F.P. Canossiano", 1, 0)
indi_ns_ss1_cPROV$prov1_INFOLINGUE <- ifelse(indi_ns_ss1_cPROV$axl_ente == "Infolingue S.a.s", 1, 0)
indi_ns_ss1_cPROV$prov1_APIND <- ifelse(indi_ns_ss1_cPROV$axl_ente == "APINDUSTRIA SERVIZI SRL", 1, 0)
indi_ns_ss1_cPROV$prov1_EUROINT <- ifelse(indi_ns_ss1_cPROV$axl_ente == "Eurointerim s.p.a.", 1, 0)

# Model estimation

rddD1_GI <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_GI,
                   score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_UMANA <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_UMANA,
                      score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_ADECCO <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_ADECCO,
                       score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_MANPOWER <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_MANPOWER,
                         score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_ENAIP <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_ENAIP,
                      score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_RANDSTAD <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_RANDSTAD,
                         score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_SYNERGIE <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_SYNERGIE,
                         score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_STAFF <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_STAFF,
                      score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_ENAC <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_ENAC,
                     score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_INFOLINGUE <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_INFOLINGUE,
                           score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_APIND <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_APIND,
                      score = indi_ns_ss1_cPROV$scoringD1_0)

rddD1_EUROINT <- rdcate(outcome = indi_ns_ss1_cPROV$prov1_EUROINT,
                        score = indi_ns_ss1_cPROV$scoringD1_0)

# Export

cream_D1 <- data.frame(
  treatment = c(rep("D1", 2)),
  GI = c(
    paste0(rddD1_GI$pointest, rddD1_GI$stars),
    rddD1_GI$ci95
  ),
  UMANA = c(
    paste0(rddD1_UMANA$pointest, rddD1_UMANA$stars),
    rddD1_UMANA$ci95
  ),
  ADECCO = c(
    paste0(rddD1_ADECCO$pointest, rddD1_ADECCO$stars),
    rddD1_ADECCO$ci95
  ),
  MANPOWER = c(
    paste0(rddD1_MANPOWER$pointest, rddD1_MANPOWER$stars),
    rddD1_MANPOWER$ci95
  ),
  ENAIP = c(
    paste0(rddD1_ENAIP$pointest, rddD1_ENAIP$stars),
    rddD1_ENAIP$ci95
  ),
  RANDSTAD = c(
    paste0(rddD1_RANDSTAD$pointest, rddD1_RANDSTAD$stars),
    rddD1_RANDSTAD$ci95
  ),
  SYNERGIE = c(
    paste0(rddD1_SYNERGIE$pointest, rddD1_SYNERGIE$stars),
    rddD1_SYNERGIE$ci95
  ),
  STAFF = c(
    paste0(rddD1_STAFF$pointest, rddD1_STAFF$stars),
    rddD1_STAFF$ci95
  ),
  ENAC = c(
    paste0(rddD1_ENAC$pointest, rddD1_ENAC$stars),
    rddD1_ENAC$ci95
  ),
  INFOLINGUE = c(
    paste0(rddD1_INFOLINGUE$pointest, rddD1_INFOLINGUE$stars),
    rddD1_INFOLINGUE$ci95
  ),
  APIND = c(
    paste0(rddD1_APIND$pointest, rddD1_APIND$stars),
    rddD1_APIND$ci95
  ),
  EUROINT = c(
    paste0(rddD1_EUROINT$pointest, rddD1_EUROINT$stars),
    rddD1_EUROINT$ci95
  )
)


write.xlsx(cream_D1, 'cream_D1.xlsx')


## D2 ---------------

providers2 <- indi_ns_ss2 %>%
  tabyl(axl_ente) %>%
  mutate(vp100 = valid_percent * 100) %>%
  arrange(desc(vp100)) %>%
  mutate(cum_vp100 = cumsum(vp100))

# Binarizing the provider variable
indi_ns_ss2_cPROV <- subset(indi_ns_ss2_c, !is.na(indi_ns_ss2_c$axl_ente))
sum(is.na(indi_ns_ss2_c$axl_ente)) # only 77 NAs out of 20,953 observations

indi_ns_ss2_cPROV$prov2_GI <- ifelse(indi_ns_ss2_cPROV$axl_ente == "GI GROUP  SPA", 1, 0)
indi_ns_ss2_cPROV$prov2_UMANA <- ifelse(indi_ns_ss2_cPROV$axl_ente == "UMANA", 1, 0)
indi_ns_ss2_cPROV$prov2_ADECCO <- ifelse(indi_ns_ss2_cPROV$axl_ente == "ADECCO ITALIA S.P.A.", 1, 0)
indi_ns_ss2_cPROV$prov2_ENAIP <- ifelse(indi_ns_ss2_cPROV$axl_ente == "ENAIP VENETO", 1, 0)
indi_ns_ss2_cPROV$prov2_MANPOWER <- ifelse(indi_ns_ss2_cPROV$axl_ente == "Manpower Srl", 1, 0)
indi_ns_ss2_cPROV$prov2_RANDSTAD <- ifelse(indi_ns_ss2_cPROV$axl_ente == "RANDSTAD ITALIA S.P.A.", 1, 0)
indi_ns_ss2_cPROV$prov2_SYNERGIE <- ifelse(indi_ns_ss2_cPROV$axl_ente == "Synergie Italia Agenzia per il Lavoro S.p.A.", 1, 0)
indi_ns_ss2_cPROV$prov2_STAFF <- ifelse(indi_ns_ss2_cPROV$axl_ente == "Staff S.p.A.", 1, 0)
indi_ns_ss2_cPROV$prov2_ENAC <- ifelse(indi_ns_ss2_cPROV$axl_ente == "Fondazione ENAC Veneto C.F.P. Canossiano", 1, 0)
indi_ns_ss2_cPROV$prov2_EUROINT <- ifelse(indi_ns_ss2_cPROV$axl_ente == "Eurointerim s.p.a.", 1, 0)
indi_ns_ss2_cPROV$prov2_APIND <- ifelse(indi_ns_ss2_cPROV$axl_ente == "APINDUSTRIA SERVIZI SRL", 1, 0)
indi_ns_ss2_cPROV$prov2_EDUFORMA <- ifelse(indi_ns_ss2_cPROV$axl_ente == "Eduforma Srl", 1, 0)
indi_ns_ss2_cPROV$prov2_INFOLINGUE <- ifelse(indi_ns_ss2_cPROV$axl_ente == "Infolingue S.a.s", 1, 0)

# Model estimation

rddD2_GI <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_GI,
                   score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_UMANA <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_UMANA,
                      score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_ADECCO <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_ADECCO,
                       score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_ENAIP <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_ENAIP,
                      score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_MANPOWER <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_MANPOWER,
                         score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_RANDSTAD <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_RANDSTAD,
                         score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_SYNERGIE <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_SYNERGIE,
                         score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_STAFF <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_STAFF,
                      score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_ENAC <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_ENAC,
                     score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_EUROINT <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_EUROINT,
                        score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_APIND <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_APIND,
                      score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_EDUFORMA <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_EDUFORMA,
                         score = indi_ns_ss2_cPROV$scoringD2_0)

rddD2_INFOLINGUE <- rdcate(outcome = indi_ns_ss2_cPROV$prov2_INFOLINGUE,
                           score = indi_ns_ss2_cPROV$scoringD2_0)





# Export

cream_D2 <- data.frame(
  treatment = c(rep("D2", 2)),
  GI = c(
    paste0(rddD2_GI$pointest, rddD2_GI$stars),
    rddD2_GI$ci95
  ),
  UMANA = c(
    paste0(rddD2_UMANA$pointest, rddD2_UMANA$stars),
    rddD2_UMANA$ci95
  ),
  ADECCO = c(
    paste0(rddD2_ADECCO$pointest, rddD2_ADECCO$stars),
    rddD2_ADECCO$ci95
  ),
  ENAIP = c(
    paste0(rddD2_ENAIP$pointest, rddD2_ENAIP$stars),
    rddD2_ENAIP$ci95
  ),
  MANPOWER = c(
    paste0(rddD2_MANPOWER$pointest, rddD2_MANPOWER$stars),
    rddD2_MANPOWER$ci95
  ),
  RANDSTAD = c(
    paste0(rddD2_RANDSTAD$pointest, rddD2_RANDSTAD$stars),
    rddD2_RANDSTAD$ci95
  ),
  SYNERGIE = c(
    paste0(rddD2_SYNERGIE$pointest, rddD2_SYNERGIE$stars),
    rddD2_SYNERGIE$ci95
  ),
  STAFF = c(
    paste0(rddD2_STAFF$pointest, rddD2_STAFF$stars),
    rddD2_STAFF$ci95
  ),
  ENAC = c(
    paste0(rddD2_ENAC$pointest, rddD2_ENAC$stars),
    rddD2_ENAC$ci95
  ),
  EUROINT = c(
    paste0(rddD2_EUROINT$pointest, rddD2_EUROINT$stars),
    rddD2_EUROINT$ci95
  ),
  APIND = c(
    paste0(rddD2_APIND$pointest, rddD2_APIND$stars),
    rddD2_APIND$ci95
  ),
  EDUFORMA = c(
    paste0(rddD2_EDUFORMA$pointest, rddD2_EDUFORMA$stars),
    rddD2_EDUFORMA$ci95
  ),
  INFOLINGUE = c(
    paste0(rddD2_INFOLINGUE$pointest, rddD2_INFOLINGUE$stars),
    rddD2_INFOLINGUE$ci95
  )
)


write.xlsx(cream_D2, 'cream_D2.xlsx')

# Is it serious the problem of creaming evidence for Synergie?
655 / (nrow(indi_ns_ss1_cPROV) + nrow(indi_ns_ss2_cPROV)) * 100
