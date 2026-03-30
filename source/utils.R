# Forest + tree + subgroup analyses by final leaves
forest_and_tree <- function(outcome, running, cutoff, df,
                            vector_y, matrix_treatrunning, matrix_x,
                            seed = 7778, treedepth = 2){
  set.seed(seed)

  # 1. Growing the forest and getting predictions
  bw <- rdrobust::rdbwselect(y = df[[outcome]],
                             x = df[[running]],
                             c = cutoff)$bws[1,1]

  distance <- abs((df[[running]] - cutoff) / bw)
  tkernel.weights <- (1 - distance) * (distance <= 1) / bw # triangular kernel
  windowsubset <- tkernel.weights > 0

  forest <- lm_forest(Y = vector_y[windowsubset],
                      W = matrix_treatrunning[windowsubset, ], # W (K=1) is treatment indicator, Z (K=2) is running variable
                      X = matrix_x[windowsubset, ],
                      sample.weights = tkernel.weights[windowsubset],
                      gradient.weights = c(1, 0),
                      seed = seed)

  message("Forest grown.\n")

  # 2. Running the model: getting predictions at individual and mean level
  tau_hat <- predict(forest)$predictions[, 1, ] # second dimension asks for predictions for the coefficient of W (K=1)
  easycost <- mean(tau_hat)
  reward <- cbind(control = -tau_hat, treat = tau_hat - easycost)

  # 3. Growing the policy tree with a training set
  matrix_xIN <- matrix_x[tkernel.weights > 0, ]

  set_train <- sample(x = 1:nrow(matrix_xIN), size = nrow(matrix_xIN) * 0.6)
  set_test <- setdiff(1:nrow(matrix_xIN), set_train)

  p_tree <- fastpolicytree(X = matrix_xIN[set_train, ],
                           Gamma = reward[set_train, ],
                           depth = treedepth)

  #plot(p_tree, leaf.labels = c("Group A", "Group B"))

  message("Policy tree grown.\n")

  # 4. Evaluating the policy tree
  ## 4.1. Gains of following the recommended assignment ("policy value")
  pi_hat_test <- predict(p_tree, matrix_xIN[set_test, ]) - 1  # Note policytree recodes the treatments to 1,2.
  # We substract one to get back to our usual encoding 0,1.

  gains_assignment <- mean((tau_hat[set_test] - easycost) * (2*pi_hat_test - 1))

  ## 4.2. Inference through proper procedures

  if(treedepth == 2){

    subgroups <- data.frame(pair = rep(NA, 2),
                            var1 = rep(NA, 2),
                            var1_thr = rep(NA, 2),
                            var2 = rep(NA, 2),
                            var2_thr = rep(NA, 2),
                            diff_estimate = rep(NA, 2),
                            diff_pvalue = rep(NA, 2),
                            Nh_left = rep(NA, 2),
                            Nh_right = rep(NA, 2))

    # Pair 1
    ## Prepare data
    firstvar_index <- p_tree$nodes[[1]]$split_variable
    firstvar <- colnames(matrix_xIN)[firstvar_index]
    firstvalue <- p_tree$nodes[[1]]$split_value

    secondvar_index <- p_tree$nodes[[2]]$split_variable
    secondvar <- colnames(matrix_xIN)[secondvar_index]
    secondvalue <- p_tree$nodes[[2]]$split_value

    df_designmatrix <- as.data.frame(cbind(vector_y, matrix_treatrunning, matrix_x))
    colnames(df_designmatrix)[1] <- outcome
    colnames(df_designmatrix)[3] <- running

    group_left_leaf <- if(p_tree$nodes[4][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups <- df_designmatrix |>
      filter(!!sym(firstvar) <= firstvalue) |>
      mutate(leaf = if_else(!!sym(secondvar) <= secondvalue, group_left_leaf, group_right_leaf))

    df_subgroups$leaf <- factor(df_subgroups$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    ## by chance, the test set could have 0 observations in one of the subgroups...
    if(sum(df_subgroups$leaf == "Control") > 0 & sum(df_subgroups$leaf == "Treat") > 0){

      # in small subgroups, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      if(any(table(df_subgroups$leaf, df_subgroups[[running]] >= 0) < 5)){

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf,
                         bw.joint = T)

      }else{

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf)

        rd_hypothesis <- c("`df_subgroups$leafTreat` - `df_subgroups$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups$diff_estimate[1] <- rdhypo$individual$estimate
        subgroups$diff_pvalue[1] <- rdhypo$individual$p_value

      }

      subgroups$Nh_left[1] <- sum(rdmodel$Nh[, 1])
      subgroups$Nh_right[1] <- sum(rdmodel$Nh[, 2])

    }

    ## Save
    subgroups$pair[1] <- 1
    subgroups$var1[1] <- firstvar
    subgroups$var1_thr[1] <- paste0("<= ", firstvalue)
    subgroups$var2[1] <- secondvar
    subgroups$var2_thr[1] <- secondvalue

    # Pair 2
    ## Prepare data
    secondvar_index <- p_tree$nodes[[3]]$split_variable
    secondvar <- colnames(matrix_xIN)[secondvar_index]
    secondvalue <- p_tree$nodes[[3]]$split_value

    group_left_leaf <- if(p_tree$nodes[6][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups <- df_designmatrix |>
      filter(!!sym(firstvar) > firstvalue) |> # Note the >
      mutate(leaf = if_else(!!sym(secondvar) <= secondvalue, group_left_leaf, group_right_leaf))

    df_subgroups$leaf <- factor(df_subgroups$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    if(sum(df_subgroups$leaf == "Control") > 0 & sum(df_subgroups$leaf == "Treat") > 0){

      # in small subgroups, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      if(any(table(df_subgroups$leaf, df_subgroups[[running]] >= 0) < 5)){

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf,
                         bw.joint = T)

      }else{

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf)

        rd_hypothesis <- c("`df_subgroups$leafTreat` - `df_subgroups$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups$diff_estimate[2] <- rdhypo$individual$estimate
        subgroups$diff_pvalue[2] <- rdhypo$individual$p_value

      }

      subgroups$Nh_left[2] <- sum(rdmodel$Nh[, 1])
      subgroups$Nh_right[2] <- sum(rdmodel$Nh[, 2])

    }

    ## Save
    subgroups$pair[2] <- 2
    subgroups$var1[2] <- firstvar
    subgroups$var1_thr[2] <- paste0("> ", firstvalue)
    subgroups$var2[2] <- secondvar
    subgroups$var2_thr[2] <- secondvalue

  }

  if(treedepth == 3){

    subgroups <- data.frame(pair = rep(NA, 4),
                            var1 = rep(NA, 4),
                            var1_thr = rep(NA, 4),
                            var2 = rep(NA, 4),
                            var2_thr = rep(NA, 4),
                            var3 = rep(NA, 4),
                            var3_thr = rep(NA, 4),
                            diff_estimate = rep(NA, 4),
                            diff_pvalue = rep(NA, 4),
                            Nh_left = rep(NA, 4),
                            Nh_right = rep(NA, 4))

    # Pair 1
    ## Prepare data
    firstvar_index <- p_tree$nodes[[1]]$split_variable
    firstvar <- colnames(matrix_xIN)[firstvar_index]
    firstvalue <- p_tree$nodes[[1]]$split_value

    secondvar_index <- p_tree$nodes[[2]]$split_variable
    secondvar <- colnames(matrix_xIN)[secondvar_index]
    secondvalue <- p_tree$nodes[[2]]$split_value

    thirdvar_index <- p_tree$nodes[[4]]$split_variable
    thirdvar <- colnames(matrix_xIN)[thirdvar_index]
    thirdvalue <- p_tree$nodes[[4]]$split_value

    df_designmatrix <- as.data.frame(cbind(vector_y, matrix_treatrunning, matrix_x))
    colnames(df_designmatrix)[1] <- outcome
    colnames(df_designmatrix)[3] <- running

    group_left_leaf <- if(p_tree$nodes[8][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups <- df_designmatrix |>
      filter(!!sym(firstvar) <= firstvalue) |>
      filter(!!sym(secondvar) <= secondvalue) |>
      mutate(leaf = if_else(!!sym(thirdvar) <= thirdvalue, group_left_leaf, group_right_leaf))

    df_subgroups$leaf <- factor(df_subgroups$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    ## by chance, the test set could have 0 observations in one of the subgroups...
    if(sum(df_subgroups$leaf == "Control") > 0 & sum(df_subgroups$leaf == "Treat") > 0){

      # in small subgroups, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      if(any(table(df_subgroups$leaf, df_subgroups[[running]] >= 0) < 5)){

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf,
                         bw.joint = T)

      }else{

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf)

        rd_hypothesis <- c("`df_subgroups$leafTreat` - `df_subgroups$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups$diff_estimate[1] <- rdhypo$individual$estimate
        subgroups$diff_pvalue[1] <- rdhypo$individual$p_value

      }

      subgroups$Nh_left[1] <- sum(rdmodel$Nh[, 1])
      subgroups$Nh_right[1] <- sum(rdmodel$Nh[, 2])
    }


    ## Save
    subgroups$pair[1] <- 1
    subgroups$var1[1] <- firstvar
    subgroups$var1_thr[1] <- paste0("<= ", firstvalue)
    subgroups$var2[1] <- secondvar
    subgroups$var2_thr[1] <- paste0("<= ", secondvalue)
    subgroups$var3[1] <- thirdvar
    subgroups$var3_thr[1] <- thirdvalue


    # Pair 2
    ## Prepare data
    thirdvar_index <- p_tree$nodes[[5]]$split_variable
    thirdvar <- colnames(matrix_xIN)[thirdvar_index]
    thirdvalue <- p_tree$nodes[[5]]$split_value

    group_left_leaf <- if(p_tree$nodes[10][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups <- df_designmatrix |>
      filter(!!sym(firstvar) <= firstvalue) |>
      filter(!!sym(secondvar) > secondvalue) |> # Note the >
      mutate(leaf = if_else(!!sym(thirdvar) <= thirdvalue, group_left_leaf, group_right_leaf))

    df_subgroups$leaf <- factor(df_subgroups$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    if(sum(df_subgroups$leaf == "Control") > 0 & sum(df_subgroups$leaf == "Treat") > 0){

      # in small subgroups, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      if(any(table(df_subgroups$leaf, df_subgroups[[running]] >= 0) < 5)){

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf,
                         bw.joint = T)

      }else{

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf)

        rd_hypothesis <- c("`df_subgroups$leafTreat` - `df_subgroups$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups$diff_estimate[2] <- rdhypo$individual$estimate
        subgroups$diff_pvalue[2] <- rdhypo$individual$p_value

      }

      subgroups$Nh_left[2] <- sum(rdmodel$Nh[, 1])
      subgroups$Nh_right[2] <- sum(rdmodel$Nh[, 2])
    }


    ## Save
    subgroups$pair[2] <- 2
    subgroups$var1[2] <- firstvar
    subgroups$var1_thr[2] <- paste0("<= ", firstvalue)
    subgroups$var2[2] <- secondvar
    subgroups$var2_thr[2] <- paste0("> ", secondvalue)
    subgroups$var3[2] <- thirdvar
    subgroups$var3_thr[2] <- thirdvalue

    # Pair 3
    ## Prepare data
    secondvar_index <- p_tree$nodes[[3]]$split_variable
    secondvar <- colnames(matrix_xIN)[secondvar_index]
    secondvalue <- p_tree$nodes[[3]]$split_value

    thirdvar_index <- p_tree$nodes[[6]]$split_variable
    thirdvar <- colnames(matrix_xIN)[thirdvar_index]
    thirdvalue <- p_tree$nodes[[6]]$split_value

    group_left_leaf <- if(p_tree$nodes[12][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups <- df_designmatrix |>
      filter(!!sym(firstvar) > firstvalue) |> # Note the >
      filter(!!sym(secondvar) <= secondvalue) |>
      mutate(leaf = if_else(!!sym(thirdvar) <= thirdvalue, group_left_leaf, group_right_leaf))

    df_subgroups$leaf <- factor(df_subgroups$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    if(sum(df_subgroups$leaf == "Control") > 0 & sum(df_subgroups$leaf == "Treat") > 0){

      # in small subgroups, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      if(any(table(df_subgroups$leaf, df_subgroups[[running]] >= 0) < 5)){

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf,
                         bw.joint = T)

      }else{

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf)

        rd_hypothesis <- c("`df_subgroups$leafTreat` - `df_subgroups$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups$diff_estimate[3] <- rdhypo$individual$estimate
        subgroups$diff_pvalue[3] <- rdhypo$individual$p_value

      }

      subgroups$Nh_left[3] <- sum(rdmodel$Nh[, 1])
      subgroups$Nh_right[3] <- sum(rdmodel$Nh[, 2])

    }


    ## Save
    subgroups$pair[3] <- 3
    subgroups$var1[3] <- firstvar
    subgroups$var1_thr[3] <- paste0("> ", firstvalue)
    subgroups$var2[3] <- secondvar
    subgroups$var2_thr[3] <- paste0("<= ", secondvalue)
    subgroups$var3[3] <- thirdvar
    subgroups$var3_thr[3] <- thirdvalue

    # Pair 4
    ## Prepare data
    thirdvar_index <- p_tree$nodes[[7]]$split_variable
    thirdvar <- colnames(matrix_xIN)[thirdvar_index]
    thirdvalue <- p_tree$nodes[[7]]$split_value

    group_left_leaf <- if(p_tree$nodes[14][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups <- df_designmatrix |>
      filter(!!sym(firstvar) > firstvalue) |> # Note the >
      filter(!!sym(secondvar) > secondvalue) |> # Note the >
      mutate(leaf = if_else(!!sym(thirdvar) <= thirdvalue, group_left_leaf, group_right_leaf))

    df_subgroups$leaf <- factor(df_subgroups$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    if(sum(df_subgroups$leaf == "Control") > 0 & sum(df_subgroups$leaf == "Treat") > 0){

      # in small subgroups, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      if(any(table(df_subgroups$leaf, df_subgroups[[running]] >= 0) < 5)){

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf,
                         bw.joint = T)

      }else{

        rdmodel <- rdhte(y = df_subgroups[[outcome]],
                         x = df_subgroups[[running]],
                         covs.hte = df_subgroups$leaf)

        rd_hypothesis <- c("`df_subgroups$leafTreat` - `df_subgroups$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups$diff_estimate[4] <- rdhypo$individual$estimate
        subgroups$diff_pvalue[4] <- rdhypo$individual$p_value

      }

      subgroups$Nh_left[4] <- sum(rdmodel$Nh[, 1])
      subgroups$Nh_right[4] <- sum(rdmodel$Nh[, 2])
    }


    ## Save
    subgroups$pair[4] <- 4
    subgroups$var1[4] <- firstvar
    subgroups$var1_thr[4] <- paste0("> ", firstvalue)
    subgroups$var2[4] <- secondvar
    subgroups$var2_thr[4] <- paste0("> ", secondvalue)
    subgroups$var3[4] <- thirdvar
    subgroups$var3_thr[4] <- thirdvalue


  }

  return(list(gains_assignment,
              subgroups,
              p_tree))

}

# Auxiliary function to check if sample size with optimal bandwidth is too low at any side of the cutoff
tiny_sample_size <- function(dataframe, outcome, running, leafvariable){

  # A. With different bandwidths at each side
  optimalbw_2 <- rdbwhte(y = dataframe[[outcome]],
                         x = dataframe[[running]],
                         covs.hte = dataframe[[leafvariable]])

  result <- if(any(optimalbw_2$Nh < 5) == T){"Tiny with two bws"}else{"Not tiny"}


  if(result == "Tiny with two bws"){

    # B. With common bandwidth
    optimalbw_1 <- rdbwhte(y = dataframe[[outcome]],
                           x = dataframe[[running]],
                           covs.hte = dataframe[[leafvariable]],
                           bw.joint = T)

    result <- if(any(optimalbw_1$Nh < 5) == T){"Tiny with one bw"}

  }

  return(result)

}

# Forest + tree + subgroup analysis at each depth level
forest_tree <- function(outcome, running, cutoff, df,
                            vector_y, matrix_treatrunning, matrix_x,
                            seed = 7778, treedepth = 2){
  set.seed(seed)

  # 1. Growing the forest and getting predictions
  bw <- rdrobust::rdbwselect(y = df[[outcome]],
                             x = df[[running]],
                             c = cutoff)$bws[1,1]

  distance <- abs((df[[running]] - cutoff) / bw)
  tkernel.weights <- (1 - distance) * (distance <= 1) / bw # triangular kernel
  windowsubset <- tkernel.weights > 0

  forest <- lm_forest(Y = vector_y[windowsubset],
                      W = matrix_treatrunning[windowsubset, ], # W (K=1) is treatment indicator, Z (K=2) is running variable
                      X = matrix_x[windowsubset, ],
                      sample.weights = tkernel.weights[windowsubset],
                      gradient.weights = c(1, 0),
                      seed = seed)

  message("Forest grown.\n")

  # 2. Running the model: getting predictions at individual and mean level
  tau_hat <- predict(forest)$predictions[, 1, ] # second dimension asks for predictions for the coefficient of W (K=1)
  y.hat.lmf <- forest$Y.hat

  #easycost <- mean(tau_hat)
  #reward <- cbind(control = -tau_hat, treat = tau_hat - easycost)

  cost <- mean(tau_hat)

  rw_df <- data.frame(actual_w = matrix_treatrunning[windowsubset, 1],
                      tau_hat = tau_hat,
                      y_hat = as.numeric(y.hat.lmf))

  rw_df <- rw_df |>
    mutate(y_hat_W0 = if_else(actual_w == 1, y_hat - tau_hat, y_hat))

  reward <- cbind(control = rw_df$y_hat_W0,
                  treat = rw_df$tau_hat - cost + rw_df$y_hat_W0)

  # 3. Growing the policy tree with a training set
  matrix_xIN <- matrix_x[tkernel.weights > 0, ]

  set_train <- sample(x = 1:nrow(matrix_xIN), size = nrow(matrix_xIN) * 0.6)
  set_test <- setdiff(1:nrow(matrix_xIN), set_train)

  p_tree <- fastpolicytree(X = matrix_xIN[set_train, ],
                           Gamma = reward[set_train, ],
                           depth = treedepth)

  #plot(p_tree, leaf.labels = c("Group A", "Group B"))

  message("Policy tree grown.\n")

  # 4. Evaluating the policy tree
  ## 4.1. Gains of following the recommended assignment ("policy value")
  pi_hat_test <- predict(p_tree, matrix_xIN[set_test, ]) - 1  # Note policytree recodes the treatments to 1,2.
  # We substract one to get back to our usual encoding 0,1.

  gains_assignment <- mean((tau_hat[set_test] - cost) * (2*pi_hat_test - 1))

  ## 4.2. Inference through proper procedures

  if(treedepth == 2){
    # DEPTH LEVEL = 1 **********************************************************
    subgroups1 <- data.frame(pair = rep(NA, 1),
                             var1 = rep(NA, 1),
                             var1_thr = rep(NA, 1),
                             var2 = rep(NA, 1),
                             var2_thr = rep(NA, 1),
                             diff_estimate = rep(NA, 1),
                             diff_pvalue = rep(NA, 1),
                             Nh_left = rep(NA, 1),
                             Nh_right = rep(NA, 1))

    # Pair
    ## Prepare data
    firstvar_index <- p_tree$nodes[[1]]$split_variable
    firstvar <- colnames(matrix_xIN)[firstvar_index]
    firstvalue <- p_tree$nodes[[1]]$split_value

    df_designmatrix <- as.data.frame(cbind(vector_y, matrix_treatrunning, matrix_x))
    colnames(df_designmatrix)[1] <- outcome
    colnames(df_designmatrix)[3] <- running

    #group_left_leaf1 <- if(p_tree$nodes[4][[1]]$action == 2){"Treat"}else{"Control"}
    #group_right_leaf1 <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups1 <- df_designmatrix |>
      #filter(!!sym(firstvar) <= firstvalue) |>
      mutate(leaf = if_else(!!sym(firstvar) <= firstvalue, "LoE", "Larger"))

    df_subgroups1$leaf <- factor(df_subgroups1$leaf, levels = c("LoE", "Larger"))

    ## Run RD analysis

    ## by chance, the test set could have 0 observations in one of the subgroups2...
    if(sum(df_subgroups1$leaf == "LoE") > 0 & sum(df_subgroups1$leaf == "Larger") > 0){

      # in small subgroups2, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      if(any(table(df_subgroups1$leaf, df_subgroups1[[running]] >= 0) < 5)){

        rdmodel <- rdhte(y = df_subgroups1[[outcome]],
                         x = df_subgroups1[[running]],
                         covs.hte = df_subgroups1$leaf,
                         bw.joint = T)

      }else{

        rdmodel <- rdhte(y = df_subgroups1[[outcome]],
                         x = df_subgroups1[[running]],
                         covs.hte = df_subgroups1$leaf)

        rd_hypothesis <- c("`df_subgroups1$leafLarger` - `df_subgroups1$leafLoE` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups1$diff_estimate[1] <- rdhypo$individual$estimate
        subgroups1$diff_pvalue[1] <- rdhypo$individual$p_value

      }

      subgroups1$Nh_left[1] <- sum(rdmodel$Nh[, 1])
      subgroups1$Nh_right[1] <- sum(rdmodel$Nh[, 2])

    }

    ## Save
    subgroups1$var1[1] <- firstvar
    subgroups1$var1_thr[1] <- paste0("<= ", firstvalue)

    # DEPTH LEVEL = 2 **********************************************************
    subgroups2 <- data.frame(pair = rep(NA, 2),
                            var1 = rep(NA, 2),
                            var1_thr = rep(NA, 2),
                            var2 = rep(NA, 2),
                            var2_thr = rep(NA, 2),
                            diff_estimate = rep(NA, 2),
                            diff_pvalue = rep(NA, 2),
                            Nh_left = rep(NA, 2),
                            Nh_right = rep(NA, 2))

    # Pair 1
    ## Prepare data
    secondvar_index <- p_tree$nodes[[2]]$split_variable
    secondvar <- colnames(matrix_xIN)[secondvar_index]
    secondvalue <- p_tree$nodes[[2]]$split_value

    #df_designmatrix <- as.data.frame(cbind(vector_y, matrix_treatrunning, matrix_x))
    #colnames(df_designmatrix)[1] <- outcome
    #colnames(df_designmatrix)[3] <- running

    group_left_leaf <- if(p_tree$nodes[4][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups2 <- df_designmatrix |>
      filter(!!sym(firstvar) <= firstvalue) |>
      mutate(leaf = if_else(!!sym(secondvar) <= secondvalue, group_left_leaf, group_right_leaf))

    df_subgroups2$leaf <- factor(df_subgroups2$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    ## by chance, the test set could have 0 observations in one of the subgroups2...
    if(sum(df_subgroups2$leaf == "Control") > 0 & sum(df_subgroups2$leaf == "Treat") > 0){

      # in small subgroups2, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      if(any(table(df_subgroups2$leaf, df_subgroups2[[running]] >= 0) < 5)){

        rdmodel <- rdhte(y = df_subgroups2[[outcome]],
                         x = df_subgroups2[[running]],
                         covs.hte = df_subgroups2$leaf,
                         bw.joint = T)

      }else{

        rdmodel <- rdhte(y = df_subgroups2[[outcome]],
                         x = df_subgroups2[[running]],
                         covs.hte = df_subgroups2$leaf)

        rd_hypothesis <- c("`df_subgroups2$leafTreat` - `df_subgroups2$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups2$diff_estimate[1] <- rdhypo$individual$estimate
        subgroups2$diff_pvalue[1] <- rdhypo$individual$p_value

      }

      subgroups2$Nh_left[1] <- sum(rdmodel$Nh[, 1])
      subgroups2$Nh_right[1] <- sum(rdmodel$Nh[, 2])

    }

    ## Save
    subgroups2$pair[1] <- 1
    subgroups2$var1[1] <- firstvar
    subgroups2$var1_thr[1] <- paste0("<= ", firstvalue)
    subgroups2$var2[1] <- secondvar
    subgroups2$var2_thr[1] <- secondvalue

    # Pair 2
    ## Prepare data
    secondvar_index <- p_tree$nodes[[3]]$split_variable
    secondvar <- colnames(matrix_xIN)[secondvar_index]
    secondvalue <- p_tree$nodes[[3]]$split_value

    group_left_leaf <- if(p_tree$nodes[6][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups2 <- df_designmatrix |>
      filter(!!sym(firstvar) > firstvalue) |> # Note the >
      mutate(leaf = if_else(!!sym(secondvar) <= secondvalue, group_left_leaf, group_right_leaf))

    df_subgroups2$leaf <- factor(df_subgroups2$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    if(sum(df_subgroups2$leaf == "Control") > 0 & sum(df_subgroups2$leaf == "Treat") > 0){

      # in small subgroups2, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      if(any(table(df_subgroups2$leaf, df_subgroups2[[running]] >= 0) < 5)){

        rdmodel <- rdhte(y = df_subgroups2[[outcome]],
                         x = df_subgroups2[[running]],
                         covs.hte = df_subgroups2$leaf,
                         bw.joint = T)

      }else{

        rdmodel <- rdhte(y = df_subgroups2[[outcome]],
                         x = df_subgroups2[[running]],
                         covs.hte = df_subgroups2$leaf)

        rd_hypothesis <- c("`df_subgroups2$leafTreat` - `df_subgroups2$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups2$diff_estimate[2] <- rdhypo$individual$estimate
        subgroups2$diff_pvalue[2] <- rdhypo$individual$p_value

      }

      subgroups2$Nh_left[2] <- sum(rdmodel$Nh[, 1])
      subgroups2$Nh_right[2] <- sum(rdmodel$Nh[, 2])

    }

    ## Save
    subgroups2$pair[2] <- 2
    subgroups2$var1[2] <- firstvar
    subgroups2$var1_thr[2] <- paste0("> ", firstvalue)
    subgroups2$var2[2] <- secondvar
    subgroups2$var2_thr[2] <- secondvalue

  }


  if(treedepth == 3){

    # DEPTH LEVEL = 1 **********************************************************
    subgroups1 <- data.frame(pair = rep(NA, 1),
                             var1 = rep(NA, 1),
                             var1_thr = rep(NA, 1),
                             var2 = rep(NA, 1),
                             var2_thr = rep(NA, 1),
                             var3 = rep(NA, 1),
                             var3_thr = rep(NA, 1),
                             diff_estimate = rep(NA, 1),
                             diff_pvalue = rep(NA, 1),
                             Nh_left = rep(NA, 1),
                             Nh_right = rep(NA, 1))

    # Pair
    ## Prepare data
    firstvar_index <- p_tree$nodes[[1]]$split_variable
    firstvar <- colnames(matrix_xIN)[firstvar_index]
    firstvalue <- p_tree$nodes[[1]]$split_value

    df_designmatrix <- as.data.frame(cbind(vector_y, matrix_treatrunning, matrix_x))
    colnames(df_designmatrix)[1] <- outcome
    colnames(df_designmatrix)[3] <- running

    df_subgroups1 <- df_designmatrix |>
      mutate(leaf = if_else(!!sym(firstvar) <= firstvalue, "LoE", "Larger"))

    df_subgroups1$leaf <- factor(df_subgroups1$leaf, levels = c("LoE", "Larger"))

    ## Run RD analysis

    ## by chance, the test set could have 0 observations in one of the subgroups2...
    if(sum(df_subgroups1$leaf == "LoE") > 0 & sum(df_subgroups1$leaf == "Larger") > 0){

      # in small subgroups2, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      tss <- tiny_sample_size(df_subgroups1, outcome, running, "leaf")
      if(tss == "Tiny with two bws"){

        rdmodel <- rdhte(y = df_subgroups1[[outcome]],
                         x = df_subgroups1[[running]],
                         covs.hte = df_subgroups1$leaf,
                         bw.joint = T)

        rd_hypothesis <- c("`df_subgroups1$leafLarger` - `df_subgroups1$leafLoE` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups1$diff_estimate[1] <- rdhypo$individual$estimate
        subgroups1$diff_pvalue[1] <- rdhypo$individual$p_value

      }

      if(tss == "Tiny with one bw"){

        message("There was a pair of subgroups with a tiny sample size!")

      }

      else{

        rdmodel <- rdhte(y = df_subgroups1[[outcome]],
                         x = df_subgroups1[[running]],
                         covs.hte = df_subgroups1$leaf)

        rd_hypothesis <- c("`df_subgroups1$leafLarger` - `df_subgroups1$leafLoE` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups1$diff_estimate[1] <- rdhypo$individual$estimate
        subgroups1$diff_pvalue[1] <- rdhypo$individual$p_value

      }

      subgroups1$Nh_left[1] <- sum(rdmodel$Nh[, 1])
      subgroups1$Nh_right[1] <- sum(rdmodel$Nh[, 2])

    }

    ## Save
    subgroups1$var1[1] <- firstvar
    subgroups1$var1_thr[1] <- paste0("<= ", firstvalue)

    # DEPTH LEVEL = 2 **********************************************************
    subgroups2 <- data.frame(pair = rep(NA, 2),
                             var1 = rep(NA, 2),
                             var1_thr = rep(NA, 2),
                             var2 = rep(NA, 2),
                             var2_thr = rep(NA, 2),
                             var3 = rep(NA, 2),
                             var3_thr = rep(NA, 2),
                             diff_estimate = rep(NA, 2),
                             diff_pvalue = rep(NA, 2),
                             Nh_left = rep(NA, 2),
                             Nh_right = rep(NA, 2))

    # Pair 1
    ## Prepare data
    secondvar_index <- p_tree$nodes[[2]]$split_variable
    secondvar <- colnames(matrix_xIN)[secondvar_index]
    secondvalue <- p_tree$nodes[[2]]$split_value

    df_subgroups2 <- df_designmatrix |>
      filter(!!sym(firstvar) <= firstvalue) |>
      mutate(leaf = if_else(!!sym(secondvar) <= secondvalue, "LoE", "Larger"))

    df_subgroups2$leaf <- factor(df_subgroups2$leaf, levels = c("LoE", "Larger"))

    ## Run RD analysis

    ## by chance, the test set could have 0 observations in one of the subgroups2...
    if(sum(df_subgroups2$leaf == "LoE") > 0 & sum(df_subgroups2$leaf == "Larger") > 0){

      # in small subgroups2, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      tss <- tiny_sample_size(df_subgroups2, outcome, running, "leaf")
      if(tss == "Tiny with two bws"){

        rdmodel <- rdhte(y = df_subgroups2[[outcome]],
                         x = df_subgroups2[[running]],
                         covs.hte = df_subgroups2$leaf,
                         bw.joint = T)

        rd_hypothesis <- c("`df_subgroups2$leafLarger` - `df_subgroups2$leafLoE` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups2$diff_estimate[1] <- rdhypo$individual$estimate
        subgroups2$diff_pvalue[1] <- rdhypo$individual$p_value

      }

      if(tss == "Tiny with one bw"){

        message("There was a pair of subgroups with a tiny sample size!")

      }

      else{

        rdmodel <- rdhte(y = df_subgroups2[[outcome]],
                         x = df_subgroups2[[running]],
                         covs.hte = df_subgroups2$leaf)

        rd_hypothesis <- c("`df_subgroups2$leafLarger` - `df_subgroups2$leafLoE` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups2$diff_estimate[1] <- rdhypo$individual$estimate
        subgroups2$diff_pvalue[1] <- rdhypo$individual$p_value

      }

      subgroups2$Nh_left[1] <- sum(rdmodel$Nh[, 1])
      subgroups2$Nh_right[1] <- sum(rdmodel$Nh[, 2])

    }

    ## Save
    subgroups2$pair[1] <- 1
    subgroups2$var1[1] <- firstvar
    subgroups2$var1_thr[1] <- paste0("<= ", firstvalue)
    subgroups2$var2[1] <- secondvar
    subgroups2$var2_thr[1] <- secondvalue

    # Pair 2
    ## Prepare data
    secondvar_index <- p_tree$nodes[[3]]$split_variable
    secondvar <- colnames(matrix_xIN)[secondvar_index]
    secondvalue <- p_tree$nodes[[3]]$split_value

    df_subgroups2 <- df_designmatrix |>
      filter(!!sym(firstvar) > firstvalue) |> # Note the >
      mutate(leaf = if_else(!!sym(secondvar) <= secondvalue, "LoE", "Larger"))

    df_subgroups2$leaf <- factor(df_subgroups2$leaf, levels = c("LoE", "Larger"))

    ## Run RD analysis

    if(sum(df_subgroups2$leaf == "LoE") > 0 & sum(df_subgroups2$leaf == "Larger") > 0){

      # in small subgroups2, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      tss <- tiny_sample_size(df_subgroups2, outcome, running, "leaf")
      if(tss == "Tiny with two bws"){

        rdmodel <- rdhte(y = df_subgroups2[[outcome]],
                         x = df_subgroups2[[running]],
                         covs.hte = df_subgroups2$leaf,
                         bw.joint = T)

        rd_hypothesis <- c("`df_subgroups2$leafLarger` - `df_subgroups2$leafLoE` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups2$diff_estimate[2] <- rdhypo$individual$estimate
        subgroups2$diff_pvalue[2] <- rdhypo$individual$p_value

      }

      if(tss == "Tiny with one bw"){

        message("There was a pair of subgroups with a tiny sample size!")

      }

      else{

        rdmodel <- rdhte(y = df_subgroups2[[outcome]],
                         x = df_subgroups2[[running]],
                         covs.hte = df_subgroups2$leaf)

        rd_hypothesis <- c("`df_subgroups2$leafLarger` - `df_subgroups2$leafLoE` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups2$diff_estimate[2] <- rdhypo$individual$estimate
        subgroups2$diff_pvalue[2] <- rdhypo$individual$p_value

      }

      subgroups2$Nh_left[2] <- sum(rdmodel$Nh[, 1])
      subgroups2$Nh_right[2] <- sum(rdmodel$Nh[, 2])

    }

    ## Save
    subgroups2$pair[2] <- 2
    subgroups2$var1[2] <- firstvar
    subgroups2$var1_thr[2] <- paste0("> ", firstvalue)
    subgroups2$var2[2] <- secondvar
    subgroups2$var2_thr[2] <- secondvalue

    # DEPTH LEVEL = 3 **********************************************************
    subgroups3 <- data.frame(pair = rep(NA, 4),
                            var1 = rep(NA, 4),
                            var1_thr = rep(NA, 4),
                            var2 = rep(NA, 4),
                            var2_thr = rep(NA, 4),
                            var3 = rep(NA, 4),
                            var3_thr = rep(NA, 4),
                            diff_estimate = rep(NA, 4),
                            diff_pvalue = rep(NA, 4),
                            Nh_left = rep(NA, 4),
                            Nh_right = rep(NA, 4))

    # Pair 1
    ## Prepare data
    firstvar_index <- p_tree$nodes[[1]]$split_variable
    firstvar <- colnames(matrix_xIN)[firstvar_index]
    firstvalue <- p_tree$nodes[[1]]$split_value

    secondvar_index <- p_tree$nodes[[2]]$split_variable
    secondvar <- colnames(matrix_xIN)[secondvar_index]
    secondvalue <- p_tree$nodes[[2]]$split_value

    thirdvar_index <- p_tree$nodes[[4]]$split_variable
    thirdvar <- colnames(matrix_xIN)[thirdvar_index]
    thirdvalue <- p_tree$nodes[[4]]$split_value

    df_designmatrix <- as.data.frame(cbind(vector_y, matrix_treatrunning, matrix_x))
    colnames(df_designmatrix)[1] <- outcome
    colnames(df_designmatrix)[3] <- running

    group_left_leaf <- if(p_tree$nodes[8][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups3 <- df_designmatrix |>
      filter(!!sym(firstvar) <= firstvalue) |>
      filter(!!sym(secondvar) <= secondvalue) |>
      mutate(leaf = if_else(!!sym(thirdvar) <= thirdvalue, group_left_leaf, group_right_leaf))

    df_subgroups3$leaf <- factor(df_subgroups3$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    ## by chance, the test set could have 0 observations in one of the subgroups...
    if(sum(df_subgroups3$leaf == "Control") > 0 & sum(df_subgroups3$leaf == "Treat") > 0){

      # in small subgroups, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      tss <- tiny_sample_size(df_subgroups3, outcome, running, "leaf")
      if(tss == "Tiny with two bws"){

        rdmodel <- rdhte(y = df_subgroups3[[outcome]],
                         x = df_subgroups3[[running]],
                         covs.hte = df_subgroups3$leaf,
                         bw.joint = T)

        rd_hypothesis <- c("`df_subgroups3$leafTreat` - `df_subgroups3$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups3$diff_estimate[1] <- rdhypo$individual$estimate
        subgroups3$diff_pvalue[1] <- rdhypo$individual$p_value

      }

      if(tss == "Tiny with one bw"){

        message("There was a pair of subgroups with a tiny sample size!")

      }

      else{

        rdmodel <- rdhte(y = df_subgroups3[[outcome]],
                         x = df_subgroups3[[running]],
                         covs.hte = df_subgroups3$leaf)

        rd_hypothesis <- c("`df_subgroups3$leafTreat` - `df_subgroups3$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups3$diff_estimate[1] <- rdhypo$individual$estimate
        subgroups3$diff_pvalue[1] <- rdhypo$individual$p_value

      }

      subgroups3$Nh_left[1] <- sum(rdmodel$Nh[, 1])
      subgroups3$Nh_right[1] <- sum(rdmodel$Nh[, 2])
    }


    ## Save
    subgroups3$pair[1] <- 1
    subgroups3$var1[1] <- firstvar
    subgroups3$var1_thr[1] <- paste0("<= ", firstvalue)
    subgroups3$var2[1] <- secondvar
    subgroups3$var2_thr[1] <- paste0("<= ", secondvalue)
    subgroups3$var3[1] <- thirdvar
    subgroups3$var3_thr[1] <- thirdvalue


    # Pair 2
    ## Prepare data
    thirdvar_index <- p_tree$nodes[[5]]$split_variable
    thirdvar <- colnames(matrix_xIN)[thirdvar_index]
    thirdvalue <- p_tree$nodes[[5]]$split_value

    group_left_leaf <- if(p_tree$nodes[10][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups3 <- df_designmatrix |>
      filter(!!sym(firstvar) <= firstvalue) |>
      filter(!!sym(secondvar) > secondvalue) |> # Note the >
      mutate(leaf = if_else(!!sym(thirdvar) <= thirdvalue, group_left_leaf, group_right_leaf))

    df_subgroups3$leaf <- factor(df_subgroups3$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    if(sum(df_subgroups3$leaf == "Control") > 0 & sum(df_subgroups3$leaf == "Treat") > 0){

      # in small subgroups3, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      tss <- tiny_sample_size(df_subgroups3, outcome, running, "leaf")
      if(tss == "Tiny with two bws"){

        rdmodel <- rdhte(y = df_subgroups3[[outcome]],
                         x = df_subgroups3[[running]],
                         covs.hte = df_subgroups3$leaf,
                         bw.joint = T)

        rd_hypothesis <- c("`df_subgroups3$leafTreat` - `df_subgroups3$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups3$diff_estimate[2] <- rdhypo$individual$estimate
        subgroups3$diff_pvalue[2] <- rdhypo$individual$p_value

      }

      if(tss == "Tiny with one bw"){

        message("There was a pair of subgroups with a tiny sample size!")

      }

      else{

        rdmodel <- rdhte(y = df_subgroups3[[outcome]],
                         x = df_subgroups3[[running]],
                         covs.hte = df_subgroups3$leaf)

        rd_hypothesis <- c("`df_subgroups3$leafTreat` - `df_subgroups3$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups3$diff_estimate[2] <- rdhypo$individual$estimate
        subgroups3$diff_pvalue[2] <- rdhypo$individual$p_value

      }

      subgroups3$Nh_left[2] <- sum(rdmodel$Nh[, 1])
      subgroups3$Nh_right[2] <- sum(rdmodel$Nh[, 2])
    }


    ## Save
    subgroups3$pair[2] <- 2
    subgroups3$var1[2] <- firstvar
    subgroups3$var1_thr[2] <- paste0("<= ", firstvalue)
    subgroups3$var2[2] <- secondvar
    subgroups3$var2_thr[2] <- paste0("> ", secondvalue)
    subgroups3$var3[2] <- thirdvar
    subgroups3$var3_thr[2] <- thirdvalue

    # Pair 3
    ## Prepare data
    secondvar_index <- p_tree$nodes[[3]]$split_variable
    secondvar <- colnames(matrix_xIN)[secondvar_index]
    secondvalue <- p_tree$nodes[[3]]$split_value

    thirdvar_index <- p_tree$nodes[[6]]$split_variable
    thirdvar <- colnames(matrix_xIN)[thirdvar_index]
    thirdvalue <- p_tree$nodes[[6]]$split_value

    group_left_leaf <- if(p_tree$nodes[12][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups3 <- df_designmatrix |>
      filter(!!sym(firstvar) > firstvalue) |> # Note the >
      filter(!!sym(secondvar) <= secondvalue) |>
      mutate(leaf = if_else(!!sym(thirdvar) <= thirdvalue, group_left_leaf, group_right_leaf))

    df_subgroups3$leaf <- factor(df_subgroups3$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    if(sum(df_subgroups3$leaf == "Control") > 0 & sum(df_subgroups3$leaf == "Treat") > 0){

      # in small subgroups3, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      tss <- tiny_sample_size(df_subgroups3, outcome, running, "leaf")
      if(tss == "Tiny with two bws"){

        rdmodel <- rdhte(y = df_subgroups3[[outcome]],
                         x = df_subgroups3[[running]],
                         covs.hte = df_subgroups3$leaf,
                         bw.joint = T)

        rd_hypothesis <- c("`df_subgroups3$leafTreat` - `df_subgroups3$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups3$diff_estimate[3] <- rdhypo$individual$estimate
        subgroups3$diff_pvalue[3] <- rdhypo$individual$p_value

      }

      if(tss == "Tiny with one bw"){

        message("There was a pair of subgroups with a tiny sample size!")

      }

      else{

        rdmodel <- rdhte(y = df_subgroups3[[outcome]],
                         x = df_subgroups3[[running]],
                         covs.hte = df_subgroups3$leaf)

        rd_hypothesis <- c("`df_subgroups3$leafTreat` - `df_subgroups3$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups3$diff_estimate[3] <- rdhypo$individual$estimate
        subgroups3$diff_pvalue[3] <- rdhypo$individual$p_value

      }

      subgroups3$Nh_left[3] <- sum(rdmodel$Nh[, 1])
      subgroups3$Nh_right[3] <- sum(rdmodel$Nh[, 2])

    }


    ## Save
    subgroups3$pair[3] <- 3
    subgroups3$var1[3] <- firstvar
    subgroups3$var1_thr[3] <- paste0("> ", firstvalue)
    subgroups3$var2[3] <- secondvar
    subgroups3$var2_thr[3] <- paste0("<= ", secondvalue)
    subgroups3$var3[3] <- thirdvar
    subgroups3$var3_thr[3] <- thirdvalue

    # Pair 4
    ## Prepare data
    thirdvar_index <- p_tree$nodes[[7]]$split_variable
    thirdvar <- colnames(matrix_xIN)[thirdvar_index]
    thirdvalue <- p_tree$nodes[[7]]$split_value

    group_left_leaf <- if(p_tree$nodes[14][[1]]$action == 2){"Treat"}else{"Control"}
    group_right_leaf <- if(group_left_leaf == "Treat"){"Control"}else{"Treat"}

    df_subgroups3 <- df_designmatrix |>
      filter(!!sym(firstvar) > firstvalue) |> # Note the >
      filter(!!sym(secondvar) > secondvalue) |> # Note the >
      mutate(leaf = if_else(!!sym(thirdvar) <= thirdvalue, group_left_leaf, group_right_leaf))

    df_subgroups3$leaf <- factor(df_subgroups3$leaf, levels = c("Control", "Treat"))

    ## Run RD analysis

    if(sum(df_subgroups3$leaf == "Control") > 0 & sum(df_subgroups3$leaf == "Treat") > 0){

      # in small subgroups3, the sample size of treated/controls in each leaf could be too small to estimate different bws at each side of the cutoff!
      tss <- tiny_sample_size(df_subgroups3, outcome, running, "leaf")
      if(tss == "Tiny with two bws"){

        rdmodel <- rdhte(y = df_subgroups3[[outcome]],
                         x = df_subgroups3[[running]],
                         covs.hte = df_subgroups3$leaf,
                         bw.joint = T)

        rd_hypothesis <- c("`df_subgroups3$leafTreat` - `df_subgroups3$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups3$diff_estimate[4] <- rdhypo$individual$estimate
        subgroups3$diff_pvalue[4] <- rdhypo$individual$p_value

      }

      if(tss == "Tiny with one bw"){

        message("There was a pair of subgroups with a tiny sample size!")

      }

      else{

        rdmodel <- rdhte(y = df_subgroups3[[outcome]],
                         x = df_subgroups3[[running]],
                         covs.hte = df_subgroups3$leaf)

        rd_hypothesis <- c("`df_subgroups3$leafTreat` - `df_subgroups3$leafControl` = 0")

        rdhypo <- rdhte_lincom(model = rdmodel, linfct = rd_hypothesis)

        # Save
        subgroups3$diff_estimate[4] <- rdhypo$individual$estimate
        subgroups3$diff_pvalue[4] <- rdhypo$individual$p_value

      }

      subgroups3$Nh_left[4] <- sum(rdmodel$Nh[, 1])
      subgroups3$Nh_right[4] <- sum(rdmodel$Nh[, 2])
    }


    ## Save
    subgroups3$pair[4] <- 4
    subgroups3$var1[4] <- firstvar
    subgroups3$var1_thr[4] <- paste0("> ", firstvalue)
    subgroups3$var2[4] <- secondvar
    subgroups3$var2_thr[4] <- paste0("> ", secondvalue)
    subgroups3$var3[4] <- thirdvar
    subgroups3$var3_thr[4] <- thirdvalue


  }

  if(treedepth == 2){
    return(list(gains_assignment,
                rbind(subgroups1, subgroups2),
                p_tree))
  }else{
    return(list(gains_assignment,
                rbind(subgroups1, subgroups2, subgroups3),
                p_tree))
  }

}
