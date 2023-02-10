# crossite_tidy.R

# preprocess_qc function
# Remove variables which meet the following criteria:
#   less than min_values non-NA values
#   TODO do we need min_records?


# Conduct one-way ANOVA for continuous variables; flag on p-values <= min_anova (default 0.001)
# Conduct chi-square test for categorical variables; flag on result <= min_chi (default 0.001)
#

library(dplyr)
library(tidyverse)
library(stddiff)


# add a default for min_values
# mark values with less than x unique values as categorical? save for next release probably
prefilter <- function(tables, min_values, convert_categorical = FALSE) {
  return(lapply(tables, function(t) {
    t <- t |>
      select_if(~ length(unique(.)) > 1) |>
      select_if(~ sum(is.na(.)) > min_values)
  }))
}

#' Site-level Outliers
#'
#' Documentation TODO when inputs is finalized
#'
#' @param tables list of dataframes that contain case report form data
#' @param table_names list of strings representing the names of the tables
#'
#' @return dataframe of outliers
#'
#' @export
multisiteQC <- function(tables,
                      table_names,
                      IDvar,
                      grouping_var,
                      min_p = 0.001,
                      min_p_cat = 0.001,
                      min_std_diff = 0.3,
                      min_std_diff_c = 0.3,
                      variables_to_ignore = NULL,
                      verbose = F,
                      include_all = F) {
  # Build Output Results Dataframe
  all_results <- data.frame(matrix(ncol=11,nrow=0))
  colnames(all_results) <- c("Table","Comparison Variable", "Comparison Value", "Mean", "Sd", "n", "Mean2", "Sd2", "n2", "pval", "stddiff")
  all_results_c <- data.frame(matrix(ncol=5,nrow=0))
  colnames(all_results_c) <- c("Table","Comparison Variable", "Comparison Value", "pval", "stddiff")

  # logging function
  logg <- function(s) {
    if (verbose) {
      print(s)
    }
  }

  selec_valid <- function(t) {
    lapply(t, function(x) {
      x <- x |>
        select_if(~ length(unique(.)) > 1)
    })
  }

  for (i in 1:length(tables)) {
    tryCatch( expr =
                {
    # make sure grouping and idvar arent removed
    variables_to_ignore = variables_to_ignore[variables_to_ignore != IDvar]
    variables_to_ignore = variables_to_ignore[variables_to_ignore != grouping_var]

    df <- as.data.frame(tables[i]) |> select(-any_of(variables_to_ignore))
    # print(colnames(df))
    tablename <- table_names[[i]]
    print(tablename)
    logg(sprintf("Analyzing Table [%s]", tablename))
    print(length(df))
    df[[grouping_var]] <- as.factor(df[[grouping_var]])

    test_variables <- df |> select(-c(all_of(IDvar), all_of(grouping_var)))

    numeric_vars <- test_variables |> select_if(~ is.numeric(.)) |> select_if(~ length(unique(.)) > 2)
    categor_vars <- test_variables |> select_if(~ !is.numeric(.)) |> select_if(~ length(unique(.)) > 1)

    numeric_data <- numeric_vars
    numeric_data[IDvar] <- df[IDvar]
    numeric_data[grouping_var] <- df[grouping_var]

    ## MoTRPAC specific join
    numeric_data <- inner_join(numeric_data, pinfo, by="pid")

    categor_data <- categor_vars
    categor_data[IDvar] <- df[IDvar]
    categor_data[grouping_var] <- df[grouping_var]

    if (ncol(numeric_vars) >= 1) {
      # compute a 1-way ANOVA for each column against grouping_id
      numeric_col_names <- colnames(numeric_vars)
      for (value in numeric_col_names) {
        # anova_data <- na.omit(numeric_data)
        anova_data <- numeric_data
        anova_data$anova_1 <- anova_data[[value]]
        anova_data$anova_2 <- as.factor(anova_data[[grouping_var]])
        anova_res_full <- unlist(summary(aov(anova_1 ~ anova_2 + age + sex + bmi, data=anova_data)))


        print(anova_res_full)
        anova_res <- anova_res_full[9]
        enough_data = sum(!is.na(anova_data$anova_1)) > 20
        if (enough_data && anova_res < min_p) {
          logg(sprintf("Variable [%s] in table [%s] flagged due to 1-way anova result of %.3f", value, tablename, anova_res))

          # get value by site
          for (l in levels(anova_data$anova_2)) {
            # logg(sprintf("Analyzing Site [%s]", l))
            G <- anova_data |> select(anova_1, anova_2) |> mutate(reference = l == anova_2) |> drop_na()
            # check that there are at least 10 measurements for the current site(l)
            n = sum(G$reference)

            # logg(n2)
            if (n > 10) {
              # in-group
              g2 <- G |> subset(reference == T)
              gsum <- mean(g2$anova_1, na.rm = T)
              gstd <- sd(g2$anova_1, na.rm = T)
              # out-group
              g3 <- G |> subset(reference == F)
              gsum_outgroup <- mean(g3$anova_1, na.rm = T)
              gstd_outgroup <- sd(g3$anova_1, na.rm = T)
              res_p = t.test(anova_1 ~ reference, data = G)[3]


              ## RMSE
              mod1 = aov(anova_1 ~ anova_2 + age + sex + bmi, data=anova_data)
              RMSE = sqrt(mod1["Residuals", "Mean Sq"])
              stddiff_num = mod1$coeff["anova_2"]
              adj_stddiff = stddiff_num / RMSE

              # keep ones that are > 0.3

              res_sd <- stddiff.numeric(data = G, gcol = which(names(G) == "reference"), vcol = which(names(G) == "anova_1"))[, "stddiff"]
              out = data.frame("Table" = tablename,
                               "Comparison Variable" = value,
                               "Comparison Value" = toString(l),
                               "Mean" = gsum,
                               "Sd" = gstd,
                               "n" = n,
                               "Mean2" = gsum_outgroup,
                               "Sd2" = gstd_outgroup,
                               "n2" = nrow(g3),
                               "pval" = sprintf("%.3f", res_p),
                               "stddiff" = sprintf("%.3f", res_sd))
              if (!is.na(res_sd)) {
                all_results <- rbind(all_results, out)
              }
            }
          }
        }
      }
   }

    # purposefully excluded
    if (ncol(categor_vars) >= 1 && FALSE) {
      # compute chi square test for each column against grouping_id
      categor_data_names <- colnames(categor_vars)
      for (value in categor_data_names) {
        chisq_data <- na.omit(categor_data)
        chisq_data$chi_1 <- chisq_data[[value]]
        chisq_data$chi_2 <- as.factor(chisq_data[[grouping_var]])
        chi_res <- chisq.test(chisq_data$chi_1, chisq_data$chi_2)$p.value
        # TODO remove TRUE
        if (chi_res < min_p_cat || TRUE) {
          logg(sprintf("Variable [%s] in table [%s] flagged due to chi square test of %.3f", value, tablename, chi_res))

          for (l in levels(chisq_data$chi_2)) {
            G <- chisq_data |> select(chi_1, chi_2) |> mutate(reference = l == chi_2)
            res_p <- chisq.test(G$chi_1, G$reference)$p.value
            res_sd <- stddiff.numeric(data = G, gcol = which(names(G) == "reference"), vcol = which(names(G) == "chi_1"))[, "stddiff"]
            out = data.frame("Table" = tablename,
                             "Comparison Variable" = value,
                             "Comparison Value" = toString(l),
                             "p val" = sprintf("%.3f", res_p),
                             "std diff" = sprintf("%.3f", res_sd))
            all_results_c <- rbind(all_results_c, out)
          }
        }
      }
    }
                },
      error = function(e) {
        print("error caught in trycatch")
        print(e)
      })
    }

  # determine if results are significant or not
  all_results <- all_results |> mutate(meaningful = (pval < min_p & stddiff > min_std_diff))
  # all_results <- subset(all_results, meaningful == TRUE, !is.na(stddiff))


  df <- all_results
  if (!include_all) {
    print('no p vals above ', min_p)
    df <- subset(all_results, meaningful == T)
  }
  if (!verbose) {
    df <- df |> select("Table", "Comparison.Variable", "Comparison.Value", "pval", "stddiff")
    df <- df |> arrange(Table, desc(stddiff))
    rownames(df) <- NULL
  }
  return(df)
}


