library(dplyr)
library(isotree)
library(plyr)


#' Individual-level Outliers
#'
#' Finds individual-level outliers using isolation forests. Results are returned on a case-report-form level. Higher 'score' indicates a higher chance the crf is an outlier.
#'
#' @param crfs list of dataframes that contain case report form data
#' @param crf_names list of strings representing the names of the crfs
#' @param n_outliers number of outliers to return; if n_outliers == -1, all data is returned
#'
#' @return List of dataframes
#'
#' @examples
#' crfs <- list(GRIP, CPET)
#' crf_names <- c("grip", "cpet")
#' print(individual_outliers(crfs, crf_names))
#'
#' @export
individual_outliers <- function(crfs, crf_names, n_outliers=-1) {
  if (length(crfs) != length(crf_names)) {
    stop("crfs and crf_names are not of equal length")
  }

  res <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(res) <- c("pid", "crf", "score")

  for (i in seq_along(crfs)) {
    crf <- crfs[[i]]
    crf_name <- crf_names[i]
    iso_results <- data.frame(pid = rep(NA, nrow(crf)), crf = rep(crf_name, nrow(crf)), score = rep(NA, nrow(crf)))

    model <- suppressWarnings(isolation.forest(crf, ndim = ncol(crf)))

    for (row in 1:nrow(crf)) {
      pid <- crf$pid[row]
      dist <- predict(model, crf[row, ])
      iso_results[row, ] <- list(pid, crf_name, dist)
    }
    res <- rbind(res, head(iso_results, n_outliers))
  }
  # if (n_outliers == -1) {
  #   return(dplyr::arrange(res,desc(crf), desc(score)))
  # } else {
  #   return(head(dplyr::arrange(res, desc(score)), n_outliers))
  # }
  return(dplyr::arrange(res, crf, desc(score)))
}
# todo fix printing




#' Finds tail-mean and tail-sd of a column
#' Gives the bounds to check for outliers todo fix docs
tail_stats <- function(vec, sds) {
  lower_bound <- quantile(vec, .1, na.rm=TRUE)
  upper_bound <- quantile(vec, .9, na.rm=TRUE)

  # tail populations
  lower_vals <- vec %>% .[.<lower_bound]
  lower_vals <- lower_vals[!is.na(lower_vals)]
  upper_vals <- vec %>% .[.>upper_bound]
  upper_vals <- upper_vals[!is.na(upper_vals)]

  return(list(tail_low  = mean(lower_vals) - sds*sd(lower_vals),
              tail_high = mean(upper_vals) + sds*sd(upper_vals),
              l_m = mean(lower_vals),
              u_m = mean(upper_vals),
              l_s = sd(lower_vals),
              u_s = sd(upper_vals)))
}

#' Univariate Participant-Level Outlier Detection
#' - values S_1(def: 8) levels from the mean are marked
#' - values S_2(def: ?) tail standard deviations away from respective tail mean are marked
#' - negatives or zeros are marked if less than L_1(def: ?) values are zeros or negative
#'
#' @export
univariate_outliers <- function(crfs, crf_names, S1 = 3, S2 = 3, L1, pass = 1) {

  outliers <- data.frame(matrix(ncol = 10, nrow = 0))
  colnames(outliers) <- c("crf", "studytype", "pid", "var", "value", "outliertype", "overallmean",
                          "sdfromoverallmean", "meanofrelevanttail", "sdfromtailmean")

  outlier_indices <- list()

  for (crf_id in 1:length(crfs)) {
    # remove non-numeric and all-NA columns
    df_raw <- crfs[[crf_id]]
    df <- dplyr::select_if(df_raw, is.numeric)
    df <- df[,colSums(is.na(df))<nrow(df)]

    # find cutoffs after which a value is an outlier
    means <- sapply(df, mean, na.rm=TRUE)
    SDs <- sapply(df, sd, na.rm=TRUE)
    upper_bounds <- means + S1 * SDs
    lower_bounds <- means - S1 * SDs

    tailstats <- sapply(df, tail_stats, sds=S2)

    # list of rows to remove if this is the first pass
    to_remove <- c()

    # print(tailstats)
    for(cid in 1:ncol(df)) {
      tailinfo <- tailstats[,cid]
      if (anyNA(tailinfo)) {
        next
      }
      for (row in 1:nrow(df)) {
        val <- df[row, cid]
        if (is.na(val)) {
          next
        }
        outlier = FALSE
        if (val > upper_bounds[cid]) {
          # outlier (high)
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           1,
                                           means[cid],
                                           abs(means[cid] - val) / SDs[cid],
                                           NA,
                                           NA)
          outlier = TRUE
        }
        else if (lower_bounds[cid] > val) {
          # outlier (low)
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           2,
                                           means[cid],
                                           abs(means[cid] - val) / SDs[cid],
                                           NA,
                                           NA)
          outlier = TRUE
        }
        else if (val < tailinfo$tail_low) {
          # lower than the tail bound
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           3,
                                           NA,
                                           NA,
                                           tailinfo$l_m,
                                           abs(tailinfo$l_m - val) / tailinfo$l_s)
          outlier = TRUE
        }
        else if (val > tailinfo$tail_high) {
          # higher than the tail bound
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           4,
                                           NA,
                                           NA,
                                           tailinfo$u_m,
                                           abs(tailinfo$u_m - val) / tailinfo$u_s)
          outlier = TRUE
        }
        if (outlier) {
          to_remove <- c(to_remove, row)
        }
      }
    }
    # print(to_remove)
    outlier_indices <- c(outlier_indices, crf_id = to_remove)
  }
  if (pass == 1) {
    return(outliers)
  }
  return(list(outliers, outlier_indices))
}

# Univariate outliers with second pass
univariate_outliers_2pass <- function(crfs, crf_names, variables_to_ignore, S1 = 3, S2 = 3, L1, pass = 1) {

  print("Running univariate outliers tuesday")

  outliers <- data.frame(matrix(ncol = 13, nrow = 0))
  colnames(outliers) <- c("crf", "type_placeholder", "pid", "var", "value", "overallmean",
                          "sdfromoverallmean", "meanofrelevanttail", "sdfromtailmean", "outliertype", "siteid", "visitcode", "date")
  outliers$date <- as.character(outliers$date)
  outlier_indices <- list()

  for (crf_id in 1:length(crfs)) {
    tryCatch(expr = {
    print(crf_id)
    df_raw <- crfs[[crf_id]]
    visit_info <- df_raw |> select(visitcode, siteID, d_visit)
    df <- dplyr::select_if(df_raw, is.numeric) |> select(-any_of(variables_to_ignore))
    df <- df[,colSums(is.na(df))<nrow(df)]


    means <- sapply(df, mean, na.rm=TRUE)
    SDs <- sapply(df, sd, na.rm=TRUE)
    if (!is.numeric(SDs)) {
      print("couldnt calculate SD of table")
      next
    }
    upper_bounds <- means + S1 * SDs
    lower_bounds <- means - S1 * SDs

    tailstats <- sapply(df, tail_stats, sds=S2)

    # list of rows to remove if this is the first pass
    to_remove <- c()

    # FIRST PASS
    for(cid in 1:ncol(df)) {
      tailinfo <- tailstats[,cid]
      if (anyNA(tailinfo)) {
        next
      }
      for (row in 1:nrow(df)) {
        val <- df[row, cid]
        visit_info_local <- visit_info[row,]
        siteID <- visit_info_local$siteID
        visitcode <- visit_info_local$visitcode
        date <- visit_info_local$d_visit
        if (is.na(val)) {
          next
        }
        outlier = FALSE
        tailmean <- if (val < means[cid]) tailinfo$l_m else tailinfo$u_m
        tailSDdelta <- if (val < means[cid]) {
            abs(tailinfo$l_m - val) / tailinfo$l_s
          } else {
            abs(tailinfo$u_m - val) / tailinfo$u_s
          }
        if (val > upper_bounds[cid]) {
          # outlier (high)

          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           means[cid],
                                           abs(means[cid] - val) / SDs[cid],
                                           tailmean,
                                           tailSDdelta,
                                           "First pass, too high from overall mean",
                                           siteID,
                                           visitcode,
                                           as.character(date))
          outlier = TRUE
        }
        else if (lower_bounds[cid] > val) {
          # outlier (low)
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           means[cid],
                                           abs(means[cid] - val) / SDs[cid],
                                           tailmean,
                                           tailSDdelta,
                                           "First pass, too low from overall mean",
                                           siteID,
                                           visitcode,
                                           as.character(date))
          outlier = TRUE
        }
        else if (val < tailinfo$tail_low) {
          # lower than the tail bound
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           means[cid],
                                           abs(means[cid] - val) / SDs[cid],
                                           tailinfo$l_m,
                                           tailSDdelta,
                                           "First pass, too low from tail mean",
                                           siteID,
                                           visitcode,
                                           as.character(date))
          outlier = TRUE
        }
        else if (val > tailinfo$tail_high) {
          # higher than the tail bound
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           means[cid],
                                           abs(means[cid] - val) / SDs[cid],
                                           tailinfo$u_m,
                                           tailSDdelta,
                                           "First pass, too high from tail mean",
                                           siteID,
                                           visitcode,
                                           as.character(date))
          outlier = TRUE
        }
        if (outlier) {
          to_remove <- c(to_remove, row)
        }
      }
    }


    # print("slicing")
    if (length(to_remove) != 0) {
      df <- slice(df, -to_remove)
    }

    # print("SECOND PASS")
    ## second pass
    S1 = 5
    S2 = 5
    means <- sapply(df, mean, na.rm=TRUE)
    SDs <- sapply(df, sd, na.rm=TRUE)
    upper_bounds <- means + S1 * SDs
    lower_bounds <- means - S1 * SDs

    tailstats <- sapply(df, tail_stats, sds=S2)

    # print(tailstats)
    for(cid in 1:ncol(df)) {
      tailinfo <- tailstats[,cid]

      if (anyNA(tailinfo)) {
        next
      }
      for (row in 1:nrow(df)) {
        val <- df[row, cid]
        visit_info_local <- visit_info[row,]
        siteID <- visit_info_local$siteID
        visitcode <- visit_info_local$visitcode
        date <- visit_info_local$d_visit
        if (is.na(val)) {
          next
        }
        outlier = FALSE
        if (val > upper_bounds[cid]) {
          # outlier (high)
          tailmean <- if (val < means[cid]) tailinfo$l_m else tailinfo$u_m
          tailSDdelta <- if (val < means[cid]) {
            abs(tailinfo$l_m - val) / tailinfo$l_s
          } else {
            abs(tailinfo$u_m - val) / tailinfo$u_s
          }
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           means[cid],
                                           abs(means[cid] - val) / SDs[cid],
                                           tailmean,
                                           tailSDdelta,
                                           "Second pass, too high from overall mean",
                                           siteID,
                                           visitcode,
                                           as.character(date))
        }
        else if (lower_bounds[cid] > val) {
          # outlier (low)
          tailmean <- if (val < means[cid]) tailinfo$l_m else tailinfo$u_m
          tailSD <- if (val < means[cid]) tailinfo$l_s else tailinfo$u_s
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           means[cid],
                                           abs(means[cid] - val) / SDs[cid],
                                           tailmean,
                                           tailSDdelta,
                                           "Second pass, too low from overall mean",
                                           siteID,
                                           visitcode,
                                           as.character(date))
        }
        else if (val < tailinfo$tail_low) {
          # lower than the tail bound
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           means[cid],
                                           abs(means[cid] - val) / SDs[cid],
                                           tailinfo$l_m,
                                           tailSDdelta,
                                           "Second pass, too low from tail mean",
                                           siteID,
                                           visitcode,
                                           as.character(date))
        }
        else if (val > tailinfo$tail_high) {
          # higher than the tail bound
          outliers[nrow(outliers)+1,] <- c(crf_names[crf_id],
                                           "SED",
                                           df_raw[row, "pid"],
                                           colnames(df)[cid],
                                           val,
                                           means[cid],
                                           abs(means[cid] - val) / SDs[cid],
                                           tailinfo$u_m,
                                           tailSDdelta,
                                           "Second pass, too high from tail mean",
                                           siteID,
                                           visitcode,
                                           as.character(date))
        }
      }
    }
    },
  error = function(e) {
    print(e)
  })}
  return(outliers)
}



subset4 <- function(dfs) {
  for (dfid in 1:length(dfs)) {
    df <- dfs[[dfid]]
    df |> slice(-c(1,2))
  }
}



# funct that sources all the target files, result is a target, returns null
