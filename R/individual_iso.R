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
individual_multivariate_outliers <- function(crfs, crf_names, id_var, site_var, n_outliers=-1) {
  if (length(crfs) != length(crf_names)) {
    stop("crfs and crf_names are not of equal length")
  }

  res <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(res) <- c("pid", "crf", "score")

  for (i in seq_along(crfs)) {
    crf <- crfs[[i]]
    rename(crf, pid = id_var)
    rename(crf, site = site_var)
    crf_name <- crf_names[i]

    columns_to_test <- crf %>% select(-c(pid, site))
    iso_results <- data.frame(pid = rep(NA, nrow(crf)), crf = rep(crf_name, nrow(crf)), score = rep(NA, nrow(crf)))

    model <- suppressWarnings(isolation.forest(columns_to_test, ndim = ncol(crf)))

    for (row in 1:nrow(crf)) {
      pid <- crf$pid[row]
      cleaned_row <- crf[row,] %>% select(-c(pid, site))
      dist <- predict(model, crf[row, ])
      iso_results[row, ] <- list(pid, crf_name, dist)
    }
    res <- rbind(res, head(iso_results, n_outliers))
  }
  return(dplyr::arrange(res, crf, desc(score)))
}
