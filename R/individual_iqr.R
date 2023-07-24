# IQR-based outlier detection

#' Individual-level Outliers
#'
#' Finds individual-level outliers using interquartile ranges
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
individual_univariate_outliers <- function(crfs, crf_names, n_outliers=-1) {
  if (length(crfs) != length(crf_names)) {
    stop("crfs and crf_names are not of equal length")
  }

  res <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(res) <- c("pid", "crf", "score")


  not_all_na <- function(x) any(!is.na(x))

  outliers <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(outliers) <- c("variable","crf","site","value")
  for(crf_id in 1:length(crfs)) {
    crf <- crfs[[crf_id]]
    sites <- crf$site
    crf_name <- crf_names[[crf_id]]
    crf <- crf |>
      select(where(not_all_na)) |>
      select(where(is.numeric))
    crf$siteID <- sites

    if(nrow(crf) < 10) next
    if(ncol(crf) < 1) next

    for(cid in 1:ncol(crf)) {
      col <- crf[[cid]]
      variable <- colnames(crf)[[cid]]
      fstQ <- summary(col, na.rm=TRUE)[[2]]
      thdQ <- summary(col,na.rm=TRUE)[[5]]
      iqr <- IQR(col,na.rm=TRUE)

      high_cutoff <- 1.5*iqr + thdQ
      low_cutoff  <- fstQ - 1.5*iqr

      for(xid in 1:nrow(crf)) {

        tryCatch(expr = {
          val <- col[[xid]]
          site <- crf[xid, ]$site

          if (!is.na(val) && (val > high_cutoff || low_cutoff > val)) {
            outliers[nrow(outliers)+1,] <- c(variable, crf_name, site, val)
          }
        },
        error=function(e) {
          print(e)
          return(e)
        }
        )
      }
    }
  }
  return(dplyr::arrange(res, crf, desc(score)))
}
