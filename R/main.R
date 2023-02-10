#' @export
QC_Summary <- function(tables,
                       table_names,
                       IDvar,
                       comparison_var,
                       grouping_vars,
                       metadata,
                       varsToIgnore,
                       min_values=10,
                       min_records=10,
                       min_p=0.01,
                       min_std_diff=0.2,
                       verbose=F) {
  crs_res <- crosssite(tables,
                       table_names,
                       IDvar,
                       comparison_var,
                       grouping_vars,
                       metadata,
                       varsToIgnore,
                       min_values=min_values,
                       min_records=min_records,
                       min_p=min_p,
                       min_std_diff=min_std_diff,
                       verbose=verbose)
  # idv_res <- individual_outliers(tables,
                                 # table_names,
                                 # n_outliers=25)
  print("Cross-Site Summary")
  print(crs_res)
  # print("Individual-Level Outliers (Higher Score = Outlier)")
  # print(idv_res)
  # return(list(crs_res,idv_res))
  return(crs_res)
}
