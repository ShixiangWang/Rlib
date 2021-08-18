
# Utils -------------------------------------------------------------------

.check_install <- function(pkg, bioc = FALSE, ...) {
  if (length(pkg) > 1) sapply(pkg, .check_install, ...)
  install_func <- if (bioc) BiocManager::install else utils::install.packages
  if (bioc) {
    .check_install("BiocManager")
  }
  if (!requireNamespace(pkg)) install_func(pkg, ...)
  message("Required package ", pkg, " has been installed.")
}

# Plot functions ----------------------------------------------------------

plot_roc <- function(
  df, markers, response = "Response", direction = "<",
  legend.position = "right",
  smooth = FALSE,
  auc = TRUE,
  ci = TRUE,
  add_ref_line = TRUE) {
  ## Plot ROC curves
  ##
  ## @param df a `data.frame`.
  ## @param markers column names specifying predictors.
  ##
  ## @examples
  ## library(pROC); data(aSAH)
  ## plot_roc(aSAH, c("s100b", "wfns"), "outcome")
  .check_install(c("dplyr", "ggplot2", "purrr", "pROC"))
  modules::import("dplyr")
  modules::import("ggplot2")
  
  nm <- length(markers)
  if (length(direction) == 1L & nm != 1L) {
    direction <- rep(direction, nm)
  }
  
  ROC_LIST <- purrr::map2(
    markers,
    direction,
    ~ pROC::roc(df[[response]], df[[.x]], direction = .y)
  ) %>% purrr::set_names(markers)
  
  p <- pROC::ggroc(
    ROC_LIST,
    aes = c("color", "linetype"),
    legacy.axes = TRUE
  ) +
    ggplot2::scale_color_discrete("") +
    ggplot2::scale_linetype_discrete("") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion()) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion())
  ggplot2::theme(legend.position = legend.position)
  
  if (add_ref_line) {
    p <- p + ggplot2::geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                          color = "grey", linetype = "dashed")
  }
  list(
    p = p,
    ROC_LIST = ROC_LIST
  )
}
