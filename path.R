
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

# Functions ---------------------------------------------------------------

rm <- function(paths) {
  ## Remove paths (files/directories) from your system
  ##
  ## @param paths A character vector with the names of the file(s) or directories to be deleted.
  ##
  ## @return A logical value indicating remove status.
  ## @export
  ##
  ## @examples
  ## x <- tempfile()
  ## file.create(x)
  ## yes <- rm_paths(x)
  ## print(yes)
  if (any(sapply(paths, dir.exists)) || any(file.exists(paths))) {
    tryCatch(
      {
        unlink(paths, recursive = TRUE, force = TRUE)
        invisible(TRUE)
      },
      error = function(e) {
        invisible(FALSE)
      }
    )
  } else {
    message("The path does not exist, no need to remove.")
    invisible(FALSE)
  }
}
