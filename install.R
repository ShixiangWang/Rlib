.check_install <- function(pkg) {
  if (!requireNamespace(pkg)) install.packages(pkg)
  message("Required package ", pkg, " has been installed.")
}

.verbose_git <- function() message("See ?remotes::install_git fore more options.")
.verbose_github <- function() message("See ?remotes::install_github fore more options.")
.verbose_bioc <- function() message("See ?BiocManager::install for more options.")
.verbose_clone <- function() message("See ?git2r::clone for more options.")

#' Install R packages
#'
#' @param pkg A package name for Git(hub) package (e.g. `ShixiangWang/ezcox`,
#' `ShixiangWang/wfun@main`) or
#' a list of packages (e.g. `c("dplyr", "maftools")`) for CRAN/BioC packages.
#' Or a file path (a directory or zip file).
#' @param gitee If `TRUE`, install package from Gitee.
#' @param ... Other arguments passing to [remotes::install_git],
#' [BiocManager::install] or [remotes::install_github] based on input.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' \dontrun{
#' install("ShixiangWang/ezcox")
#' install("ShixiangWang/tinyscholar", gitee = TRUE)
#' install(c("ggplot2", "Biobase"))
#' }
install <- function(pkg, gitee = FALSE, ...) {
  if (file.exists(pkg) || dir.exists(pkg)) {
    if (dir.exists(pkg)) {
      message("Installing package from a directory...")
      utils::install.packages(pkg, repos = NULL)
    } else {
      message("Installing package from a zip file")
      .install_zip(pkg)
    }
    return(invisible(0))
  }

  .check_install("remotes")

  stopifnot(is.logical(gitee))
  if (gitee) {
    if (!grepl("/", pkg)) stop("Invalid package name, should be name/repo format.")
    .verbose_git()
    remotes::install_git(paste0("https://gitee.com/", pkg), ...)
  } else {
    if (any(grepl(":", pkg))) {
      .verbose_git()
      remotes::install_git(pkg, ...)
    } else {
      if (any(grepl("/", pkg))) {
        tryCatch(
          {
            message("Installing from GitHub mirrors...")
            if (grepl("@", pkg)) {
              pkg2 <- unlist(strsplit(pkg, "@"))
              pkg <- pkg2[1]
              ref <- pkg2[2]
            } else ref = "HEAD"
            .install_zip_gh(pkg, ref = ref)
          },
          error = function(e) {
            message("Install error when using GitHub mirror, roll back to official GitHub.")
            .verbose_github()
            remotes::install_github(pkg, ...)
          }
        )
      } else {
        .verbose_bioc()
        .check_install("BiocManager")
        BiocManager::install(pkg, ...)
      }
    }
  }
  
  return(invisible(0))
}

#' Clone a Git repository
#'
#' @inheritParams git2r::clone
#' @param gitee If `TRUE`, clone repository from Gitee.
#' @param reset_remote if `TRUE`, reset GitHub repository remote url.
#' @param ... Other arguments passing to [git2r::clone]
#'
#' @return Nothing
#' @export
#'
#' @examples
#' \donttest{
#' x <- file.path(tempdir(), "ezcox")
#' if (dir.exists(x)) rm_paths(x)
#' clone("ShixiangWang/ezcox", x, reset_remote = TRUE)
#'
#' y <- file.path(tempdir(), "tinyscholar")
#' if (dir.exists(y)) rm_paths(y)
#' clone("ShixiangWang/tinyscholar", y, gitee = TRUE)
#' }
clone <- function(url, local_path, gitee = FALSE, reset_remote = FALSE, ...) {
  stopifnot(length(url) == 1L)
  .check_install("git2r")

  local_path <- path.expand(local_path)
  if (!grepl(":", url)) {
    if (gitee) {
      url <- paste0("https://gitee.com/", url)
    } else {
      message("Treat input as a GitHub repo.")
      url <- paste0("https://hub.fastgit.org/", url)
    }
  }
  url <- sub("github.com", "hub.fastgit.org", url, fixed = TRUE)
  git2r::clone(url, local_path, bare = FALSE, ...)
  if (reset_remote) {
    if (grepl("fastgit", url)) {
      url <- sub("hub.fastgit.org", "github.com", url, fixed = TRUE)
      message("Reset remote url to ", url)
      git2r::remote_set_url(local_path,
        name = git2r::remotes(local_path),
        url
      )
    } else {
      message("No need to reset.")
    }
  }
}

#' Download GitHub/Gitee repo release or archive file
#'
#' @param repo A GitHub/Gitee repo in the format of `username/repo`, e.g. `ShixiangWang/tinyscholar`.
#' @param destdir A target path to save the file.
#' @param release Set to the version number (e.g. `v1.0`) for downloading release file.
#' @param gitee If `TRUE`, download from Gitee repo instead of GitHub repo.
#' @param ... Other arguments passing to [utils::download.file].
#'
#' @return Nothing
#' @export
#'
#' @examples
#' \donttest{
#' x <- tempdir()
#' download("ShixiangWang/tinyscholar", destdir = x)
#' dir(x)
#' }
download <- function(repo, destdir, release = NULL, gitee = FALSE, ...) {
  stopifnot(
    length(repo) == 1L, length(destdir) == 1L,
    is.character(repo), is.character(destdir)
  )

  if (!grepl("/", repo)) stop("Repo should in format of username/repo")

  if (!dir.exists(destdir)) dir.create(destdir, recursive = TRUE)
  if (is.null(release)) {
    message("Downloading repo archive...")
    url <- sprintf(
      "https://%s/%s/archive/master.zip",
      if (gitee) "gitee.com" else "download.fastgit.org",
      repo
    )
  } else {
    url <- sprintf(
      "https://%s/%s/releases/download/%s/%s.tar.gz",
      if (gitee) "gitee.com" else "download.fastgit.org",
      repo, release
    )
  }
  utils::download.file(url, destfile = file.path(destdir, basename(url)), ...)
}


# Copy from yu utils ------------------------------------------------------
# source: https://github.com/YuLab-SMU/yulab.utils/blob/master/R/.install_zip.R

##' install github package
##'
##' it download the zip file first and use `.install_zip` to install it
##' @title .install_zip_gh
##' @param repo github repo
##' @param ref github branch, default is master
##' @param args argument to build package
##' @return NULL
##' @export
##' @author Guangchuang Yu
.install_zip_gh <- function(repo, ref = "master", args = "--no-build-vignettes") {
  url <- paste0("https://codeload.github.com/", repo, "/zip/", ref)
  f <- tempfile(fileext = ".zip")
  method <- "auto"
  if (.Platform$OS.type == "windows") method <- "curl"
  utils::download.file(url, destfile = f, method = method)
  .install_zip(f, args = args)
}

##' install R package from zip file of source codes
##'
##'
##' @title .install_zip
##' @param file zip file
##' @param args argument to build package
##' @return NULL
##' @export
##' @author Guangchuang Yu
.install_zip <- function(file, args = "--no-build-vignettes") {
  .check_install("pkgbuild")

  dir <- tempfile()
  utils::unzip(file, exdir = dir)
  fs <- list.files(path = dir, full.names = T)
  if (length(fs) == 1 && dir.exists(fs)) {
    dir <- fs
  }

  pkg <- pkgbuild::build(dir, args = args)
  utils::install.packages(pkg, repos = NULL)
}


# Save the script to local ------------------------------------------------

save <- function() {
  message("Downloading script...")
  url <- "https://biosisyphus.github.io/Rlib/install.R"
  if (!dir.exists(path.expand("~/.R"))) dir.create(path.expand("~/.R"))
  tryCatch(
    utils::download.file(url, destfile = file.path(path.expand("~/.R"), "install.R")),
    error = function(e) {
      message("Try downloading from gitee...")
      url <- "https://gitee.com/BioSisyphus/Rlib/raw/master/install.R"
      utils::download.file(url, destfile = file.path(path.expand("~/.R"), "install.R"))
    }
  )
  if (!file.exists(path.expand("~/.Rprofile"))) file.create(path.expand("~/.Rprofile"))
  
  message("Writing to config file ~/.Rprofile...")
  write(paste0("\nsource(\"", file.path(path.expand("~/.R"), "install.R\""), ", local = TRUE)"), path.expand("~/.Rprofile"), append = TRUE)
  message("Done. Please restart your R session and see if you can directly use `install()` function.")
}
