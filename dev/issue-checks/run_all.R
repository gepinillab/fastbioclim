# Runs every dev/issue-checks/issue-NN.R script and prints a PASS/FAIL table.
#
# Usage, from the package root:
#   Rscript dev/issue-checks/run_all.R              # installs this checkout into a temporary library
#   Rscript dev/issue-checks/run_all.R --installed  # uses the fastbioclim already installed
#   Rscript dev/issue-checks/run_all.R --lib=PATH   # uses the fastbioclim installed in PATH
#
# Each script runs in its own R process and must print a line "RESULT: PASS" or
# "RESULT: FAIL". A script that errors or prints neither counts as ERROR.

args <- commandArgs(trailingOnly = TRUE)
options(width = 200)
check_dir <- file.path("dev", "issue-checks")
if (!file.exists("DESCRIPTION") || !dir.exists(check_dir)) {
  stop("Run this script from the package root (the folder with DESCRIPTION).")
}

lib <- NULL
lib_arg <- grep("^--lib=", args, value = TRUE)
if (length(lib_arg) > 0) {
  lib <- normalizePath(sub("^--lib=", "", lib_arg[1]))
} else if (!"--installed" %in% args) {
  lib <- tempfile("fastbioclim_lib_")
  dir.create(lib)
  src <- tempfile("fastbioclim_src_")
  dir.create(src)
  # Copy the package files only, so a dirty dev/ or output folder is not installed.
  pkg_files <- c("DESCRIPTION", "NAMESPACE", "R", "inst", "man")
  pkg_files <- pkg_files[file.exists(pkg_files)]
  file.copy(pkg_files, src, recursive = TRUE)
  cat("Installing this checkout into", lib, "...\n")
  out <- system2(file.path(R.home("bin"), "R"),
                 c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load",
                   paste0("--library=", shQuote(lib)), shQuote(src)),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cat(out, sep = "\n")
    stop("Installation failed.")
  }
}

env <- character()
if (!is.null(lib)) env <- paste0("R_LIBS=", paste(c(lib, .libPaths()), collapse = .Platform$path.sep))
version <- system2(file.path(R.home("bin"), "Rscript"),
                   c("-e", shQuote("cat(as.character(packageVersion('fastbioclim')))")),
                   env = env, stdout = TRUE, stderr = TRUE)
cat("fastbioclim version under test:", version[length(version)], "\n\n")

scripts <- list.files(check_dir, pattern = "^issue-[0-9]+\\.R$", full.names = TRUE)
scripts <- scripts[order(as.integer(gsub("\\D", "", basename(scripts))))]

rows <- lapply(scripts, function(s) {
  t0 <- Sys.time()
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(s),
                                  env = env, stdout = TRUE, stderr = TRUE))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res_line <- grep("^RESULT:", out, value = TRUE)
  result <- if (length(res_line) == 0) "ERROR" else trimws(sub("^RESULT:", "", res_line[length(res_line)]))
  title_line <- grep("^Issue #[0-9]+:", out, value = TRUE)
  title <- if (length(title_line) == 0) "" else sub("^Issue #[0-9]+: *", "", title_line[1])
  if (result == "ERROR") {
    cat("---", basename(s), "output ---\n")
    cat(utils::tail(out, 20), sep = "\n")
  }
  data.frame(issue = paste0("#", gsub("\\D", "", basename(s))),
             result = result,
             seconds = round(secs, 1),
             title = substr(title, 1, 70),
             stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE, right = FALSE)
cat("\n", sum(tab$result == "PASS"), "PASS,", sum(tab$result == "FAIL"), "FAIL,",
    sum(tab$result == "ERROR"), "ERROR\n")
