# Installs two versions of fastbioclim into separate temporary libraries and runs
# the same issue check script against each.
#
# Usage, from the package root:
#   Rscript dev/issue-checks/compare.R <issue> <base_ref> <pr_ref>
#   Rscript dev/issue-checks/compare.R 22 main batch/1-wrong-results
#   Rscript dev/issue-checks/compare.R all main batch/1-wrong-results
#
# <issue> is an issue number (its script is dev/issue-checks/issue-<issue>.R),
# several numbers separated by commas, or "all". The scripts are taken from the
# current checkout, so they can be run against a ref that does not contain them.
# <base_ref> and <pr_ref> are any git refs (branches, tags or commits).

args <- commandArgs(trailingOnly = TRUE)
options(width = 200)
if (length(args) != 3) {
  stop("Usage: Rscript dev/issue-checks/compare.R <issue|all> <base_ref> <pr_ref>")
}
check_dir <- file.path("dev", "issue-checks")
if (!file.exists("DESCRIPTION") || !dir.exists(check_dir)) {
  stop("Run this script from the package root (the folder with DESCRIPTION).")
}
issues <- args[1]
refs <- c(base = args[2], pr = args[3])

if (identical(issues, "all")) {
  scripts <- list.files(check_dir, pattern = "^issue-[0-9]+\\.R$", full.names = TRUE)
} else {
  nums <- trimws(strsplit(issues, ",")[[1]])
  scripts <- file.path(check_dir, paste0("issue-", nums, ".R"))
  missing <- scripts[!file.exists(scripts)]
  if (length(missing) > 0) stop("No check script: ", paste(missing, collapse = ", "))
}
scripts <- normalizePath(scripts[order(as.integer(gsub("\\D", "", basename(scripts))))])

r_bin <- file.path(R.home("bin"), "R")
rscript_bin <- file.path(R.home("bin"), "Rscript")

install_ref <- function(ref) {
  sha <- system2("git", c("rev-parse", "--verify", "--quiet", shQuote(paste0(ref, "^{commit}"))),
                 stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(sha, "status"))) stop("Unknown git ref: ", ref)
  src <- tempfile("fastbioclim_src_")
  lib <- tempfile("fastbioclim_lib_")
  dir.create(src)
  dir.create(lib)
  tarball <- tempfile(fileext = ".tar")
  st <- system2("git", c("archive", "--format=tar", paste0("--output=", shQuote(tarball)), sha))
  if (st != 0) stop("git archive failed for ", ref)
  utils::untar(tarball, exdir = src)
  cat(sprintf("Installing %s (%s) into %s ...\n", ref, substr(sha, 1, 7), lib))
  out <- system2(r_bin, c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load",
                          paste0("--library=", shQuote(lib)), shQuote(src)),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cat(out, sep = "\n")
    stop("Installation of ", ref, " failed.")
  }
  list(lib = lib, sha = substr(sha, 1, 7))
}

installed <- lapply(refs, install_ref)

run_script <- function(script, lib) {
  env <- paste0("R_LIBS=", paste(c(lib, .libPaths()), collapse = .Platform$path.sep))
  out <- suppressWarnings(system2(rscript_bin, shQuote(script), env = env,
                                  stdout = TRUE, stderr = TRUE))
  res_line <- grep("^RESULT:", out, value = TRUE)
  result <- if (length(res_line) == 0) "ERROR" else trimws(sub("^RESULT:", "", res_line[length(res_line)]))
  list(out = out, result = result)
}

rows <- list()
for (s in scripts) {
  issue <- paste0("#", gsub("\\D", "", basename(s)))
  res <- list()
  for (nm in names(refs)) {
    cat(sprintf("\n===== %s on %s (%s) =====\n", basename(s), refs[[nm]], installed[[nm]]$sha))
    res[[nm]] <- run_script(s, installed[[nm]]$lib)
    cat(res[[nm]]$out, sep = "\n")
  }
  rows[[length(rows) + 1]] <- data.frame(
    issue = issue,
    base = res$base$result,
    pr = res$pr$result,
    fixed = res$base$result != "PASS" && res$pr$result == "PASS",
    stringsAsFactors = FALSE
  )
}
tab <- do.call(rbind, rows)
names(tab) <- c("issue", paste0("base (", refs[["base"]], ")"), paste0("pr (", refs[["pr"]], ")"), "fixed")
cat("\n")
print(tab, row.names = FALSE, right = FALSE)
