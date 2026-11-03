# Check for issue #13: calculate_average()/calculate_sum() (terra) mislabel
# layers when index is not sorted.
#
# Run with: Rscript dev/issue-checks/issue-13.R

suppressPackageStartupMessages({
  library(fastbioclim)
  library(terra)
})
cat("Issue #13: calculate_average()/calculate_sum() (terra) mislabel layers when index is not sorted\n")
cat("fastbioclim version:", as.character(packageVersion("fastbioclim")), "\n")

# One cell, six layers. Groups 3, 1 and 2 appear in that order.
r <- rast(nrows = 1, ncols = 1, nlyrs = 6)
values(r) <- matrix(c(30, 31, 10, 11, 20, 21), 1)
index <- c(3, 3, 1, 1, 2, 2)

avg <- calculate_average(r, index = index, output_names = "x", method = "terra",
                         output_dir = tempfile(), verbose = FALSE)
smm <- calculate_sum(r, index = index, output_names = "s", method = "terra",
                     output_dir = tempfile(), verbose = FALSE)
got_avg <- unlist(avg[1])
got_sum <- unlist(smm[1])

exp_avg <- c(x_01 = 10.5, x_02 = 20.5, x_03 = 30.5)
exp_sum <- c(s_01 = 21, s_02 = 41, s_03 = 61)

cat("EXPECTED: average", paste(names(exp_avg), exp_avg, sep = "=", collapse = " "),
    "| sum", paste(names(exp_sum), exp_sum, sep = "=", collapse = " "), "\n")
cat("GOT:      average", paste(names(got_avg), got_avg, sep = "=", collapse = " "),
    "| sum", paste(names(got_sum), got_sum, sep = "=", collapse = " "), "\n")
ok <- isTRUE(all.equal(got_avg, exp_avg)) && isTRUE(all.equal(got_sum, exp_sum)) &&
  identical(names(avg), names(exp_avg)) && identical(names(smm), names(exp_sum))
cat("RESULT:", if (ok) "PASS" else "FAIL", "\n")
