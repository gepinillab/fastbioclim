fastbioclim 0.4.2.9002
=============
- `derive_bioclim()` with the terra backend passes each static index only to the bio that uses it, so supplying an index meant for another bio no longer stops with `unused argument` (#48).
- `derive_statistics()` with the terra backend accepts static indices (`max_unit`, `min_period`, `max_interactive`, ...). It compares their extent, resolution and CRS with the variable, not their number of layers, and checks that each has one layer (#49).

fastbioclim 0.4.2.9001
=============
- Tiled backends with a `user_region` no longer write values from cells beyond the edge of the cropped grid into other cells, and they now drop the cells outside the region that the terra backend masks (#11).
- `calculate_average()` and `calculate_sum()` with the terra backend now order and name the output layers by sorted `index` values when `index` is not sorted (#13).
- `derive_bioclim()` and `derive_statistics()` apply one NA policy on both backends: an output is NA in a cell when any temporal unit of any variable it is computed from is NA there. The terra backend used to ignore missing units, and the tiled backend only looked at the first layer (#21).
- `derive_statistics(stats = "stdev")` returns the sample standard deviation on the terra backend, as on the tiled backend and in bio04/bio31 (#22).
- The tiled `derive_statistics()` no longer falls back to the full extent when all values inside `user_region` are <= 0 (#23).
- Tiled backends read the selected layers of a multi-band file, e.g. `r[[13:24]]`. `derive_bioclim()` and `derive_statistics()` used to stop with a length error, and the `calculate_*()` functions used the first layers of the file (#27).
- Tiled period-based bios return NA for a static period index beyond the number of periods, as the terra backend does, instead of a period number (#28).
- Tiled outputs keep the CRS of the input, including its EPSG code (#31).
- Add scripts in dev/issue-checks/ that reproduce each of these issues and compare two versions of the package.

fastbioclim 0.4.2.9000
=============
- Require terra >= 1.8-15 because the package calls mem_info(print = FALSE) (#9, closes #15).
- Remove LazyData from DESCRIPTION and document internal function arguments (#9).
- Add the CRAN installation instructions to the README and remove committed .DS_Store files (#9).
- Add tests for terra and tiled backend equivalence, hand-computed bio01-bio19 values and edge cases. Five tests for known bugs are skipped for now (#10).

fastbioclim 0.4.2
=============
- CITATION files and text added.

fastbioclim 0.4.1
=============
- Update testthat after aggregate_ instead of average_

fastbioclim 0.4.0
=============
- New calculate_sum()
- Improved documentation
- Standarized parameters
- Updated bionames()
- Updated README and DESCRIPTION files
- New Logo

fastbioclim 0.3.0
=============
- Changing to qs2 instead of qs package
- Update bionames
- Fix name template
- Fix progressbars in tutorials
- Verbose parameter
- Consistent warning messages
- Unit tests implemented (with the help of Gemini Pro 2.5)

fastbioclim 0.2.1
=============
- Delete unused functions (e.g., mismatch_NA)
- Verbose added.

fastbioclim 0.2.1
=============
- Small changes for making the package ready for CRAN.
- Update example data in README file.

fastbioclim 0.2.0
=============
- New package based on bioclima (https://github.com/gepinillab/bioclima)
- New workflow to handle calculations based on Rfast and exactextract.
- New functions to calculate averages and roll averages.
- Tutorials in Spanish and English.
- Canned data available in egdata.fastbioclim (just in GitHub)

bioclima 0.1.0
=============
- Easier names for internal functions
- Option to add bioclimatic variables bio24-35
- bionames() to get names easily
- Option to add static periods. Useful when you wnat to use same periods in time-series
- New stats_clima() to make your own summaries of variables (e.g., evapotraspiraion, cloudiness)

bioclima 0.0.4
=============
- Fix bug in mismatch_NA().
- Some small changes in documenting the package.

bioclima 0.0.3
=============
- Changing pipe and adding call to terra functions and internal functions. Thanks to the contribution of [@francisvolh](https://github.com/francisvolh)
- Updated node in pkgdown.

bioclima 0.0.2
=============
- Vignette in Spanish to be show in the "Modelamiento de Nichos Ecológicos" facebook group
- Add stdev for population in bio04 and bio15
- No round of variables
- More infomation in README file

bioclima 0.0.1
=============
- First code public