# Changelog

## fastbioclim 0.4.2.9000

- Require terra \>= 1.8-15 because the package calls mem_info(print =
  FALSE) ([\#9](https://github.com/gepinillab/fastbioclim/issues/9),
  closes [\#15](https://github.com/gepinillab/fastbioclim/issues/15)).
- Remove LazyData from DESCRIPTION and document internal function
  arguments ([\#9](https://github.com/gepinillab/fastbioclim/issues/9)).
- Add the CRAN installation instructions to the README and remove
  committed .DS_Store files
  ([\#9](https://github.com/gepinillab/fastbioclim/issues/9)).
- Add tests for terra and tiled backend equivalence, hand-computed
  bio01-bio19 values and edge cases. Five tests for known bugs are
  skipped for now
  ([\#10](https://github.com/gepinillab/fastbioclim/issues/10)).

## fastbioclim 0.4.2

CRAN release: 2026-05-08

- CITATION files and text added.

## fastbioclim 0.4.1

CRAN release: 2026-03-10

- Update testthat after aggregate\_ instead of average\_

## fastbioclim 0.4.0

- New calculate_sum()
- Improved documentation
- Standarized parameters
- Updated bionames()
- Updated README and DESCRIPTION files
- New Logo

## fastbioclim 0.3.0

CRAN release: 2025-11-26

- Changing to qs2 instead of qs package
- Update bionames
- Fix name template
- Fix progressbars in tutorials
- Verbose parameter
- Consistent warning messages
- Unit tests implemented (with the help of Gemini Pro 2.5)

## fastbioclim 0.2.1

CRAN release: 2025-09-08

- Delete unused functions (e.g., mismatch_NA)
- Verbose added.

## fastbioclim 0.2.1

CRAN release: 2025-09-08

- Small changes for making the package ready for CRAN.
- Update example data in README file.

## fastbioclim 0.2.0

- New package based on bioclima
  (<https://github.com/gepinillab/bioclima>)
- New workflow to handle calculations based on Rfast and exactextract.
- New functions to calculate averages and roll averages.
- Tutorials in Spanish and English.
- Canned data available in egdata.fastbioclim (just in GitHub)
