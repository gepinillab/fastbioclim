# fastbioclim code review (`R/`, v0.4.2)

Scope: all files in `R/`. No package code was changed.

## How the findings were checked

- **Environment:** R 4.3.3 on Ubuntu 24.04, with terra 1.7-65, sf, exactextractr 0.10.0, future 1.33.1, future.apply 1.11.1, purrr, glue and progressr.
- **Missing dependencies:** the sandbox's network policy blocks CRAN and its mirrors, and Ubuntu doesn't package **`Rfast` or `qs2`**. The package therefore couldn't be installed. Instead, every file in `R/` was loaded with `sys.source()` and the real functions were run.
  - Everything up to the first `Rfast`/`qs2` call ran for real. That covers the whole terra backend, the wrappers, the `auto` logic, input validation, and the tiled setup (geometry, region, tiling and cell translation).
  - Where a tiled result depends on code after that point, the finding is marked **(by reading)**. In those cases the R lines involved were run verbatim where possible.
- **The examples below** use `library(fastbioclim)` / `fastbioclim:::` so they can be pasted into a normal install. Each one was run as shown, through the sourced code.
- **Severity labels:**
  - **bug** means wrong results or a crash on valid input.
  - **risk** means fragile behaviour likely to cause problems.
  - **improvement** means performance or maintainability.

---

## Summary

The terra backend's formulas are correct, including quarters that wrap around the year end. Every bio01–bio19 value matched an independent dismo-style reference on a southern-hemisphere test cell (finding 0). Most serious problems are in input handling, in the tiled backend's bookkeeping, and in the two backends disagreeing with each other.

| # | Severity | Short description | Confirmed by |
|---|---|---|---|
| 1 | bug | Tiled + `user_region` writes pixels from outside the region into wrong cells | run |
| 2 | bug | `calculate_average()`/`calculate_sum()` tiled without `user_region` always fails (`noMap` undefined) | run |
| 3 | bug | `calculate_average()`/`calculate_sum()` (terra) mislabel layers when `index` is not sorted | run |
| 4 | bug | `bioclim_terra()` picks up global variables named like static indices (`exists()`) | run |
| 5 | bug | `method = "auto"` (the default everywhere) fails on terra versions allowed by DESCRIPTION | run |
| 6 | bug | Static period indices (e.g. `wettest_period`) can't replace the data they stand in for, in either backend | run |
| 7 | bug | bio26/27/34/35 crash in terra unless another tavg-based bio is also requested | run |
| 8 | bug | terra path rejects `tavg`-only input for bio01/04/10/11… | run |
| 9 | bug | Tiled bio05/bio06 require both `tmin` and `tmax` | run |
| 10 | bug | Exported `check_rasters()`: geometry check fails for ≥3 inputs; NA check always errors | run |
| 11 | bug | NA handling: terra silently sums the months it has; tiled returns NA; tiled NA mask comes from one layer | run (terra) / reading (tiled) |
| 12 | bug | `stats = "stdev"`: population SD (terra) vs sample SD (tiled) | run (terra) / reading (tiled) |
| 13 | bug | `stats_fast()`: a region whose first-layer values are all ≤ 0 is treated as "empty" | run |
| 14 | bug | `user_region` in another CRS: terra path and `auto` fail; tiled reprojects | run |
| 15 | bug | `derive_statistics()`: invalid `inter_stats` accepted silently (copy-paste typo) | run |
| 16 | bug | `period_length` not validated in terra path (13 → double-counted month; 2.5 → wrong mean) | run |
| 17 | bug | Tiled backends read the wrong layers for layer-subset inputs and reject multi-band files | partly run |
| 18 | bug | Out-of-range static period index returns the `min_idx`/`max_idx` column as data (bios_fast) | reading |
| 19 | bug | `sf::sf_use_s2(FALSE)` changed globally and never restored | reading |
| 20 | bug | `calculate_roll()` tiled finds no outputs when `name_template` doesn't start with `{prefix}` | reading |
| 21 | risk | Tiled outputs lose the CRS identity (EPSG:3035/2193 not `same.crs` as input) | run |
| 22 | risk | Intermediate `*_qs_*` directories left in `output_dir` after any failure | run |
| 23 | risk | Overwrite handling: pre-check misses `cv_cli`; `write_layers()` ignores `overwrite`; default `tempdir()` collides | run |
| 24 | risk | Weak validation and error messages (`bios` 36, duplicate `bios`, `call. = FALSE` printed as "FALSE") | run |
| A1–A4 | risk/improvement | `auto` method: estimate, preconditions, CRS | run/reading |
| P1–P5 | improvement | Tiled performance and parallel scaling | measured/reading |
| M1–M6 | improvement | Maintainability and API consistency | reading |

---

## 0. Formula check (no defect)

Test cell: monthly tavg `25 24 20 15 10 6 5 7 11 16 21 26` (warmest quarter is Dec–Jan–Feb). Monthly prcp `90 40 30 20 10 5 5 10 20 50 100 120` (wettest quarter is Nov–Dec–Jan). The terra backend (`circular = TRUE`) matched an independent implementation for bio01, 04, 08–19:

```
            bio01    bio04 bio08 bio09 bio10 bio11 bio12 bio13 bio14   bio15 bio16 bio17 bio18 bio19
fastbioclim  15.5 766.9301    24     7    25     6   500   120     5 94.0602   310    20   250    20
reference    15.5 766.9301    24     7    25     6   500   120     5 94.0602   310    20   250    20
```

The tiled formulas are algebraically the same, by reading `bios_fast.R`:
- bio15 = `sd(p) / (1 + mean(p)) * 100`, which equals `sd(p+1)/mean(p+1) * 100`.
- bio04 uses the sample SD.
- `compute_periods()` handles wrap-around the same way as `get_window()`.

Cells that are all NA give NA in every bio. The remaining formula-level differences between backends are in finding 11 (NA), finding 12 (stdev) and bio03 when bio07 = 0: terra gives NaN/Inf, while `bio03_fast` returns 0 (`bios_fast.R`, `bio03_fast`).

---

## Confirmed bugs

### 1. [bug] Tiled + `user_region`: pixels outside the region are written into wrong cells
**Where:** `R/helper_functions.R:284` (`define_translate`), used by `bioclim_fast.R:313`, `stats_fast.R:262`, `aggregate_fast.R:126`, `roll_fast.R:106`. Assembled by `write_layers.R:126-127`.

**What:**
1. The target grid is `terra::crop(ref, region, mask = TRUE)`. Its default `snap = "near"` can drop a cell that the polygon only partly covers.
2. `exact_extract()` still returns that cell, because it returns every cell with `coverage_fraction > 0`.
3. `define_translate()` only checks `row_tgt >= 1 & col_tgt >= 1`. It never checks `col_tgt <= ncol_tgt` or `row_tgt <= nrow_tgt`.
4. So a right-edge cell with `col_tgt = ncol_tgt + 1` wraps to column 1 of the next row. A bottom-edge cell lands past `ncell`, and `rvals[cellID] <- …` silently grows the vector.

The wrapped value survives only if no correct value overwrites it later. That happens when the true target cell is outside the polygon (non-rectangular regions), or when the tile holding the wrapped cell is assembled after the tile holding the true cell. Tile files are read in `list.files()` order: `_1, _10, _11, _2…`.

The tiled path also keeps edge cells whose centre lies outside the polygon, while the terra path masks by cell centre. The two backends' outputs therefore differ along every region boundary even without the wrap.

**Evidence:** real terra, sf and exactextractr calls, plus the package's `define_translate`. The last step reproduces `write_layers`' assignment.
```r
library(terra); library(sf)
ref <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = "EPSG:4326")
values(ref) <- 1:100                                   # value = source cell id
poly <- st_as_sf(st_union(                             # L-shaped region
  st_as_sfc(st_bbox(c(xmin = 2.3, ymin = 4, xmax = 6.3, ymax = 6.3), crs = st_crs(4326))),
  st_as_sfc(st_bbox(c(xmin = 4,   ymin = 2, xmax = 6.3, ymax = 4),   crs = st_crs(4326)))))
tgt <- crop(ref, vect(poly), mask = TRUE)              # as in bioclim_fast.R:257
tr  <- fastbioclim:::define_translate(ncol(ref), ncol(tgt),
         rowFromY(ref, ymax(tgt) - res(tgt)[2]/2) - 1L, colFromX(ref, xmin(tgt) + res(tgt)[1]/2) - 1L)
ex <- exactextractr::exact_extract(ref, poly, include_cell = TRUE)[[1]]
ex$tgt <- tr(ex$cell); ex <- ex[!is.na(ex$tgt), ]
rvals <- rep(NA_real_, ncell(tgt)); rvals[ex$tgt] <- ex$value          # write_layers.R:126-127
cbind(terra = values(tgt)[, 1], tiled = rvals[1:ncell(tgt)])[c(9, 13), ]
#>      terra tiled
#> [1,]    NA    57     <- value of source cell 57 (x = 6.5, outside the region's column range)
#> [2,]    NA    67
```
7 of the 19 extracted cells were translated to the wrong target cell. Three of them were cells with `coverage_fraction` around 1e-16, which are union slivers.

**Fix:**
- In `define_translate()`, take `nrow_tgt` as well and use `valid <- row_tgt >= 1L & col_tgt >= 1L & row_tgt <= nrow_tgt & col_tgt <= ncol_tgt`.
- For parity with the terra path, keep only cells whose centre is inside the region. One way is to build the target mask once, with `crop(..., mask = TRUE)` plus `cells()`, and drop translated cells that are NA in it. Another is a small `coverage_fraction` threshold, but the mask is exact.
- Add a parity test: terra vs tiled on a small raster with a non-rectangular `user_region`.

### 2. [bug] `calculate_average()` / `calculate_sum()` with `method = "tiled"` and no `user_region` always fail
**Where:** `R/aggregate_fast.R:184-191`

**What:** `noMap` is only created inside `if (!is.null(user_region))`, but line 191 always uses `df[nonaID[noMap], …]`. Every tile throws `object 'noMap' not found`, and `future_lapply` passes the error on.

Large rasters with no region are exactly the case `auto` sends to this path. The existing tests mock `aggregate_fast()`, so this branch never runs in CI.

**Evidence:** lines 180–191 run verbatim, with `user_region <- NULL`:
```
Error in g() : object 'noMap' not found
```
**Fix:** set `noMap <- rep(TRUE, length(cell_ids_source))` before the `if`, or index with `nonaID` in the `else` branch as `roll_fast.R:158-168` does.

### 3. [bug] `calculate_average()` / `calculate_sum()` (terra path) attach the wrong names when `index` is not sorted
**Where:** `R/aggregate_terra.R:10-11`, with the names built in `calculate_average.R` / `calculate_sum.R` from `sort(unique(index))`.

**What:** `terra::tapp()` returns groups in order of first appearance. The code then assigns names generated from the *sorted* unique values. Custom seasons or water years (e.g. index `c(10,10,11,11,12,12,1,1,…)`) end up silently mislabelled, and the files are named wrongly too. The tiled path uses sorted order, so the two backends disagree.

**Evidence:**
```r
r <- terra::rast(nrows = 1, ncols = 1, nlyrs = 6); terra::values(r) <- matrix(c(30, 31, 10, 11, 20, 21), 1)
a <- calculate_average(r, index = c(3, 3, 1, 1, 2, 2), output_names = "x", method = "terra",
                       output_dir = tempfile(), verbose = FALSE)
unlist(a[1])
#>  x_01 x_02 x_03
#>  30.5 10.5 20.5       # expected 10.5 20.5 30.5
```
**Fix:** reorder after `tapp`, e.g. `agg_rast <- agg_rast[[order(unique(index))]]`, or pass `index = factor(index, levels = sort(unique(index)))`. Better still, derive the names from the result's own group labels.

### 4. [bug] `bioclim_terra()` uses objects from the user's workspace as static indices
**Where:** `R/bioclim_terra.R:62, 66, 80, 85, 163, 166, 179, 182, 219, 222` (`exists("…")`)

**What:** `exists("warmest_period")` searches the function's enclosing environments. That chain reaches the global environment and every attached package. If the user has any object called `warmest_period`, `coldest_period`, `wettest_period`, `driest_period`, `high_mois_period` or `low_mois_period`, it is silently used instead of the computed index.

**Evidence:**
```r
ok <- derive_bioclim(bios = 10, tmin = tmin, tmax = tmax, method = "terra", output_dir = tempfile(), verbose = FALSE)
warmest_period <- tmin[[1]] * 0 + 1           # unrelated object in the global env
bad <- derive_bioclim(bios = 10, tmin = tmin, tmax = tmax, method = "terra", output_dir = tempfile(), verbose = FALSE)
c(ok[1][[1]], bad[1][[1]])
#> 10.684 1.396     (bio10 silently becomes the mean of the *first* quarter)
```
**Fix:** declare every static index as a formal argument defaulting to `NULL` and test with `is.null()`. Alternatively, keep the `dot_args` list and test `!is.null(dot_args$warmest_period)`. Never use `exists()` without `inherits = FALSE`.

### 5. [bug] Default `method = "auto"` fails on terra versions that DESCRIPTION allows
**Where:**
- `derive_bioclim.R:242-245, 257-260`
- `derive_statistics.R:172, 187`
- `calculate_average.R:112-113`
- `calculate_sum.R:111-112`
- `calculate_roll.R:157-158`
- `DESCRIPTION` (`terra (>= 1.7-0)`)

**What:** the code calls `terra::mem_info(x, n, print = FALSE)` and reads `["fits_mem"]`.
- In terra 1.7-65, the version Ubuntu 24.04 ships, `mem_info()` has signature `(x, n = 1)`.
- Its return value has no `fits_mem` name: the 5th element is unnamed.
- So all five exported functions fail with their default arguments.

**Evidence:**
```
> args(terra::mem_info)        # terra 1.7-65
function (x, n = 1)
> derive_bioclim(bios = 1, tmin = x, tmax = x, output_dir = tempfile())
Error in terra::mem_info(full_raster_stack[[1]], ... : unused argument (print = FALSE)
> mem_info(r, 12)              # return value
   needed  available   memfrac  chunksize   <NA>
```
**Fix:**
- Raise the requirement to the first terra release whose `mem_info()` has `print` and a named `fits_mem` (check terra's NEWS).
- Or don't rely on it: compute `needed <- ncell * nlyr_total * 8` and compare it with `terra::free_RAM() * 1024 * terra::terraOptions(print = FALSE)$memfrac`.
- Add a CI job on the oldest supported terra.

### 6. [bug] Static period indices can't replace the data they stand in for, in either backend
**Where:**
- terra: `bioclim_terra.R:62-69` drops prcp from `req_prcp` when `wettest_period` is supplied, but `:161-162` still calls `get_window(prcp, …)`.
- tiled: `bioclim_fast.R:65-79` (`needs`), `:115`, `:573-579` always require prcp, or tavg, for bio08/09/18/19/24–27/34/35.

**What:** the "Static Indices" section of `?derive_bioclim` presents indices as overriding the dynamic calculation. With `wettest_period` supplied and no `prcp`, the terra path crashes with an unrelated S4 error, and the tiled path stops with `prcp_path required`. The same holds for `warmest_period` / `coldest_period` with bio26/27/34/35, which still need tavg.

**Evidence:**
```r
wp <- terra::writeRaster(tmin[[1]] * 0 + 2, tempfile(fileext = ".tif"))
derive_bioclim(bios = 8, tmin = tmin, tmax = tmax, wettest_period = wp, method = "terra", output_dir = tempfile())
#> Error ... unable to find an inherited method for function 'nlyr' for signature '"NULL"'
fastbioclim:::bioclim_fast(bios = 8, n_units = 12, tmin_path = terra::sources(tmin),
  tmax_path = terra::sources(tmax), wettest_period_path = terra::sources(wp), output_dir = tempfile())
#> Error in bioclim_fast(...) : prcp_path required.
```
**Fix:** compute `wet` / `tmp` only when a bio needs an index that wasn't supplied, or needs the period values of that variable. Drive both backends from one dependency table (see M3).

### 7. [bug] bio26/27/34/35 crash on the terra path unless another tavg-based bio is also requested
**Where:** `R/bioclim_terra.R:124` (`if (any(c(1, 4, 8, 9, 10, 11, 18, 19) %in% bios))`)

**What:** `tavg` is only derived from tmin/tmax for the listed bios. bio26, 27, 34 and 35 also need `tmp` (line 178) to find the warmest/coldest period, so `get_window(NULL, …)` is called.

**Evidence:**
```r
derive_bioclim(bios = 26, tmin = tmin, tmax = tmax, srad = srad, method = "terra", output_dir = tempfile())
#> Error ... unable to find an inherited method for function 'nlyr' for signature '"NULL"'
```
`bios = c(10, 26)` works, which hides the problem in most tests.

**Fix:** add 26, 27, 34 and 35 to the list on line 124, or compute tavg whenever `tmp` is needed (line 177).

### 8. [bug] terra path rejects `tavg`-only input for bios that need only tavg
**Where:** `R/bioclim_terra.R:90-100`

**What:** `if (is.null(tmin) & is.null(tmax)) stop("… require(s) tmin and tmax")` runs before checking whether `tavg` was supplied. The tiled path accepts it (`req_tavg_load`).

**Evidence:**
```r
derive_bioclim(bios = 1, tavg = tavg, method = "terra", output_dir = tempfile())
#> Error ... Bio01 require(s) tmin and tmax.
```
**Fix:** only require tmin/tmax for bio02/03/05/06/07, and for tavg-based bios when `tavg` is `NULL`.

### 9. [bug] Tiled bio05/bio06 require both `tmin` and `tmax`
**Where:** `R/bioclim_fast.R:67` (`needs$tavg` contains 5 and 6) and `:557-562`

**What:** bio05 is the max of tmax and bio06 is the min of tmin, as the package's own documentation says. The tiled path still treats them as tavg-dependent, so `bios = 5` with only tmax fails. The terra path works.

**Evidence:**
```r
fastbioclim:::bioclim_fast(bios = 5, n_units = 12, tmax_path = terra::sources(tmax), output_dir = tempfile())
#> Error in bioclim_fast(...) : tmin_path required.
```
**Fix:** remove 5 and 6 from `needs$tavg`, and remove the `!tavg_available` condition on lines 558-562.

### 10. [bug] Exported `check_rasters()` is broken in both modes
**Where:** `R/check_rasters.R:46`, `:61-66`

**What:**
- `purrr::reduce(raster_list, terra::compareGeom)` returns `TRUE` after the first pair. The next step then calls `compareGeom(TRUE, r3)`. Any three or more inputs are therefore always reported as "FAILURE: Geometries are inconsistent".
- `terra::lapp(r, fun = function(x) …)` is the wrong tool: `lapp` passes each layer as a separate argument, and the function isn't vectorised. The NA check therefore always errors.

**Evidence:**
```r
check_rasters(tmin = sources(a), tmax = sources(b), prcp = sources(c3), check_nas = FALSE)$geom_report
#> "FAILURE: Geometries are inconsistent. unable to find an inherited method for function 'compareGeom' for signature '"logical", "SpatRaster"'"
check_rasters(tmin = sources(a), prcp = sources(b2))
#> Error: [lapp] cannot use 'fun'. The number of values returned is less than the number of input cells.
```
**Fix:**
- Geometry: use `all(vapply(raster_list[-1], terra::compareGeom, logical(1), x = raster_list[[1]], stopOnError = FALSE))`.
- NAs: build `terra::app(r, function(v) { k <- rowSums(is.na(v)); ifelse(k == 0, 2, ifelse(k == ncol(v), 0, 1)) })`, or more simply compare `terra::app(r, "anyNA")` across inputs.

### 11. [bug] NA handling differs between backends, and both behaviours are questionable
**Where:**
- terra: `bios_terra.R` (`app(..., na.rm = TRUE)`, `stdev(na.rm = TRUE)`); `helper_functions.R:64` (`get_window`: `sum, na.rm = TRUE`)
- tiled: `bioclim_fast.R:393-396`; `Rfast::rowmeans/rowsums/rowMaxs/rowVars`, which have no `na.rm`

**What:**
- **terra path:** a cell with one missing month silently returns the sum of the other months. bio12 is too low, and the wettest period's total is understated. It may also select the wrong period, because a period with a missing month sums lower.
- **tiled path:** the NA mask comes from the *first layer of the first variable only* (`nonaID <- which(!is.na(df[[ref_col_name]]))`). Missing values in other months or variables pass straight into Rfast, which gives NA, or an unspecified result for `rowMaxs(value = FALSE)`.
- **Net effect:** the same input gives different outputs depending on what `auto` happens to pick.

**Evidence:** terra path, with December prcp set to NA:
```
bio12 bio16 bio17
  380   190    20      # true total with Dec = 500; wettest quarter Nov-Dec-Jan = 310
```
The tiled side is **(by reading)**: Rfast wasn't available to run it.

**Fix:** choose one policy and apply it in both backends. A good default is NA when any month is NA, i.e. drop `na.rm = TRUE`. Optionally add an `na_policy` argument. In the tiled path, compute `nonaID` from `rowSums(is.na(climate_mat)) == 0`, or apply the same policy explicitly.

### 12. [bug] `derive_statistics(stats = "stdev")`: population SD on the terra path, sample SD on the tiled path
**Where:** `R/stats_terra.R:104` (`app(variable, "std")`) vs `R/stats_fast.R:576` (`Rfast::rowVars(std = TRUE)`). bio04 and bio31 use the sample SD.

**Evidence:** terra path on the 12-value series above gives `var_stdev = 7.3428`, which is the population SD. The sample SD is `7.6693`. The tiled side is **(by reading)**: Rfast `rowVars` uses n − 1.

**Fix:** use `terra::stdev(variable, pop = FALSE, na.rm = …)` in `stats_terra`, matching bio04.

### 13. [bug] `stats_fast()`: a region whose first-layer values are all ≤ 0 is treated as "empty"
**Where:** `R/stats_fast.R:208` (`sum(terra::values(target_template_rast, na.rm = TRUE)) > 0`)

**What:**
- Temperatures below 0 °C, or zero precipitation, make the check fail. The output then silently uses the *full* input extent instead of the region, with a misleading warning.
- `terra::values()` also loads the whole cropped first layer into RAM, in the backend that exists to avoid that.

**Evidence:**
```r
neg <- <12 file-backed layers with values in [-20, -1]>
derive_statistics(variable = neg, stats = "mean", user_region = terra::vect(terra::ext(2, 6, 2, 6), crs = "EPSG:4326"),
                  method = "tiled", output_dir = tempfile())
#> Warning: Cropping resulted in empty or all-NA raster, using original geometry.
```
**Fix:** use `terra::ncell(x) > 0 && terra::global(!is.na(x), "max")[[1]] == 1`, or just the `ncell` check that the other three `*_fast` files use.

### 14. [bug] A `user_region` in a different CRS works only on the tiled path
**Where:**
- `derive_bioclim.R:253, 256, 282`
- `derive_statistics.R:183-186, 211`
- `calculate_*`: `crop(x[[1]], user_region)` and `terra::crop(x, user_region, mask = TRUE)`

**What:** the `*_fast` backends reproject the region with `sf::st_transform`. The terra path and the `auto` overlap check compare raw extents. A region in, say, EPSG:3857 over EPSG:4326 rasters therefore:
- stops with "does not overlap" in `auto`;
- errors in `method = "terra"`;
- works with `method = "tiled"`.

**Evidence:**
```r
reg <- terra::project(terra::vect(terra::ext(2, 6, 2, 6), crs = "EPSG:4326"), "EPSG:3857")
derive_bioclim(bios = 12, prcp = p, user_region = reg, method = "terra", output_dir = tempfile())
#> Error : [crop] extents do not overlap
```
**Fix:** normalise `user_region` once in the wrapper. Convert it to a SpatVector, check `terra::is.valid`, project it to `terra::crs(x)` when the CRS differs, and pass that object to both backends.

### 15. [bug] `derive_statistics()`: invalid `inter_stats` values are silently ignored
**Where:** `R/stats_terra.R:63-64`

**What:** `inv_inter_stats` is computed, but the condition tests `inv_period_stats`. A typo such as `inter_stats = "max_interr"` returns only the other statistics, without any warning. The wrapper doesn't validate `inter_stats` either.

**Evidence:** `derive_statistics(variable = v, stats = "mean", inter_variable = p, inter_stats = "max_interr", method = "terra")` returns only `var_mean`, with no error.

**Fix:** change line 64 to `if (length(inv_inter_stats) > 0)`. Better, validate `stats`, `inter_stats` and `period_stats` in `derive_statistics()` itself, before any work, with `match.arg(several.ok = TRUE)`.

### 16. [bug] `period_length` isn't validated on the terra path
**Where:** `R/helper_functions.R:45-68` (`get_window`), `bioclim_terra.R`, `stats_terra.R`. The tiled `compute_periods` checks only `> 0` and `<= n_units`.

**What:**
- **13 on 12 layers:** each window silently includes one month twice.
- **2.5:** windows hold 2 months, but `tmp` is divided by 2.5, so the "mean" is 20% low. The tiled path averages over 2, so the backends differ.
- **0:** "subscript out of bounds".

**Evidence:** `get_window(r, 13, TRUE)` gives a first-window sum of 790, which is 780 for all twelve months plus January again. `get_window(r, 2.5, TRUE)` sums 2 layers.

**Fix:** in every exported wrapper, require `period_length` to be a single whole number with `1 <= period_length <= nlyr`. Pass `circular = FALSE` only when `period_length < nlyr`.

### 17. [bug] Tiled backends assume one single-band file per layer (partly run)
**Where:**
- `derive_bioclim.R:294`
- `derive_statistics.R:233`
- `calculate_average.R:147`, `calculate_sum.R:147`, `calculate_roll.R:195` (`paths = terra::sources(x)`)
- workers do `terra::rast(paths)`

**What:** `sources()` returns *files*, not layers. Two problems follow:
- **derive_\*:** a 12-layer view of a multi-band file, such as the package's own `prcp.tif` example, gives 1 path. The run then stops with "prcp_path length error". `auto` can route such input to tiled on its own (A2).
- **calculate_\*:** `rast(paths)` reopens all bands. For `x <- r[[13:24]]` from a 60-band file, `aggregate_fast` / `roll_fast` get 60 columns. `index` then refers to columns 1–12, which are the *wrong* layers, with no error.

**Evidence:**
```r
m <- terra::rast(system.file("extdata", "prcp.tif", package = "fastbioclim"))[[1:12]]
terra::nlyr(m); basename(terra::sources(m)); terra::nlyr(terra::rast(terra::sources(m)))
#> 12  "prcp.tif"  60
derive_bioclim(bios = 12, prcp = m, method = "tiled", output_dir = tempfile())
#> Error ... prcp_path length error: expected 12, received 1.
```
The silent wrong-layer outcome for `calculate_*` is **(by reading)**: it needs Rfast to complete.

**Fix:** pass `(file, band)` pairs. Use `terra::sources(x, bands = TRUE)` in recent terra, or build `terra::rast(file, lyrs = band)`. In workers, build the stack from those pairs. Alternatively, accept any SpatRaster by writing a VRT, `terra::vrt()`, that holds exactly the selected layers, and pass the VRT path.

### 18. [bug] An out-of-range static period index returns the `min_idx`/`max_idx` column as data (by reading)
**Where:** `R/bios_fast.R`: `bio08_fast` … `bio35_fast` period extractors, e.g. `tperiod[cbind(seq_len(nrow(tperiod)), pperiod_max_idx)]`

**What:**
- `var_periods()` appends `min_idx` and `max_idx` as the last two columns.
- The extractors only *warn* when an index is greater than `num_period_cols`, then index anyway.
- **Example:** a static `warmest_period` made with `circular = TRUE` (values 1–12) used in a `circular = FALSE` run (10 periods). Index 11 returns the *period number* in `min_idx` as if it were a temperature.
- The terra path (`selectRange`) returns NA for the same input, and `stats_fast` masks invalid indices correctly.

**Fix:** mask invalid indices to NA, as `stats_fast.R:589-597` does. Factor that into one helper used by all period extractors.

### 19. [bug] `sf::sf_use_s2(FALSE)` changes the user's global sf setting (by reading)
**Where:** `bioclim_fast.R:330`, `stats_fast.R:279`, `aggregate_fast.R:147`, `roll_fast.R:122`

**What:** after any tiled run, S2 spherical geometry stays off for the rest of the session. That changes the results of the user's own `st_intersection`, `st_area` and similar calls on geographic data. (The tiled runs here stopped at `qs2`, before this line.)

**Fix:** `old <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(old), add = TRUE)`. Better, build the tiles with terra (see P2) and drop the sf geometry step.

### 20. [bug] `calculate_roll()` tiled returns nothing when `name_template` doesn't start with `{prefix}` (by reading)
**Where:** `R/calculate_roll.R:129, 211`; `write_layers.R:68-73`

**What:** `write_layers()` keeps only intermediate files whose names start with `file_pattern = output_prefix`. A valid template such as `"{fun}_{prefix}_{idx_unit}"` produces files that don't match, so the call stops with "No intermediate '.qs2' files found matching…". The terra path works.

**Fix:** the temp directory belongs to this call alone. Give `write_layers()` the list of expected names (`all_output_names`) instead of a prefix.

---

## Risks

### 21. [risk] Tiled outputs lose the CRS identity
**Where:** `*_fast.R` (`crs_txt <- terra::crs(ref, proj = TRUE)`), `write_layers.R:95-100`

**What:** the template is rebuilt from a PROJ.4 string, which drops the WKT2 definition, the authority code and datum details.

**Evidence:**
```
EPSG:3035 -> code after round-trip: NA  same.crs: FALSE
EPSG:2193 -> code after round-trip: NA  same.crs: FALSE
EPSG:4326 -> code after round-trip: NA  same.crs: TRUE
```
Outputs may then fail `compareGeom()` against their own inputs, or against terra-path outputs.

**Fix:** store `terra::crs(ref)`, which is WKT2.

### 22. [risk] Intermediate directories are left in `output_dir` after any failure
**Where:** `bioclim_fast.R:209`, `stats_fast.R:163`, `aggregate_fast.R:4`, `roll_fast.R:4`

**What:** these directories are deleted only at the end of `write_layers()`. A failed or interrupted run leaves `bioclima_qs_*` (etc.) in the user's output folder. A large run can leave many gigabytes there.

**Evidence:** after a failed `method = "tiled"` run, `list.files(output_dir)` returns `"bioclima_qs_25bed92657d"`.

**Fix:** create the directory under `tempdir()`, or register `on.exit(unlink(dir, recursive = TRUE), add = TRUE)` in the wrapper. Also, `derive_statistics.R:252-253` creates a `stats_dir` and passes it, but `stats_fast()` ignores it and makes its own.

### 23. [risk] Overwrite and output-directory behaviour
- **Missed pre-check.** `derive_statistics.R:126` builds the expected names from `stats`. `"cv_cli"` is written as `<prefix>_cv.tif`, so the pre-check misses an existing file. The run then does all the work and fails inside `writeRaster` ("file exists"). Confirmed.
- **`overwrite` ignored.** `write_layers.R:143` always uses `overwrite = TRUE`, so on the tiled path the `overwrite` argument does nothing past the wrapper's name check.
- **Failed writes still returned.** `write_layers()` only warns on a failed write, but still returns that file's path. The wrapper's `terra::rast(...)` then fails with a less useful error.
- **Default `tempdir()` collides.** With the default `output_dir = tempdir()`, two analyses in one session collide. The second call, even with different data, errors with "Output files already exist" (confirmed). The cause is that every output (`bio12.tif`, `var_mean.tif`, …) goes into the root of the session temp directory.

**Fix:**
- Build the expected file names with the same helper that names the outputs.
- Pass `overwrite` through to `write_layers`.
- Make the default `file.path(tempdir(), "fastbioclim", <unique id>)`, or require `output_dir`.

### 24. [risk] Input validation and error messages
- **`bios` isn't validated in the terra path.**
  - `bios = 36` fails with `value for 'bio36' not found`.
  - `bios = c(12, 12)` fails with `[writeRaster] duplicate filenames`.
  - Both are confirmed. `bioclim_fast` validates `bios` but the wrapper doesn't.
  - Fix: `stopifnot(is.numeric(bios), all(bios %in% 1:35)); bios <- sort(unique(bios))` in `derive_bioclim()`.
- **`message(..., call. = FALSE)`.**
  - At `derive_bioclim.R:223-227`, `message()` has no `call.` argument, so "FALSE" is appended to the text.
  - The message also ignores `verbose = FALSE`. Confirmed output: `...vs 'Warmest Month').FALSE`.
- **Errors raised through `do.call()`** (`derive_bioclim.R:286`, `derive_statistics.R:224`).
  - These print the *entire deparsed function body* as the error call, often hundreds of lines.
  - Use `stop(..., call. = FALSE)` in the backends, or call them directly instead of through `do.call`.
- **`stop(c(..., "i" = ...))`** (e.g. `derive_bioclim.R:175-179`) uses cli-style names that base `stop()` ignores. Either use `cli::cli_abort` or plain strings.
- **`tile_degrees` is in CRS units.** For projected data it means metres: 5 m tiles over a 100 km UTM raster would be `st_make_grid` with 4×10⁸ cells. Validate it against `terra::is.lonlat()`, or express tile size in cells.
- **`inMemory` checks disagree.** `derive_*` use `any(inMemory)`, but `calculate_*` use `all(inMemory(x))`. A raster that is partly in memory passes `calculate_*`, and its `sources()` then contains `""`.

---

## The `auto` method

### A1. [bug] It fails on supported terra versions
See finding 5.

### A2. [risk] It doesn't check tiled preconditions, and has no fallback
**Where:** `derive_bioclim.R:231-268, 287-292` (same pattern in the other wrappers)

**What:** when the estimate says "doesn't fit", `auto` switches to tiled, and only then checks what tiled requires:
- file-backed inputs;
- for `derive_*`, one single-band file per layer (finding 17).

So `auto` *errors* on large in-memory or multi-band inputs, which are common. Meanwhile the terra path would have run, because terra already processes out-of-core in chunks with temporary files when data don't fit.

**Fix:** evaluate whether tiled is feasible first. If tiled isn't possible, keep terra and say so with `message()`. Choose tiled only when it's both feasible and expected to help.

### A3. [improvement] The memory estimate leaves out the real peaks
**What:**
- **Inputs only.** `n = length(inputs) * nlyr(inputs[[1]])` counts input layers. It leaves out the period stacks the terra path builds (`wet`, `tmp`, `prad`, `pmois`, about +1× layers per variable) and the derived `tavg`.
- **Tiled-path peaks ignored.** The estimate doesn't model the tiled path's own peaks:
  - `write_layers()` allocates a full `ncell`-length double vector per output (P1);
  - each future worker holds several copies of a tile matrix (P4).
- **Not a true fork.** Because terra itself degrades gracefully, `fits_mem == 0` doesn't mean "terra will fail". It means "terra will use temporary files".

**Fix:** decide on the expected runtime and feasibility, not only on `fits_mem`. At minimum:
- include the period layers and `future::nbrOfWorkers()` in the estimate;
- document that `auto` → tiled needs about `ncell * 8` bytes free per output during assembly.

### A4. [bug] It ignores the CRS of `user_region`
See finding 14. `terra::relate(ext(user_region), ext(raster), "intersects")` compares extents that may be in different CRSs.

---

## Tiled path: performance and parallel scaling

### P1. [improvement] Assembly is neither out-of-core nor parallel
**Where:** `write_layers.R:115-128`

**What:** each output variable needs `rvals <- rep(NA_real_, n_target_cells)`, which means reading every tile's `.qs2` into one full-size vector. For a global 30″ grid (43200 × 21600) that is about 7.5 GB per variable, assembled one variable at a time and single-threaded. The rows are then written 16 at a time. This step is the RAM ceiling of the "memory-safe" backend.

**Fix:** use tiles aligned to the cell grid (P2) and write each tile's result straight into its window of the output. Two ways:
- write one small GeoTIFF per tile, then `terra::vrt()` plus one `writeRaster`;
- or use horizontal strips that can be streamed with `writeStart` / `writeValues` in row order.

Either way, the `.qs2` round-trip (P3) and the full-size vector both go away.

### P2. [improvement] Polygon tiles cause duplicated cells and costly extraction
**Where:** `bioclim_fast.R:330-339` (same in the 3 others), plus `exact_extract()` per tile.

**What:**
- **Duplicates.** Tiles come from `st_make_grid(cellsize = tile_degrees)` anchored at the bbox, so tile edges don't follow the raster's cell edges. `exact_extract` returns every partly covered cell, so boundary cells are read, computed and written two to four times. It also computes coverage fractions that are never used.
- **Measured.** On a 0.3° grid with 5° tiles, 22 155 cell extractions were made for 20 000 cells, **10.8 % duplicates**. Extraction took 1.13 s with polygon tiles, against 0.47 s for reading the same windows with `terra::crop(r, ext)`.
- **With `user_region`.** Tiles become `st_intersection(region, grid)` polygons, which are more complex and slower, and they also cause finding 1.

**Fix:** define tiles as row/column windows of the raster, e.g. `terra::getTileExtents(r, n)` or computed row ranges, and read them with `terra::crop(r, ext)` / `readValues(row, nrows, col, ncols)`. Apply `user_region` as a mask rasterised once (`terra::rasterize`, cell-centre rule, same as the terra path). Tiles that contain no region cells are then skipped cheaply.

### P3. [improvement] Every result is written twice, in a bulky format
**What:**
- Each requested variable gets one `.qs2` per tile, holding `cbind(value, cell)` as doubles. That doubles the size, because cell ids are stored as doubles.
- Everything is later read back and written again as GeoTIFF.
- With 35 bios and 300 tiles that's 10 500 small files, all in `output_dir`.

**Fix:** use P1. Failing that, store `cell` as integer, or write one file per tile holding all variables.

### P4. [improvement] Parallel scaling
**Where:** `future_lapply(seq_len(ntiles), …)` in the 4 `*_fast` files.

**What:**
- **Load imbalance.** The default `future.scheduling = 1` cuts tiles into `nbrOfWorkers()` contiguous chunks. Tiles are ordered by grid row, and ocean tiles return almost at once. One worker can get all the empty southern tiles while another gets all the land.
- **Peak memory grows with worker count.** Each worker holds, for one tile:
  - the `exact_extract` data.frame;
  - `data.frame.to_matrix`;
  - per-variable `check_evar` copies;
  - period matrices.

  That is about 4–5 copies of `ncell_tile × nlyr_total × 8` bytes. At 30″ resolution and 5° tiles with 4 variables, it comes to roughly 0.6–0.8 GB per worker. Neither `auto` nor the documentation mentions this.
- **Small overheads.** Each tile calls `gc()` (a full collection) and reopens every input file with `terra::rast(paths)`.
- **No progress shown.** The `progressr::progressor()` is created, but `with_progress()` is commented out (`bioclim_fast.R:366`). Users see nothing unless they set `handlers(global = TRUE)`, which isn't documented in the function help.

**Fix:**
- use `future.scheduling = Inf`, or `future.chunk.size = 1`, for dynamic load balancing;
- shuffle the tile order;
- drop the per-tile `gc()`;
- size tiles from `available RAM / nbrOfWorkers()` instead of a fixed `tile_degrees`;
- document the progress handler.

### P5. [improvement] Repeated work across tiles
`compute_periods()` is called again for every tile (`bioclim_fast.R:506`, `stats_fast.R`), but it depends only on `n_units`, `period_length` and `circular`. Compute it once and send it to the workers. `var_periods()` builds its matrix with `sapply` over periods and then `cbind`s the result. That's fine for 12 periods, but it allocates two copies.

---

## Maintainability and API consistency

### M1. [improvement] Four copies of the tiling, region and template setup
`bioclim_fast.R:186-339`, `stats_fast.R:140-290`, `aggregate_fast.R:8-155` and `roll_fast.R:8-130` hold about 150 lines each of the same steps: geometry check, region normalisation, crop template, offsets, `define_translate`, `template_info.qs2`, `st_make_grid`. The copies have already drifted apart:
- only `stats_fast` has the `sum(values) > 0` check (finding 13);
- `stats_fast` exports a non-existent global `"period"` and omits `"period_length"` (`stats_fast.R:332`). That works only because the closure environment happens to carry it.
- `aggregate_fast` exports `"path_variables"` and `"translate_cell"`, which may not exist.
- the `noMap` handling differs (finding 2).

**Fix:** one internal `prepare_tiles(paths, user_region, tile_spec)` that returns `{template_info, tiles, translate}`, plus one `run_tiles(tiles, worker_fun)`.

### M2. [improvement] Duplicated wrappers
- `derive_bioclim()` and `derive_statistics()` repeat the overwrite pre-check, the `auto` block, the CRS/region handling and the dispatch almost line for line.
- `calculate_average()` and `calculate_sum()` are about 100 identical lines that differ only in `aggregation_type`.
- **Fix:** an internal `choose_method()`, `normalise_region()` and `check_outputs()`, and one `calculate_aggregate(type)` behind the two exported functions.

### M3. [improvement] Bio dependency logic is written twice, with different rules
- `bioclim_terra.R:59-132` uses ad-hoc `req_*` vectors plus `exists()`.
- `bioclim_fast.R:65-117` and `:552-582` use a `needs` list, a fixed-point loop and a 30-line `can_calculate` cascade.
- The two disagree, which is the cause of findings 6–9. `bioclim_fast` is a single ~750-line function containing a 140-line `if/else if` chain over `bio_num`.
- **Fix:** one declarative table:

  ```r
  bio_spec <- list(bio08 = list(vars = "tavg", period_of = "prcp", index = "wettest_period", which = "max"), ...)
  ```

  Then use a small resolver that both backends call, and dispatch with `switch()` or a function lookup instead of the `if` chain.

### M4. [improvement] 70 near-identical per-bio functions
- `bios_terra.R` has 35 functions that call `app(...)` or `selectRange(...)` and set a name.
- `bios_fast.R` has 35 functions in which the select-by-unit-index block (bio05/06/13/14/21/22/29/30) and the select-by-period block (16 functions) are pasted copies.
- **Fix:** three or four generic helpers, `summarise_units(x, fun)`, `select_unit(x, idx, fun)` and `select_period(periods, idx)`, plus a name table. This also fixes finding 18 in one place.

### M5. [improvement] Wrappers depend on implicit partial argument matching
`derive_bioclim()` passes `tmin = <paths>` to `bioclim_fast(tmin_path = …)`, and `derive_statistics()` passes `variable = <paths>` to `stats_fast(variable_path = …)`. Both work only through R's partial matching of argument names. This breaks silently if another argument that starts with `tmin`… is ever added, and it triggers `R CMD check --as-cran` partial-match warnings when that option is enabled. Name the arguments explicitly.

### M6. [improvement] API consistency
- **Static-index rules.** Static indices are SpatRasters in `...` for both `derive_*` functions. `bioclim_fast` then silently ignores multi-file sources, because it requires `length(path) == 1` and only *warns* "invalid/not found. Ignoring". If a user asked for a static index, failing to use it should be an error.
- **`user_region` types.** The terra path accepts anything `crop()` accepts. The tiled path accepts only `sf`/`sfc`/`SpatVector`.
- **Output names vs argument values.** Statistic names don't match the values users pass (`"cv_cli"` gives `_cv`). `calculate_*` add `_NN` suffixes, while `calculate_roll` uses a glue template.
- **`verbose` not always respected.** See finding 24.
- **Return values.** `derive_*` and `calculate_*` return SpatRasters, but layer order on the tiled path comes from `write_layers`' sort (numeric for `"bio"`, alphabetical otherwise), not the order requested.

### Tests
Every test that touches a backend replaces it with `mockery` stubs (`tests/testthat/*`), so neither the tiled code nor the terra and tiled numerics is ever executed. One small parity test on a 20×20×12 raster would have caught findings 1, 2, 3, 9, 11 and 12. It would run both `method = "terra"` and `method = "tiled"` and compare with `all.equal()`, with and without a non-rectangular `user_region` and with some NA cells. It would take seconds and should be added before any refactor (M1–M4).
