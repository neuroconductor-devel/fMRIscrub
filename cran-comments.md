## Test environments

* Linux x86_64-pc-linux-gnu, R 4.0.4
* Mac x86_64-apple-darwin17.0, R 4.1.1
* Windows x86_64-w64-mingw32, R 4.1.1

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

`templateICAr` works fine with this new version.

## Tests

Passes all the tests in `tests/testthat.R`.

## Submission after being archived

This is the second submission since `fMRIscrub` was removed from CRAN on May 2, 2022 due to a package in Suggests being no longer available. The current submission resolves this issue. Since the first submission after being archived, we have also removed or replaced problematic code from the package.

## Resubmission

The first submission had three NOTES:

  Possibly misspelled words in DESCRIPTION:
      Afyouni (28:17)
      DVARS (27:60, 33:42)
      Muschelli (31:29)
      Pham (26:68)
      VARianceS (28:5)
      aCompCor (30:53)
      al (26:76, 29:73, 31:42)
      denoising (26:31)
      detrending (32:46)
      et (26:73, 29:70, 31:39)

These words are correctly spelled.

    Suggests or Enhances not in mainstream repositories:
      glmgen

  Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-x86_64
  Check: package dependencies, Result: NOTE
    Package suggested but not available for checking: 'glmgen'

We have removed the `glmgen` dependency.

  Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-x86_64
  Check: Rd files, Result: NOTE
    checkRd: (-1) fusedPCA.Rd:49: Escaped LaTeX specials: \$ \$

We have remoted the escaped LaTeX specials.

# Resubmission 2

The previous submission failed without suggested packages:

  * checking examples ... [3s] ERROR
  Running examples in 'fMRIscrub-Ex.R' failed
  The error most likely occurred in:

  > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
  > ### Name: bandstop_filter
  > ### Title: Bandstop filter
  > ### Aliases: bandstop_filter
  >
  > ### ** Examples
  >
  > library(gsignal)
  Error in library(gsignal) : there is no package called 'gsignal'
  Execution halted

We have fixed the issue by checking for availability of suggested packages prior to running examples.