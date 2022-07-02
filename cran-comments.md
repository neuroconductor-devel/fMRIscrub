## Test environments

* Linux x86_64-pc-linux-gnu, R 4.0.4
* Mac x86_64-apple-darwin17.0, R 4.1.1

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

`templateICAr` works fine with this new version.

## Tests

Passes all the tests in `tests/testthat.R`.

## Submission after being archived

This is the second submission since `fMRIscrub` was removed from CRAN on May 2, 2022 due to a package in Suggests being no longer available. The current submission resolves this issue. Since the first submission after being archived, we have also removed or replaced problematic code from the package.