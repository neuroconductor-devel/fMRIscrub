test_that("pscrub works", {
  psx <- testthat::expect_warning(fMRIscrub:::pscrub_multi(
    Dat1,
    projection = c("ICA", "ICA_kurt"), # "all"
  ))
  myplot <- fMRIscrub:::plot.scrub_projection_multi(psx)

  psx <- testthat::expect_warning(fMRIscrub:::pscrub_multi(
    Dat2,
    projection = c("PCA", "ICA_kurt"),
    kurt_quantile = .90,
    cutoff = 5,
    verbose = TRUE
  ))
  myplot <- fMRIscrub:::plot.scrub_projection_multi(psx, title = "My Plot")

  psx <- testthat::expect_warning(pscrub(Dat1))
  plot(psx)

  # psx <- pscrub(
  #   matrix(rnorm(9000), nrow=30), "fusedPCA", 1, center=FALSE, PESEL=FALSE, kurt_quantile=.8,
  #   full_PCA = TRUE, cutoff=5, verbose=TRUE
  # )
  # plot(psx)

  psx <- testthat::expect_warning(pscrub(
    matrix(rnorm(10000), ncol=50)
  ))

  psx <- testthat::expect_warning(pscrub(
    matrix(rnorm(10000), nrow=100) + 100, nuisance=fMRItools::dct_bases(100, 2)
  ))

  psx <- testthat::expect_warning(pscrub(
    Dat2, projection="PCA", nuisance=cbind(1, fMRItools::dct_bases(nrow(Dat2), 12)),
    comps_mean_dt=2, comps_var_dt=2, get_dirs=TRUE, get_outliers=FALSE
  ))
  plot(psx)

  # scrub_xifti(ciftiTools::read_xifti("../Data/rfMRI_REST1_LR_Atlas_MSMAll_6k.dtseries.nii", brainstructures="left"))

})