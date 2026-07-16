load_test_spectra_matrix <- function(max_points = 5000) {

  data_file <- testthat::test_path(
    "test-data",
    "loaded_spectra.RData"
  )

  load(data_file)

  # object stored in RData is called spectra
  NMRData1D_list <- spectra


  intensity <- lapply(
    NMRData1D_list,
    function(x) {

      Re(x@processed$intensity)

    }
  )


  ppm <- NMRData1D_list[[1]]@processed$direct.shift


  mat <- do.call(
    rbind,
    intensity
  )


  rownames(mat) <- names(NMRData1D_list)


  colnames(mat) <- ppm


  # reduce size for tests
  if(ncol(mat) > max_points){

    keep <- seq(
      1,
      ncol(mat),
      length.out = max_points
    )

    mat <- mat[, keep, drop = FALSE]

  }


  mat

}
