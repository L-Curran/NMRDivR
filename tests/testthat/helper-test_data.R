# tests/testthat/helper-test_data.R

load_test_spectra_matrix <- function(n = 6) {

  data_file <- testthat::test_path(
    "test-data",
    "loaded_spectra.RData"
  )

  load(data_file)

  # object saved as "spectra"
  NMRData1D_list <- spectra


  # subset spectra
  NMRData1D_list <- NMRData1D_list[seq_len(min(n,length(NMRData1D_list)))]


  ppm <- purrr::map_dbl(
    NMRData1D_list[[1]]@processed$direct.shift,
    identity
  )


  intensity <- purrr::map_dfc(
    NMRData1D_list,
    function(x){

      Re(
        x@processed$intensity
      )

    }
  )


  intensity <- as.matrix(intensity)

  intensity <- t(intensity)

  colnames(intensity) <- ppm

  rownames(intensity) <- paste0(
    "sample_",
    seq_len(nrow(intensity))
  )


  intensity
}
