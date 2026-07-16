# Broadly-used utility functions

#------------------------------------------------------------------------
#' Safely read data from file
#'
#' Reads data from file if it exists and is readable. Outputs sane
#' error message if either of these conditions is not met.
#'
#' @param path Path to file.
#' @param data Character determining read function to use: 'char' for
#'             readChar or 'bin' for readBin.
#' @param ... Arguments passed to underlying read function.
#
#' @return Character or numeric vector depending on input data.
#'
#' @export
safe_read <- function(path, data, ...) {

  # Checking if path is a directory
  if (dir.exists(path)) {
    msg <- sprintf('Path "%s" points to a directory not a file', path)
    stop(msg, call. = FALSE)
  }

  # Checking if path exists
  if (! file.exists(path)) {
    msg <- sprintf('Path "%s" does not point to a file', path)
    stop(msg, call. = FALSE)
  }

  # Choosing read function
  if (data == 'char') {
    if (!'nchars' %in% list(...)) { 
      f_read <- function(path, ...) {
        nchars <- file.info(path)$size
        readChar(path, nchars, ...)
      }
    }
    else {
      f_read <- function(path, ...) {
        readChar(path, ...)
      }
    }
  }
  else if (data == 'bin') {
    f_read <- function(path, ...) {
      readBin(path, ...)
    }
  }
  else {
    msg <- '"data" argument must be either "char" or "bin"'
    stop(msg, .call = FALSE)
  }
  
  out <- f_read(path, ...)

  return(out)
}

#------------------------------------------------------------------------
#' Validate and extract acqus/procs parameters
#'
#' Checks to see whether required acqus or procs parameters are present in 
#' main list or appriate sublist. Outputs the required parameters if
#' all are present and issues an error otherwise. This function is mainly
#' intended for internal use.
#'
#' @param param.list A list possibly containing parameters.
#' @param required.param A vector of required parameters.
#' @param sublist Name of sublist that parameters may be stored in
#'
#' @return A list containing only the required parameters
#'
#' @export
validate_param <- function(param.list, required.param, sublist = NULL) {

    missing <- ! required.param %in% names(param.list)

    # If all parameters are missing, check to see if the parameters
    # are nested within a sublist.
    if (all(missing) && !is.null(sublist)) {
      param.list <- param.list[[sublist]]
      missing <- ! required.param %in% names(param.list)
    }

    # If any parameters are still missing, issue error
    if (any(missing)) {
      msg <- sprintf('The following parameters are missing: %s',
                     paste(required.param[missing], collapse=', '))
      stop(msg)
    }

    # Otherwise, returning the required parameters
    param.list[required.param]
}

#------------------------------------------------------------------------
#' Find the indeces of closest matching values
#'
which_approx <- function(x, y) {
  n <- length(y)
  out <- rep(0, n) 

  for (i in 1:n) {
    abs.diff <- abs(x - y[i])
    out[i] <- which(abs.diff == min(abs.diff))[1]
  }
  out
}

#------------------------------------------------------------------------
#' Find the indeces of closest matching values within a matrix
#'
which_approx_2d <- function(x1, x2, y1, y2) {
  if ( (length(x1) != length(x2)) || (length(y1) != length(y2)) ) {
    msg <- 'The length of x1 should match x2, similarly with y1 and y2.'
    stop(msg, call. = FALSE)
  }

  n <- length(y1)
  out <- rep(0, n) 

  # Picking an arbitrary x and y range
  for (i in 1:n) {
    diff1 <- abs(x1 - y1[i])
    diff2 <- abs(x2 - y2[i])
    abs.diff <- diff1*diff2

    out[i] <- which(abs.diff == min(abs.diff))
  }
  out
}
