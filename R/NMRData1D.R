# Definition of a class structure for 1D NMR data.

#------------------------------------------------------------------------
#' A class combining NMR data and scan parameters
#'
#' Although it's likely to change in the near future, this class currently
#' serves to combine processed 1D data from the pdata folder of an experiment
#' with the acqus and procs parameter files.
#'
#' @slot processed 1r/1i data stored as a data.frame with direct.shift and
#'                 intensity as columns.
#' @slot acqus A list of acqus parameters.
#' @slot procs A list of procs parameters.
#' @slot product A vector representation of an apodization function that was 
#'               used to generate the intensity data (in this context, 
#'               zero-filling is considered to be the equivalent of 
#'               apodization by a step function). To facilitate convolution, 
#'               the length of this vector is set to 2*n - 1 where n is the 
#'               number of spectrum points.
#' @slot convolution The Fourier transform of the product slot, suitable
#'                   for convolution with the spectrum intensity data.
#'
#' @name NMRData1D-class
#' @export
NMRData1D <- setClass("NMRData1D", 
                      contains = "NMRData",
                      slots = c(product = 'numeric',
                                convolution = 'numeric'))

# Validity testing consists of simply checking the processed data.frame columns
validNMRData1D <- function(object) {

  valid <- TRUE
  msg <- c()

  processed <- object@processed
  product <- object@product
  convolution <- object@convolution

  # Checking column names
  valid.columns <- c('direct.shift', 'intensity')
  new.msg <- sprintf('Processed data must have two columns: %s',
                     paste(valid.columns, collapse = ', '))

  if (! identical(valid.columns, colnames(processed)) ) {
    valid <- FALSE
    msg <- c(msg, new.msg)
  }

  # Checking convolution length
  new.msg <- paste('If not empty, "convolution" and "product" slots must',
                   'have 2*n - 1 elements where n is the number of points in',
                   'the "processed" data.frame')

  if ( (length(product) > 0) || (length(convolution) > 0) ) {

    n <- nrow(processed)
    if ( (length(product) != (2*n-1)) || (length(convolution) != (2*n-1)) ) {
      valid <- FALSE
      msg <- c(msg, new.msg)
    }

  }

  if ( valid ) TRUE
  else msg
}

setValidity("NMRData1D", validNMRData1D)

#------------------------------------------------------------------------
#' Constructor for generating an NMRData1D object
#'
#' Loads data in JCAMP file or Bruker directory and uses it to initialize an 
#' NMRData1D object.
#'
#' @param path Path to a Bruker scan directory or JCAMP file.
#' @param procs.number Specifies processing file number to use when loading
#'                     from a Bruker scan directory. Defaults to the smallest 
#'                     number in the pdata directory. Ignored if loading from
#'                     a JCAMP file.
#' @param blocks.number Specifies block number to use when loading from a
#'                      JCAMP file. Defaults to the first block encountered.
#'                      Ignored if loading from Bruker directory.  
#' @param ntuples.number Specifies ntuple entry number to use when loading 
#'                       from a JCAMP file. Defaults to the first ntuple 
#'                       entry encountered. Ignored if loading from Bruker 
#'                       directory.  
#' 
#' @return An NMRData1D object containing the 1r/1i processed data as well
#'         as the procs and acqus parameters.
#'
#' @export
nmrdata_1d <- function(path, procs.number = NA, 
                       blocks.number = NA, ntuples.number = NA) {

  # If file exists at path, treating import as JCAMP, otherwise, Bruker scan
  if ( file.exists(path) && !dir.exists(path) ) {
    # Load jcamp
    jcamp <- read_jcamp(path, process.tags = TRUE, process.entries = TRUE)

    # If block number specified, check if it's valid
    if (! is.na(blocks.number) ) {
      # Double check that specified block exists
      msg <- sprintf("Specified block number not found in %s", path)
      if ( blocks.number > length(jcamp$blocks) ) stop(msg)
    }
    else {
      blocks.number <- 1
    }

    jcamp$blocks <- jcamp$blocks[blocks.number]

    # Check to make sure that ntuples exist
    msg <- 'Import from JCAMP file currently limited to NTUPLES entries'
    if (! 'ntuples' %in% names(jcamp$blocks[[1]]) ) stop(msg)

    # If ntuple number specified, check if it's valid
    if (! is.na(ntuples.number) ) {
      # Double check that specified block exists
      msg <- sprintf('Specified NTUPLES number not found in block %i of %s', 
                     blocks.number, path)
      if ( ntuples.number > length(jcamp$blocks[[1]]) ) stop(msg)
    }
    else {
      ntuples.number <- 1
    }

    jcamp$blocks[[1]]$ntuples <- jcamp$blocks[[1]]$ntuples[ntuples.number]

    # Flattening file
    jcamp.flat <- flatten_jcamp(jcamp)

    # Extracting processed data from ntuples
    descriptors <- jcamp.flat$ntuples$descriptors
    pages <- jcamp.flat$ntuples$pages

    # Checking variable names
    variables <- tolower(descriptors$var.name)
    real.index <- which(str_detect(variables, 'spectrum.*real'))
    imag.index <- which(str_detect(variables, 'spectrum.*imag'))

    # If both real and imaginary data isn't there, abort
    msg <- 'Import from JCAMP file currently limited to real/imaginary spectra'
    if ( length(c(real.index, imag.index)) < 2 ) stop(msg)

    # Double check that the first entry is frequency
    msg <- 'Import from JCAMP file currently limited to frequency abscissa'
    if (! str_detect(variables[1], 'freq') ) stop(msg)

    # Proceed to extract data
    real.data <- pages[[real.index - 1]]
    imag.data <- pages[[imag.index - 1]]

    # Checking that frequency is the same
    real.frequency <- real.data[, 1]
    imag.frequency <- imag.data[, 1]
    msg <- 'Mismatch in real and imaginary frequency, likely parsing error'
    if ( any(abs(real.frequency - imag.frequency) > 1e-6) ) stop(msg)

    # Starting with raw values
    frequency <- real.frequency
    real.data <- real.data[, 2]
    imag.data <- imag.data[, 2]

    # Scaling if factors is non zero
    scale <- descriptors$factor[1]
    if (scale > 1e-6) frequency <- frequency*scale

    scale <- descriptors$factor[real.index]
    if (scale > 1e-6) real.data <- real.data*scale

    scale <- descriptors$factor[imag.index]
    if (scale > 1e-6) imag.data <- imag.data*scale

    # Offsetting if max-min difference is non zero
    d.max <- descriptors$max[1]
    d.min <- descriptors$min[1]
    if ( (d.max - d.min) > 1e-6 ) {
      frequency <- frequency - max(frequency) + d.max
    }

    d.max <- descriptors$max[real.index]
    d.min <- descriptors$min[real.index]
    if ( (d.max - d.min) > 1e-6 ) {
      real.data <- real.data - max(real.data) + d.max
    }

    d.max <- descriptors$max[imag.index]
    d.min <- descriptors$min[imag.index]
    if ( (d.max - d.min) > 1e-6 ) {
      imag.data <- imag.data - max(imag.data) + d.max
    }

    # Doing one final check on the direct shift to check on offset
    direct.shift <- frequency/jcamp.flat$sf
    
    delta <- jcamp.flat$offset - max(direct.shift)
    if ( abs(delta) > 1e-6 ) direct.shift <- direct.shift + delta

    # Finally, combine the data
    intensity <- complex(real = real.data, imaginary = imag.data)
    processed <- data.frame(direct.shift = direct.shift,
                            intensity = intensity)

    # Returning class object
    new("NMRData1D", processed = processed, parameters = jcamp.flat,
                     procs = list(), acqus = list())       
  }
  else {
    # First, loading procs parameters
    procs <- read_procs_1d(path, procs.number)

    # Using the procs file to load the processed data
    processed <- read_processed_1d(path, procs, procs.number)

    # Finally, loading the general acquisition parameters
    acqus <- read_acqus_1d(path)

    # Returning class object
    new("NMRData1D", processed = processed, parameters = list(),
                     procs = procs, acqus = acqus)
  }
}

#========================================================================>
# Basic conversions and slot access functions (inherited from NMRData)
#========================================================================>

#------------------------------------------------------------------------

#' @rdname processed 
#' @export
setMethod("processed", "NMRData1D", 
  function(object) callNextMethod())

#' @rdname processed-set
#' @export
setReplaceMethod("processed", "NMRData1D", 
  function(object) callNextMethod())

#------------------------------------------------------------------------

#' @rdname parameters
#' @export
setMethod("parameters", "NMRData1D", 
  function(object) callNextMethod())

#' @rdname parameters-set
#' @export
setReplaceMethod("parameters", "NMRData1D", 
  function(object) callNextMethod())

#------------------------------------------------------------------------

#' @rdname procs 
#' @export
setMethod("procs", "NMRData1D", 
  function(object) callNextMethod())

#' @rdname procs-set
#' @export
setReplaceMethod("procs", "NMRData1D", 
  function(object) callNextMethod())

#------------------------------------------------------------------------

#' @rdname acqus 
#' @export
setMethod("acqus", "NMRData1D", 
  function(object) callNextMethod())

#' @rdname acqus-set
#' @export
setReplaceMethod("acqus", "NMRData1D", 
  function(object) callNextMethod())



#========================================================================>
# Defining list and data.frame like behaviour
#========================================================================>

#' @rdname as.list.NMRData
#' @export
setMethod("as.list", "NMRData1D",  
  function(x) {
    callNextMethod()
  })

#------------------------------------------------------------------------

#' @rdname as.data.frame.NMRData
#' @export
setMethod("as.data.frame", "NMRData1D",  
  function(x) {
    callNextMethod()
  })



#========================================================================>
# Non-inherited methods
#========================================================================>

#------------------------------------------------------------------------
#' Get parameter from NMRData1D object
#'
#' A safe method of getting a parameter value from either acqus or procs
#' list with a fallback to generic parameters list.
#'
#' @name get_parameter 
#' @export
setGeneric("get_parameter", 
           function(object, ...) standardGeneric("get_parameter"))

#------------------------------------------------------------------------
#' @param object An NMRData1D object.
#' @param param.name Name of parameter entry. 
#' @param list.name Name of slot list to acqus (basically either 'acqus' or
#'                  'procs')
#' @param error Issue an error if the required parameters can't be found.
#'
#' @rdname get_parameter
#' @export
setMethod("get_parameter", "NMRData1D", 
  function(object, param.name, list.name, error = FALSE) {

    value <- slot(object, list.name)[[param.name]]

    if ( is.null(value) ) value <- object@parameters[[param.name]]

    msg <- sprintf('Parameter "%s" not found in "%s" or parameters list',
                   param.name, list.name)
    if ( error && is.null(value) ) stop(msg)

    value
  })

#------------------------------------------------------------------------
#' Display NMRData1D object
#'
#' Display a quick summary of the processed data.frame.
#'
#' @name show
#' @export
setMethod("show", "NMRData1D", 
  function(object) {
    cat('An object of NMRData1D class\n')
    print(summary(object@processed))
    cat('\n')
  })

#------------------------------------------------------------------------
#' Set convolution applied to processed data
#'
#' Although it's possible to fit the processed data based on its final
#' appearance, it may sometimes be more useful to try and fit the underlying
#' lineshape with an explicit consideration of the convolution that it has
#' undergone. This convolution includes any apodization as well as the process 
#' of zero-filling, which can be considered as the convolution of a 
#' hypothetically longer experiment by a step function. By default, the 
#' applied apodization is read from the procs parameters, but these can
#' be overwritten by either a list of parameters or a custom function.
#'
#' @param object An NMRData1D object.
#' @param trim Number of point to truncate from the end of the fid's 
#'             acquisition time (not including zero fill). It's common to 
#'             circle shift Bruker's group delay artefact from the start of 
#'             a collected fid to the last 50-100 points. While most 
#'             apodization methods will passively get rid of or diminish the
#'             influence of these points, they may have a stronger impact
#'             with minimal or no apodization.
#' @param param A list of lists combining apodization names and parameters,
#'              (see apodize_signal) for a full list of apodizations and their
#'              parameters. E.g. list('exponential'=list('lb'=1)).
#' @param f A function f(n) capable of generating the convolution vector
#'          of length n (ranging from 0 to the acquisition time). If provided, 
#'          f overrides the param list.
#'
#' @return An NMRData1D object with modified convolution vector. 
#'
#' @name set_convolution
#' @export
setGeneric("set_convolution", 
  function(object, trim = 0, param = list(), f = NULL) {
    standardGeneric("set_convolution")
  })

#' @rdname set_convolution
#' @export
setMethod("set_convolution", "NMRData1D",
  function(object, trim = 0, param = list(), f = NULL) {

    processed <- object@processed

    # Erroring out for unevenly sampled data
    differences <- sort(unique(diff(processed$direct.shift)))
    if ( (length(differences) > 1) && 
         any(abs(diff(differences)) > 1e-10) ) {
      msg <- 'Convolution vectors can only be set for evenly sampled data.'
      stop(msg, call. = TRUE)
    }

    si <- get_parameter(object, 'si', 'procs')
    td <- get_parameter(object, 'td', 'acqus')
    sw <- get_parameter(object, 'sw.h', 'acqus')

    # Correcting group delay
    if ( trim > 0 ) {
      processed <- processed[order(processed$direct.shift), ]
      signal <- ifft(ifftshift(processed$intensity))
      index <- 1:(round(td/2) - trim)
      signal[-index] <- 0
      processed$intensity <- fftshift(fft(signal))

      object@processed <- processed
    }

    # Defining what fraction of overall points are signal (vs. zero-fill)
    frac <- (td/2 - trim)/si
    aq <- 1/(sw)*(td/2 - trim)

    # If param is empty, building up the list from the procs parameters
    if ( length(param) == 0 ) {

      lb <- get_parameter(object, 'lb', 'procs')
      gb <- get_parameter(object, 'gb', 'procs')
      ssb <- get_parameter(object, 'ssb', 'procs')
      tm1 <- get_parameter(object, 'tm1', 'procs')
      tm2 <- get_parameter(object, 'tm2', 'procs')

      if ( lb != 0 ) {
        if ( lb > 0 ) param$exponential <- list(lb = lb)
        else param$gaussian <- list(lb = lb, gb = gb)
      }

      if ( ssb != 0 ) {
        msg <- paste('Assuming that "ssb" parameter pertains to sine rather',
                     'than qsin apodization. Use "param" to override.')
        warning(msg, call. = TRUE)

        param$sine <- list(ssb = ssb)
      }

      if ( (tm1 !=0 ) || (tm2 != 0) ) {
        param$trapezoid <- list(tm1 = tm1, tm2 = tm2)
      }
    }

    # If f is not defined, building f from param list
    if ( is.null(f) ) {

      valid.names <- c('exponential', 'gaussian', 'sine', 'sine2', 'trapezoid')
      if ( any(! names(param) %in% valid.names) ) {
        msg <- sprintf('Names of "param" list must only include: %s',
                       paste(valid.names, collapse = ', '))
        stop(msg, call. = TRUE)
      }

      f <- function(n) {
        
        direct.time <- seq(0, aq, length.out = n)
        signal <- rep(1, n)

        # Looping through the apodizations
        for ( method in names(param) ) {
          apod.args <- list(direct.time = direct.time, fid = signal)
          apod.args <- c(apod.args, method = method, as.list(param[[method]]))

          signal <- do.call(apodize_signal, apod.args)
        }

        signal
      }
    }

    # Filling in values
    n.out <- 2*nrow(processed) - 1
    n.signal <- round(frac*n.out)
    n.fill <- n.out - n.signal
    product <- c(f(n.signal), rep(0, n.fill))
    object@product <- product/sum(abs(product))
    object@convolution <- Re(fftshift(fft(object@product)))

    object
  })

#------------------------------------------------------------------------
#' Filter 1D chemical shift data by selecting specific regions
#'
#' Filters processed data to include only that which is contained between
#' a set of lower and upper bounds on chemical shift.
#'
#' @param object An NMRData1D object.
#' @param lower A lower bound for chemical shift (in the direct dimension).
#' @param upper An upper bound for chemical shift (in the direct dimension).
#' @param round.up True to round up the total number of points to a power of 2.
#'              This is useful to ensure rapid convolution using the fft when
#'              convolution is included in the fit.
#'
#' @return An NMRData1D object with filtered processed data.
#'
#' @name filter_1d
#' @export
setGeneric("filter_1d", 
  function(object, lower = NULL, upper = NULL, round.up = FALSE, ...) {
    standardGeneric("filter_1d")
  })

#' @rdname filter_1d 
#' @export
setMethod("filter_1d", "NMRData1D", 
  function(object, lower, upper, round.up) {

    x <- object@processed$direct.shift
  
    if ( is.null(lower) ) lower <- min(x)
    if ( is.null(upper) ) upper <- max(x)

    # Checking bounds
    if ( lower > upper ) {
      msg <- 'The lower bound must be smaller than or equal to the upper bound.'
      stop(msg)
    }

    if ( length(c(lower, upper)) != 2 ) {
      msg <- 'Multiple range definition is not currently supported.'
      stop(msg)
    }

    # If x happens to be empty, don't do anything else
    if ( length(x) == 0 ) return(object)

    # Warning if bounds are outside range of chemical shifts 
    if ( (lower > max(x)) | (upper < min(x)) ) {
      msg <- 'Some bounds are entirely out of chemical shift range'
      warning(msg, call. = FALSE)
    }

    # Generating index
    index <- (x > lower) & (x < upper)

    # Rounding if necessary
    if ( round.up ) {
      n <- sum(index)
      pwr = log(n)/log(2)

      n.new <- 2^ceiling(pwr)
      n.extra <- n.new - n

      # Figuring out if there is space to add new points equally on either side
      rle.index <- rle(index)


      # If there are three entries, they are FALSE, TRUE, FALSE
      if ( length(rle.index$value) == 3 ) {
        left.count <- rle.index$lengths[1]
        i.left <- left.count + 1

        right.count <- rle.index$lengths[3]
        i.right <- sum(rle.index$lengths[1:2])
        
        n.left <- min(left.count, floor(n.extra/2))
        index[(i.left - n.left):(i.left - 1)] <- TRUE
        n.extra <- n.extra - n.left

        n.right <- min(right.count, n.extra)
        index[(i.right + 1):(i.right + n.right)] <- TRUE
        n.extra <- n.extra - n.right
      # If there are two entries, they are TRUE FALSE, or FALSE TRUE
      } else if ( length(rle.index$value) == 2 ) {
        if ( rle.index$value[1] ) {
          right.count <- rle.index$lengths[2]
          i.right <- rle.index$lengths[1:2]
          n.right <- min(right.count, n.extra)
          index[(i.right + 1):(i.right + n.right)] <- TRUE
          n.extra <- n.extra - n.right
        } else {
          left.count <- rle.index$lengths[1]
          i.left <- left.count + 1
          n.left <- min(left.count, n.extra)
          index[(i.left - n.left):(i.left - 1)] <- TRUE
          n.extra <- n.extra - n.left
        }
      } 
    
      # Checking if there are any more indexes unaccounted for
      if ( n.extra > 0 ) {
        msg <- 'It was not possible to round desired range to a power of 2.'
        warning(msg, call. = FALSE)
      }
    
    }

    object@processed <- object@processed[index, ]

    # Adjusting length of convolution vectors
    if ( length(object@product) > 0 ) {
      n <- 2*sum(index) - 1
      product <- object@product
      object@product <- approx(1:length(product), product, n = n)$y
      object@product <- object@product/sum(object@product)
      object@convolution <- Re(fftshift(fft(object@product)))
    }

    object
  })

#========================================================================>
# Defining a plot function
#========================================================================>

#------------------------------------------------------------------------
#' Plot NMRData1D object
#'
#' Generates an interactive plot object using the plotly package.
#'
#' @param x NMRData1D object.
#' @param components One of either 'r', 'i', or 'r/i' to include real,
#'                   imaginary or both components. If both components
#'                   are selected, they are displayed in separate subplots.
#' @param nrows Max number of rows to display subplots.
#' @param reverse TRUE to order x-axis with large values on the left and 
#'                small values on the right like typical NMR plots.
#'
#' @export
plot.NMRData1D <- function(x, components = 'r', nrows = 2, reverse = TRUE) {
  d <- processed(x)

  plots <- list()
  legend.opts <- list(orientation = 'h', xanchor = "center", x = 0.5)

  if ( grepl('r', components) ) {
    p <- plot_ly(x = d$direct.shift, y = Re(d$intensity), 
                 color = I('black'), name = 'Real') %>%
           add_lines()  %>%
           layout(legend = legend.opts)

    if ( reverse ) p <- p %>% layout(xaxis = list(autorange = "reversed"))

    plots$r <- p
  }

  if ( grepl('i', components) ) {
    p <- plot_ly(x = d$direct.shift, y = Im(d$intensity), 
                 color = I('grey'), name = 'Imaginary') %>%
           add_lines()  %>%
           layout(legend = legend.opts)

    if ( reverse ) p <- p %>% layout(xaxis = list(autorange = "reversed"))

    plots$i <- p
  }

  if ( length(plots) == 0 ) NULL
  else subplot(plots, shareX = TRUE, shareY = TRUE, 
               nrows = min(length(plots), nrows))
}

setMethod("plot", "NMRData1D", plot.NMRData1D)
