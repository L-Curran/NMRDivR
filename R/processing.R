# Functions for processing spectral data

#========================================================================>
# Misc. fft functions

ifft <- function(x) {
  fft(x, inverse = TRUE)/length(x)
}

fftshift <- function(x) {
  n <- length(x)
  i <- ceiling(n/2)

  x[c((i + 1):n, 1:i)]
}

ifftshift <- function(x) {
  n <- length(x)
  i <- floor(n/2)

  x[c((i + 1):n, 1:i)]
}

#========================================================================>
# Window/weighing functions

#------------------------------------------------------------------------
#' Apply apodization function
#'
#' A thin wrapper around a number of common apodization functions.
#'
#' @param direct.time A vector of time points of the signal data.
#' @param signal A vector of complex signal data.
#' @param method One of either 'exponential', 'gaussian', 'sine', 'sine2',
#'               or 'trapezoid' (where sine2 stands for sine squared -- qsin).
#' @param output.weights Set to TRUE to output the apodization function rather
#'                       than the weighted signal.
#' @param ... Shape parameters that depend on method,
#'              lb -- Line broadening (in Hz) used for exponential and
#'                    gaussian apodization (negative for gaussian)
#'              gb -- Gaussian broadening factor (0 - 1)
#'              ssb -- Sine shift bell used for sine and sine2 apodization
#'              tm1 -- Leftmost edge of trapezoid (0 - 1)
#'              tm2 -- Rightmost edge of trapezoid (0 - 1)
#'
#' @return Corrected signal vector (or vector of weights if output.weights
#'         is TRUE).
#'
#' @export
apodize_signal <- function(direct.time, fid, method = 'exponential', 
                           ..., output.weights = FALSE) {

  if (is.null(method)) {
    return(fid)
  }

  functions <- list(apodize_exponential, apodize_gaussian, 
                      apodize_sine, apodize_sine2, apodize_trapezoid)
  names(functions) <- c('exponential', 'gaussian', 'sine', 'sine2', 'trapezoid')

  if (! method %in% names(functions) ) {
    msg <- sprintf('"method" must be one of %s',
                   paste(names(functions), sep = ', '))
    stop(msg)
  } else {
    f <- functions[[method]]
  }

  f(direct.time, fid, output.weights = output.weights, ...)
}

#------------------------------------------------------------------------
#' Apply exponential apodization function
#'
#' Apodizes signal data according to an exponential function.
#'
#' @param direct.time A vector of time points of the signal data.
#' @param signal A vector of complex signal data.
#' @param lb Line broadening factor (in Hz), where the apodization function
#'           is calculated as exp(-lb*pi*t).
#' @param output.weights Set to TRUE to output the apodization function rather
#'                       than the weighted signal.
#'
#' @return Corrected signal vector (or vector of weights if output.weights
#'         is TRUE).
#'
#' @export
apodize_exponential <- function(direct.time, signal, lb = 1,
                                output.weights = FALSE) {

  weight <- exp(-lb*pi*direct.time)
  weight <- weight/max(weight)

  if ( output.weights ) {
    weight
  } else {
    signal*weight
  }
}

#------------------------------------------------------------------------
#' Apply Gaussian apodization function
#'
#' Apodizes signal data according to a Gaussian function.
#'
#' @param direct.time A vector of time points of the signal data.
#' @param signal A vector of complex signal data.
#' @param lb Line broadening factor (in Hz) component of the apodization (must
#'           be negative to conform with Topspin parameters).
#' @param gb Gaussian broadening factor (0 - 1), corresponding to the position of
#'           the function maximum as a fraction of acquisition time.
#' @param output.weights Set to TRUE to output the apodization function rather
#'                       than the weighted signal.
#'
#' @return Corrected signal vector (or vector of weights if output.weights
#'         is TRUE).
#'
#' @export
apodize_gaussian <- function(direct.time, signal, lb = 1, gb = 0.5,
                             output.weights = FALSE) {
  
  a <- -pi*lb
  b <- a/(2*gb*max(direct.time))
  
  weight <- exp(-a*direct.time + b*direct.time^2)
  weight <- weight/max(weight)
  
  if ( output.weights ) {
    weight
  } else {
    signal*weight
  }
}

#------------------------------------------------------------------------
#' Apply sine apodization function
#'
#' Apodizes signal data according to a sine function.
#'
#' @param direct.time A vector of time points of the signal data.
#' @param signal A vector of complex signal data.
#' @param ssb Sine bell position parameter -- 1 for pure sine, 2 for pure cosine.
#' @param output.weights Set to TRUE to output the apodization function rather
#'                       than the weighted signal.
#'
#' @return Corrected signal vector (or vector of weights if output.weights
#'         is TRUE).
#'
#' @export
apodize_sine <- function(direct.time, signal, ssb = 1,
                         output.weights = FALSE) {

  if ( ssb == 1 ) {
    weight <- sin(pi*direct.time/max(direct.time))
  } else if ( ssb == 2 ) {
    weight <- cos(pi/2*direct.time/max(direct.time))
  } else {
    a <- pi/ssb
    weight <- sin((pi - a)*(direct.time/max(direct.time)) + a)
  }
  
  weight <- weight/max(weight)

  if ( output.weights ) {
    weight
  } else {
    signal*weight
  }
}

#------------------------------------------------------------------------
#' Apply sine squared apodization function
#'
#' Apodizes signal data according to a squared sine function.
#'
#' @param direct.time A vector of time points of the signal data.
#' @param signal A vector of complex signal data.
#' @param ssb Sine bell position parameter -- 1 for pure sine, 2 for pure cosine.
#' @param output.weights Set to TRUE to output the apodization function rather
#'                       than the weighted signal.
#'
#' @return Corrected signal vector (or vector of weights if output.weights
#'         is TRUE).
#'
#' @export
apodize_sine2 <- function(direct.time, signal, ssb = 1,
                          output.weights = FALSE) {
  
  if ( ssb < 2 ) {
    weight <- (sin(pi*x/aq))^2
  } else {
    a <- pi/ssb
    weight <- (sin((pi - a)*(direct.time/max(direct.time)) + a))^2
  }
  
  weight <- weight/max(weight)
  
  if ( output.weights ) {
    weight
  } else {
    signal*weight
  }
}

#------------------------------------------------------------------------
#' Apply trapezoid apodization function
#'
#' Apodizes signal data according to a trapezoid function that ramps up from
#' 0 to tm1 and ramps down from tm2 until the end of the acquisition time, where
#' tm1 and tm2 are relative values from 0 to 1.
#'
#' @param direct.time A vector of time points of the signal data.
#' @param signal A vector of complex signal data.
#' @param tm1 The leftmost edge of the trapezoid as a fraction of 
#'            acquisition time.
#' @param tm2 The rightmost edge of the trapezoid as a fraction of 
#'            acquisition time.
#' @param output.weights Set to TRUE to output the apodization function rather
#'                       than the weighted signal.
#'
#' @return Corrected signal vector (or vector of weights if output.weights
#'         is TRUE).
#'
#' @export
apodize_trapezoid <- function(direct.time, signal, tm1 = 0, tm2 = 0.5,
                              output.weights = FALSE) {

  x <- direct.time/max(direct.time)
  
  # Applying the ramp down
  weight <- ifelse(x > tm2, (tm2-x)/(1-tm2) + 1, 1)
  
  # Ramp up if needed
  if (tm1 > 0) weight <- ifelse((x < tm1) & (x < tm2), x/tm1, weight)
  
  weight <- weight/max(weight)
  
  if ( output.weights ) {
    weight
  } else {
    signal*weight
  }
}

#========================================================================>
# Zero-fill

#------------------------------------------------------------------------
#' Zero-fill signal
#'
#' Extends the signal by a specified number of time by adding zeros at
#' the end.
#'
#' @param direct.time A vector of time points of the signal data.
#' @param signal A vector of complex signal data.
#' @param times The number of times to multiply original data length.
#'
#' @return Corrected signal vector.
#'
#' @export
zero_fill <- function(direct.time, fid, times = 1) {

  n <- length(direct.time)
  min.time <- direct.time[1]
  max.time <- direct.time[n]*(1 + times)
  delta <- direct.time[2] - direct.time[1]

  new.time <- seq(min.time, max.time, delta)
  new.n <- length(new.time) - n

  data.frame(direct.time = new.time,
             signal = c(fid, rep(complex(re = 0, im = 0), new.n)))
}


#========================================================================>
# Phasing and shifting

#------------------------------------------------------------------------
#' Decode echo-antiecho signals
#'
#' Decodes echo-antiecho pulses by converting them into real/imaginary
#' sequence pairs.
#'
#' @param indirect.time A vector of delay times for the signal.
#' @param direct.time A vector of sampling times for the signal.
#' @param signal A vector of complex fid data.
#'
#' @return Corrected fid data.
#'
#' @export
corr_antiecho <- function(indirect.time, direct.time, signal) {

  # Generating data.frame
  d <- data.frame(indirect.time, direct.time, signal)

  # Ensuring proper order
  d <- arrange(d, indirect.time, direct.time)

  # Counting
  d.n <- d %>%
           group_by(indirect.time) %>%
           summarize(n.signal = n())
  d.n <- d.n$n.signal

  # Adding pairs
  n <- length(d.n)
  n.pairs <- d.n[seq(1, n, by = 2)] + d.n[seq(2, n, by = 2)]
  i.pairs <- rep(seq(1, n/2), n.pairs)

  d$i.pair <- i.pairs

  # Defining decoding function
  f_decode <- function(x, i) {
    n <- n.pairs[i]
    
    x1 <- x[1:(n/2)]
    x2 <- x[(n/2 + 1):n]

    x[1:(n/2)] <- (x1 + x2)/2
    x[(n/2 + 1):n] <- (x1 - x2)/complex(imaginary = 2)

    x
  }

  # Performing the decoding
  d <- d %>%
         group_by(i.pair) %>%
         mutate(signal = f_decode(signal, i.pair[1]))

  d$signal
}


#------------------------------------------------------------------------
#' Correct Bruker group delay artefact
#'
#' Corrects Bruker group delay artefact based on process of Westler and 
#' Abildgaard (1996): 
#' https://ucdb.googlecode.com/hg/application/ProSpectND/html/
#'         dmx_digital_filters.html
#'
#' See the following for details about artefact:
#' http://nmr-analysis.blogspot.ca/2008/02/
#'        why-arent-bruker-fids-time-corrected.html
#'
#' @param signal A vector of complex fid data.
#' @param acqus.list A list of acqus parameters that contains 'grpdly'
#'                   or 'dspfvs' and 'decim' entries. This list can be 
#'                   generated using read_acqus_1d() or through other means.
#'                   These parameters can also be nested within a list
#'                   item called 'acqus' if multiple dimensions are read
#'                   at once.
#'
#' @return Corrected fid data.
#'
#' @export
corr_group_delay <- function(signal, acqus.list) {

  # Checking for required entries
  required.acqus <- c('grpdly', 'dspfvs', 'decim')
  missing <- ! required.acqus %in% names(acqus.list)

  # If all parameters are missing, check to see if the acqus parameters
  # are neseted within acqus.list.
  if (all(missing)) {
    acqus.list <- acqus.list$acqus
  }

  if (!'grpdly' %in% names(acqus.list)) {
    required.acqus <- c('dspfvs', 'decim')
    missing <- ! required.acqus %in% names(acqus.list)

    if (any(missing)) {
      msg <- sprintf('Either %s or one of %s must be provided in acqus.list',
                     'grpdly', paste(required.acqus, collapse=', '))
      stop(msg)
    }
  }

  # Correction lookup table (unrounded table from Python's nmrpipe package)
  correction_table <- list(
    '10'=list('2'=44.75, '3'=33.5, '4'=66.625, '6'=59.083333333333333, 
              '8'=68.5625, '12'=60.375, '16'=69.53125, '24'=61.020833333333333, 
              '32'=70.015625, '48'=61.34375, '64'=70.2578125, 
              '96'=61.505208333333333, '128'=70.37890625, '192'=61.5859375, 
              '256'=70.439453125, '384'=61.626302083333333, '512'=70.4697265625,
              '768'=61.646484375, '1024'=70.48486328125, 
              '1536'=61.656575520833333, '2048'=70.492431640625),
    '11'=list('2'=46., '3'=36.5, '4'=48., '6'=50.166666666666667, '8'=53.25,
              '12'=69.5, '16'=72.25, '24'=70.166666666666667, '32'=72.75, 
              '48'=70.5, '64'=73., '96'=70.666666666666667, '128'=72.5, 
              '192'=71.333333333333333, '256'=72.25, '384'=71.666666666666667, 
              '512'=72.125, '768'=71.833333333333333, '1024'=72.0625, 
              '1536'=71.916666666666667, '2048'=72.03125),
    '12'=list('2'=46., '3'=36.5, '4'=48., '6'=50.166666666666667, '8'=53.25, 
              '12'=69.5, '16'=71.625, '24'=70.166666666666667, '32'=72.125, 
              '48'=70.5, '64'=72.375, '96'=70.666666666666667, '128'=72.5, 
              '192'=71.333333333333333, '256'=72.25, '384'=71.666666666666667, 
              '512'=72.125, '768'=71.833333333333333, '1024'=72.0625, 
              '1536'=71.916666666666667, '2048'=72.03125),
    '13'=list('2'=2.75, '3'=2.8333333333333333, '4'=2.875, 
              '6'=2.9166666666666667, '8'=2.9375, '12'=2.9583333333333333, 
              '16'=2.96875, '24'=2.9791666666666667, '32'=2.984375, 
              '48'=2.9895833333333333, '64'=2.9921875, '96'=2.9947916666666667)) 

  # Reading off phase correction
  if ('grpdly' %in% names(acqus.list)) {
    total.phase <- acqus.list[['grpdly']]
  } else {
    dspfvs <- as.character(acqus.list[['dspfvs']])
    decim <- as.character(acqus.list[['decim']])

    total.phase <- correction_table[[dspfvs]][[decim]]
  }

  # Circle shift by whole number
  truncated <- floor(total.phase)
  n <- length(signal)
  shift <- exp(2*pi*seq(0,n-1)/n * truncated * complex(imaginary=1))
  signal <- ifft(fft(signal)*shift)

  # Touching up with phase correction
  remainder <- total.phase - truncated
  phase_signal(signal, c(0, remainder*360))
}

#------------------------------------------------------------------------
#' Correct signal phase
#'
#' This function applies a phase correction directly to the signal in the 
#' time domain once the correct phase angles have been identified.
#'
#' @param signal A vector of complex signal data.
#' @param phase A vector of phase corrections, from zero-order and up.
#' @param degrees TRUE if the phase angles are in degrees, FALSE if in radians.
#
#' @return Corrected signal data.
#'
#' @export
phase_signal <- function(signal, phase, degrees=TRUE) {

  # Changing phase angle units if required
  if (degrees) phase <- phase*pi/180

  # Summing up phase components
  n.phase <- length(phase)
  if ( n.phase == 0 ) return(signal)

  n <- length(signal)
  phase.total <- rep(phase[1], n)

  xf <- seq(0, n-1, n)
  x <- rep(1, n)

  if ( n.phase > 1 ) {
    x <- x*xf

    for (i in 2:n.phase) {
      phase.total <- phase.total + phase[i]*x
    }
  }

  # Correcting
  apod <- exp(phase.total * complex(imaginary=1))
  signal <- signal*apod
}

#------------------------------------------------------------------------
#' Correct spectrum phase
#'
#' This function applies a phase correction to the fourier transformed
#' spectrum once the correct phase angles have been identified.
#'
#' @param spectrum A vector of complex spectrum data.
#' @param phase A vector of phase corrections, from zero-order and up.
#' @param degrees TRUE if the phase angles are in degrees, FALSE if in radians.
#
#' @return Corrected spectrum.
#'
#' @export
phase_spectrum <- function(spectrum, phase, degrees=TRUE) {

  # Summing up phase components
  n.phase <- length(phase)
  if ( n.phase == 0 ) return(spectrum)

  # Changing phase angle units if required
  if (degrees) phase <- phase*pi/180

  n <- length(spectrum)
  phase.total <- rep(phase[1], n)

  xf <- seq(0, n-1, n)
  x <- rep(1, n)

  if ( n.phase > 1 ) {
    x <- x*xf

    for (i in 2:n.phase) {
      phase.total <- phase.total + phase[i]*x
    }
  }  

  im <- Im(spectrum)
  re <- Re(spectrum)

  complex(real = re * cos(phase.total) + im * sin(phase.total),
          imaginary = im * cos(phase.total) - re * sin(phase.total))
}

#========================================================================>
# Signal to spectrum

#------------------------------------------------------------------------
#' Fourier tranform signal data 
#'
#' Performs fourier transform of complex signal data (i.e., in the direct
#' dimension).
#'
#' @param direct.time A vector of sampling times for the signal.
#' @param signal A vector of complex signal data.
#' @param acqus.list A list of acqus parameters that contains 'o1', 'sfo1', 
#'                   and 'sw' entries. This list can be 
#'                   generated using read_acqus_1d() or through other means.
#'                   These parameters can also be nested within a list
#'                   item called 'acqus' if multiple dimensions are read
#'                   at once.
#
#' @return A data.frame made of two columns -- "direct.shift" containing the 
#'         chemical shift and "intensity" containing the real valued 
#'         spectrum data from the real component of the signal.
#'
#' @export
fft_signal <- function(direct.time, signal, acqus.list) {

  # Checking for required acqus.entries
  required.acqus <- c('o1', 'sfo1', 'sw')
  missing <- ! required.acqus %in% names(acqus.list)

  # If all parameters are missing, check to see if the acqus parameters
  # are nested within acqus.list.
  if (all(missing)) {
    acqus.list <- acqus.list$acqus
    missing <- ! required.acqus %in% names(acqus.list)
  }

  # If any parameters are still missing, issue error
  if (any(missing)) {
    msg <- sprintf('The following parameters are missing from acqus.list: %s',
                   paste(required.acqus[missing], collapse=', '))
    stop(msg, call.=FALSE)
  }

  # Extracting parameters
  o1 <- acqus.list$o1
  sw <- acqus.list$sw
  sfo1 <- acqus.list$sfo1

  n <- length(signal)
  center <- o1/sfo1

  direct.shift <- seq(center - sw/2, center + sw/2, length.out = n)
  intensity <- fftshift(fft(signal))

  data.frame(direct.shift = direct.shift, intensity = intensity)
}

#------------------------------------------------------------------------
#' Fourier tranform intensity data 
#'
#' Performs fourier transform of complex intensity data (i.e., in the 
#' indirect dimension). Although the procedure is very similar to
#' fft_signal, real/imaginary data processing will differ depending on
#' the mode of indirect data acquisition.
#'
#' @param indirect.time A vector of delay times for the signal.
#' @param intensity A vector of complex intensity data.
#' @param acqus.list A list of lists containing acqus parameters with 
#'                   'o1', 'sfo1', 'sw', and 'fnmode' entries. 
#'                   This list can be generated using read_acqus() or 
#'                   through other means. These parameters can also be nested 
#'                   within a list item called 'acqus2' if multiple dimensions
#'                   are read at once.
#' @param hypercomplex TRUE to output full quadrature components
#'                     (rr, ri, ir, ii), FALSE to omit imaginary components
#'                     in the direct dimension (ir, ii). 
#
#' @return A data.frame made of two columns -- "indirect.shift" containing
#'         the indirect dimension chemical shift and "intensity" containing 
#'         the real component of the fourier transformed direct dimension data 
#'         (corresponding to rr and ri as a complex pair). If the hypercomplex 
#'         option is set to TRUE, the ir and ii data is added as another
#'         complex column labelled "dispersive".
#'
#' @export
fft_spectrum <- function(indirect.time, intensity, acqus.list, 
                         hypercomplex = FALSE) {
  
  # Checking for required acqus.entries
  required.acqus <- c('o1', 'sfo1', 'sw', 'fnmode')
  missing <- ! required.acqus %in% names(acqus.list)

  # If all parameters are missing, check to see if the acqus parameters
  # are nested within acqus.list.
  if (all(missing)) {
    acqus.list <- acqus.list$acqu2s
    missing <- ! required.acqus %in% names(acqus.list)
  }

  # If any parameters are still missing, issue error
  if (any(missing)) {
    msg <- sprintf('The following parameters are missing from acqus.list: %s',
                   paste(required.acqus[missing], collapse=', '))
    stop(msg, call.=FALSE)
  }

  # Extracting parameters
  o1 <- acqus.list$o1
  sw <- acqus.list$sw
  sfo1 <- acqus.list$sfo1
  fnmode <- acqus.list$fnmode

  # The specific processing method depends on fnmode
  if ( (fnmode == 1) || (fnmode == 'qf') ) {

    n <- length(intensity)
    center <- o1/sfo1

    indirect.shift <- seq(center + sw/2, center - sw/2, length.out = n)

    r <- fftshift(fft(intensity))
    i <- NA

  } else if ( (fnmode == 6) || (fnmode == 'echo-antiecho') ) {

    n <- length(intensity)/2
    center <- o1/sfo1

    indirect.shift <- seq(center + sw/2, center - sw/2, length.out = n)

    odd <- seq(1, n*2, by=2)
    even <- seq(2, n*2, by=2)

    r <- complex(real = Re(intensity[odd]), 
                 imaginary = Re(intensity[even]))
    r <- fftshift(fft(r))

    i <- complex(real = Im(intensity[odd]), 
                 imaginary = Im(intensity[even]))
    i <- fftshift(fft(i))
  } else {
    msg <- '"%s" fnmode not implemented'
    stop(msg)
  }

  d <- data.frame(indirect.shift = indirect.shift, intensity = r)

  if ( hypercomplex ) d$dispersive <- i
  d
}
