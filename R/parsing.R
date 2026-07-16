# Parsing coupling patterns

#------------------------------------------------------------------------
#' A convenience function for coupling syntax error messages.
#' 
check_syntax <- function(invalid.syntax, syntax, hint) {
  if ( any(invalid.syntax) ) {
     msg <- sprintf('The following multiplets could not be parsed (%s): %s',
                    hint, paste(syntax[invalid.syntax], collapse = ', '))
     stop(msg, call. = FALSE)
   }
}
#------------------------------------------------------------------------
#' Parse coupling strings
#'
#' Converts a coupling specification of the form '3.0 d 1.0' into the
#' chemical shift, thenumber of peaks involved and the coupling between them. 
#' Currently, the supported codes include s, d, t, q, quint, sext, sept, oct, 
#' dd, dt, td, and tt. The latter 4 expect two coupling values. The exact 
#' syntax of the specification is quite flexible, as the function tries to 
#' remove all brackets and punctuation and convert common spacers like - and _
#' into whitespace.
#'
#' @param coupling.string A character vector with elements of the form 
#'                        '3.0 d 1.0'
#'
#' @return A list of three lists -- chemical.shift, peaks, and coupling
#'         respectively. Both peaks and couplng lists may have one or two
#'         length vectors as elements. Singlets are assigned a coupling of NA.
#'
#' @export
parse_coupling <- function(coupling.string) {

  # Preserving multiplet names and original values
  ids <- names(coupling.string)
  blank <- ids == ''
  
  # Numbering blank ids
  if ( is.null(ids) ) {
    ids <- 1:length(coupling.string)
  } else {
    numbers <- 1:length(coupling.string)
    ids[blank] <- numbers[blank]
  }

  ids <- factor(ids, levels = ids)
  original <- coupling.string 

  # Dealing with numeric elements, allowing for peak.list to be a list or vector
  coupling.string <- lapply(coupling.string, 
                     function(x) ifelse(is.numeric(x), paste(x, 's'), x))
  coupling.string <- unlist(coupling.string)

  # Ensure lowercase
  coupling.string <- tolower(coupling.string)

  # Remove any and all brackets
  coupling.string <- str_replace_all(coupling.string, '[(){}]|\\[\\]', ' ')

  # Replace all possible separators with a single whitespace
  coupling.string <- str_replace_all(coupling.string, '[ ,;#$%_]+', ' ')

  # Remove whitespace at beginning or end
  coupling.string <- str_trim(coupling.string, 'both')

  #------------------------------------------------------------
  # Dealing with satellite specification

  logic <- grepl('j', coupling.string) 
  if ( any(logic) ) {
    sat.mat <- str_split(coupling.string, 'j', n = 2, simplify = TRUE)

    # The left portion is processed as normal
    coupling.string <- str_trim(sat.mat[, 1], 'both')

    # The right portion is parsed as satellites
    sat.string <- str_trim(sat.mat[logic, 2], 'both')
    sat.ids <- ids[logic]
    levels(sat.ids) <- levels(ids)

    # Splitting by colon to identify offsets
    codes <- str_split(sat.string, ' ')
    sat.numbers <- unlist(lapply(codes, length))
    names(sat.numbers) <- sat.ids

    # Flattening to parse colons
    flat.codes <- unlist(codes)
    flat.mat <- str_split(flat.codes, ':', n = 2, simplify = TRUE)

    # Preparing to reset proper lengths
    split.groups <- rep(1:length(sat.numbers), sat.numbers)

    flat.coupling <- as.numeric(flat.mat[, 1])
    flat.offsets <- flat.mat[, 2]
    flat.offsets <- as.numeric(ifelse(flat.offsets == '', '0', flat.offsets))

    flat.invalid <- is.na(flat.coupling) | is.na(flat.offsets)

    sat.invalid <- unlist(lapply(split(flat.invalid, split.groups), any))
    invalid.codes <- rep(FALSE, length(original))
    invalid.codes[logic] <- sat.invalid
    
    check_syntax(invalid.codes, original,
                 'likely invalid satellite specification') 

    sat.offsets <- split(flat.offsets, split.groups)
    names(sat.offsets) <- sat.ids
    sat.coupling <- split(flat.coupling, split.groups)
    names(sat.coupling) <- sat.ids

  } else {
    sat.numbers <- NULL
    sat.offsets <- NULL
    sat.coupling <- NULL
    sat.ids <- NULL
  }

  #------------------------------------------------------------
  # Moving on to multiplet codes

  # Split into three columns
  coupling.mat <- str_split(coupling.string, ' ', n = 3, simplify = TRUE)

  # Converting shift to numeric
  shifts <- as.numeric(coupling.mat[, 1])
  
  check_syntax(is.na(shifts), original, 
               'likely invalid chemical shift')

  # Checking codes for overall validity
  codes <- coupling.mat[, 2]

  s <- '^s$'
  m1 <- '^d$|^t$|^q$|^(quint)$|^(sext)$|^(sept)$|^(oct)$'
  m2 <- '^(dd)$|^(dt)$|^(td)$|^(tt)$'

  is.s <- grepl(s, codes)
  is.m1 <- grepl(m1, codes)
  is.m2 <- grepl(m2, codes)

  invalid.codes <- ! ( is.s | is.m1 | is.m2 )
  check_syntax(invalid.codes, original, 
               'likely invalid multiplet syntax') 

  # Converting codes to numbers
  peak.numbers <- list('s' = 1, 'd' = 2, 't' = 3, 'q' = 4, 
                       'quint' = 5, 'sext' = 6, 'sept' = 7, 'oct' = 8,
                       'dd' = c(2, 2), 'dt' = c(2, 3), 
                       'td' = c(3, 2), 'tt' = c(3, 3))
  peak.numbers <- peak.numbers[codes]
  names(peak.numbers) <- NULL

  # Converting couplings to numbers
  coupling <- lapply(str_split(coupling.mat[, 3], ' '), as.numeric)
  coupling.length <- lapply(coupling, length)

  # Checking coupling lengths
  invalid.couplings <- rep(FALSE, length(coupling))

  invalid.couplings[is.s] <- ! is.na(coupling[is.s])
  invalid.couplings[is.m1] <- coupling.length[is.m1] != 1
  invalid.couplings[is.m2] <- coupling.length[is.m2] != 2

  check_syntax(invalid.couplings, original,
              'likely invalid coupling specification') 

  list(ids = ids, chemical.shift = shifts, 
       peaks = peak.numbers, coupling = coupling,
       sat.ids = sat.ids, sat.offsets = sat.offsets,
       sat.numbers = sat.numbers, sat.coupling = sat.coupling)
}

#------------------------------------------------------------------------
#' Generate splitting pattern based on coupling numbers and constants. 
#'
#' Converts a vector of peak numbers and coupling constants into a pattern
#' of single peaks with accompanying area/height ratios as well as chemical
#' shift offsets. An input of c(2, 2) and c(1, 2) would correspond to a
#' double doublet with coupling constants of 1 and 2 respectively. Although
#' units aren't kept track of, it makes most sense to convert coupling
#' constants to ppm.
#'
#' @param peaks A vector of peaks per coupling constant, e.g., c(2, 2) for dd.
#' @param coupling A vector of coupling constants.
#'
#; @return A data.frame of chemical shift offsets and height/area ratios in
#'         relation to the center.
#' 
#' @export
split_coupling <- function(peaks, coupling) {

  if (length(peaks) != length(coupling)) {
    msg <- 'Vector of peaks must match vector of coupling constants.'
    stop(msg)
  }

  # Handling the trivial case
  if ( (length(peaks) == 1) && (peaks == 1) ) {
    out <- data.frame(offsets = 0, ratios = 1)
    return(out)
  }

  offsets <- list()

  for ( i in 1:length(peaks) ) {
    if ( peaks[1] == 1 ) {
      offsets[[i]] <- 1
    } else if ( (peaks[i] %% 2) == 0 ) {
      offsets[[i]] <- c(-0.5, 0.5)
      n.remainder <- (peaks[i] - 2)/2
      i.remainder <- seq(1.5, 0.5 + n.remainder, length.out = n.remainder)
      offsets[[i]] <- c(-rev(i.remainder), offsets[[i]], i.remainder)
    } else {
      offsets[[i]] <- c(0)
      n.remainder <- (peaks[i] - 1)/2
      i.remainder <- seq(1, n.remainder, length.out = n.remainder)
      offsets[[i]] <- c(-rev(i.remainder), offsets[[i]], i.remainder)
    }
  }

  # Expanding and multiplying
  offsets <- expand.grid(rev(offsets))
  offsets <- rowSums(offsets * rep(rev(coupling), each = nrow(offsets)))

  # Iteratively calculate ratios. First, calculate the ratio associated with
  # the last coupling, then repeat the process for the first if there are
  # more than one. The following assumes that n can only be 1 or 2
  n <- length(peaks)
  ratios <- choose(peaks[n]-1, 0:(peaks[n]-1))
  if (n > 1) {
    ratios <- rep(ratios, prod(peaks[-n]))
    new.ratio <- choose(peaks[1]-1, 0:(peaks[1]-1))
    ratios <- ratios * rep(new.ratio, each = peaks[2])
  }

  # Combining and adding up peaks that fall on each other
  out <- cbind(offsets, ratios)
  out <- out[order(offsets), ]

  r <- rle(offsets)
  logic <- r$lengths > 1
  n <- sum(logic)

  i.end <- cumsum(r$lengths)[logic]
  i.start <- i.end - r$lengths[logic] + 1

  for (i in seq(n, 1, length.out = n)) {
    out[i.start[i], 2] <- sum(out[i.start[i]:i.end[i], 2])
    out <- out[-((i.start[i]+1):i.end[i]), ]
  }

  out
}
