# -------------------------------------------------------------------------
# Vendored from: rnmrfit (version 1.0.0)
# Original Author/Copyright: Stanislav Sokolenk", email = stanislav@sokolenko.net, author/creator
# License: Apache (>= 2.0)
#
# Description: Support functions to import and read spectra from bruker format

# Functions for importing and processing jcamp data

#========================================================================>
# JCAMP import functions

#------------------------------------------------------------------------
#' Import data from generic JCAMP file
#'
#' Reads data from JCAMP file in the most general fashion possible. No
#' checks are performed on file validity. Tag and entry conversion is
#' done by passing output directly into process_jcamp().
#'
#' @param path Path to file.
#' @param process.tags TRUE to convert entry tags to lowercase with dots
#'                     in place of spaces.
#' @param process.entries FALSE to keep all block elements as characters.
#'                        TRUE to attempt some basic conversion.
#' @param ... Arguments passed into process_jcamp().
#
#' @return A list with the following entries: 'header', 'blocks',
#'         and 'comments'. Header elements are all those contained between
#'         the 'TITLE', 'OWNER', and 'ORIGIN', and may vary between files.
#'         All bon-header entries (including notes) are stored
#'         in the 'blocks' sublist. If no block information is specified,
#'         all information is stored in block 1. comments' contains a list of
#'         comment blocks indexed by corresponding JCAMP element of the line
#'         number on which they were found (if there is no corresponding
#'         element). If either convert.tags or contvert.entries is true,
#'         the parsed output is also processed using the process_jcamp()
#'         function. See ?process_jcamp for more details.
#'
#'
read_jcamp <- function(path, process.tags = TRUE,
                       process.entries = TRUE, ...) {

  # Reading file
  file.string <- safe_read(path, 'char')

  #------------------------------------------------------------------------
  # First, splitting by line
  by.line <- str_split(file.string, '[ \t\r]*\n[ \t]*')[[1]]

  # Removing lines that look like breaks or comments
  remove <- grepl('^[.*=-]{3,}', by.line)
  remove <- remove | (by.line == '')

  by.line <- by.line[!remove]

  # Initializing comments
  comments <- list()

  # Parsing all comments that start on blank line
  line.comments <- str_extract(by.line, '(?<=^((\\$\\$)|(##=))).*')
  line.numbers <- which(!is.na(line.comments))

  comments <- as.list(str_trim(line.comments[line.numbers], 'both'))
  names(comments) <- line.numbers

  # Recombining lines without the comments
  if (length(line.numbers) > 0) {
    file.string <- paste(by.line[-line.numbers], collapse = '\n')
  }

  #------------------------------------------------------------------------
  # Splitting by element indicator
  by.element <- str_split(file.string, '[\n]*##[.$]*')[[1]][-1]

  # Splitting by equals sign
  tag.matrix <- str_split_fixed(by.element, '=', n = 2)
  tag.matrix[, 1] <- str_trim(tag.matrix[, 1], 'both')

  # Parsing comments within lines
  inline.comments <- str_extract_all(tag.matrix[, 2], '(?<=\\$\\$).*')
  inline.comments <- lapply(inline.comments, str_trim, 'both')
  inline.index <- sapply(inline.comments, function (x) { length(x) > 0 })

  new.comments <- as.list(inline.comments[inline.index])
  names(new.comments) <- tag.matrix[inline.index, 1]
  comments <- c(comments, new.comments)

  # Removing comment strings from entries
  tag.matrix[, 2] <- str_replace_all(tag.matrix[, 2], '\\$\\$.*', '')
  tag.matrix[, 2] <- str_trim(tag.matrix[, 2], 'both')

  #------------------------------------------------------------------------
  # Looping through entries to check where header information extends,
  # breaking iteration if there are repeats
  last.item <- 1
  required.tags <- c('TITLE', 'ORIGIN', 'OWNER')
  found.tags <- rep(FALSE, 3)
  names(found.tags) <- required.tags

  # Restricting check to the first 10 rows
  for (i in 1:10) {
    tag <- tag.matrix[i, 1]
    if (tag %in% required.tags) {

      if (found.tags[tag] == TRUE) break
      else found.tags[tag] <- TRUE

      last.item <- i
    }
  }

  header <- as.list(tag.matrix[1:last.item, 2])
  names(header) <- tag.matrix[1:last.item, 1]

  # Removing header from tag matrix
  tag.matrix <- tag.matrix[-(1:last.item), ]

  # If more TITLE tags follow, then the header is extended to include
  # everything until the next TITLE (with the assumption that the
  # subsequent titles refer to blocks).
  if ('TITLE' %in% tag.matrix[, 1]) {
    last.item <- which('TITLE' == tag.matrix[, 1])[1] - 1

    if (last.item != 0) {
      new.header <- as.list(tag.matrix[1:last.item, 2])
      names(new.header) <- tag.matrix[1:last.item, 1]
      header <- c(header, new.header)

      # Removing header from tag matrix
      tag.matrix <- tag.matrix[-(1:last.item), ]
    }
  }

  #------------------------------------------------------------------------
  # Looping through blocks
  blocks <- list()
  block.number <- 1
  current.index <- 1

  while (TRUE) {
    if (current.index >= nrow(tag.matrix)) break
    if (tag.matrix[current.index, 1] == 'END') break

    # Looking for next END tag
    end.indexes <- which(tag.matrix[-(1:current.index), 1] == 'END')

    if (length(end.indexes) == 0) end.index <- nrow(tag.matrix)
    else end.index <- end.indexes[1] + current.index - 1

    # Subsetting block tags
    s.tag.matrix <- tag.matrix[current.index:end.index, ]
    current.index <- end.index + 2

    # If block contains NTUPLES, they must be extracted and reformatted
    ntuples <- list()
    if (sum(str_detect(s.tag.matrix[, 1], 'NTUPLES')) > 0) {

      # First check the number of NTUPLES and END NTUPLES statements
      logic.start <- str_detect(s.tag.matrix[, 1], '^NTUPLES')
      names.start <- s.tag.matrix[logic.start, 2]
      len.start <- sum(logic.start)

      logic.end <- str_detect(s.tag.matrix[, 1], '^END NTUPLES')
      names.end <- s.tag.matrix[logic.end, 2]
      len.end <- sum(logic.end)

      if (len.start != len.end) {
        msg <- 'NTUPLES sections are unbounded'
        stop(msg)
      }

      if (!all(names.start == names.end)) {
        msg <- 'Mismatch in NTUPLES section names'
        stop(msg)
      }

      # If everything is fine, loop through ntuples tags
      for (i in 1:len.start) {

        name <- names.start[[i]]
        start <- which(logic.start)[i]
        end <- which(logic.end)[i]

        temp.matrix <- s.tag.matrix[(start + 1):(end - 1), ]
        s.tag.matrix <- s.tag.matrix[-(start:end), ]

        # Finding pages
        page.index <- which(str_detect(temp.matrix[, 1], 'PAGE'))
        page.index <- page.index

        # Last page is assumed to go until end of tuples section
        page.index <- c(page.index, nrow(temp.matrix) + 1)

        # Adding everything before pages
        ntuples[[name]] <- as.list(temp.matrix[1:(page.index[1] - 1), 2])
        names(ntuples[[name]]) <- temp.matrix[1:(page.index[1] - 1), 1]

        n.pages <- length(page.index) - 1

        if (n.pages > 0) {

          ntuples[[name]][['PAGES']] <- list()
          current.page <- 1

          # Looping through pages
          for (j in 1:n.pages) {
            start <- page.index[j] + 1
            end <- page.index[j + 1] - 1
            page.content <- temp.matrix[start:end, 2]
            page.names <- temp.matrix[start:end, 1]
            ntuples[[name]][['PAGES']][[current.page]] <- as.list(page.content)
            names(ntuples[[name]][['PAGES']][[current.page]]) <- page.names
            current.page <- current.page + 1
          }
        }
      }
    }

    # Storing raw string data
    current.block <- as.list(s.tag.matrix[, 2])
    names(current.block) <- s.tag.matrix[, 1]

    # If NTUPLES detected, tacking them on
    if (length(ntuples) > 0) current.block[['NTUPLES']] <- ntuples

    blocks[[block.number]] <- current.block
    block.number <- block.number + 1
  }

  if (current.index < (nrow(tag.matrix) - 1)) {
    msg <- 'Block parsing stopped before end of file'
    warning(msg)
  }

  #------------------------------------------------------------------------
  # Combining data
  out <- list(header = header, blocks = blocks, comments = comments)

  # Feeding through the process function
  process_jcamp(out, tags = process.tags, entries = process.entries)
}

#------------------------------------------------------------------------
#' Process imported JCAMP data
#'
#' Converts string key-value pairs from parsed JCAMP file into numeric and
#' matrix formats where applicable.
#'
#' @param jcamp_list The output of read_jcamp().
#' @param tags TRUE to convert entry tags. All tags are converted to
#'             lowercase with dots replacing spaces. Some tags known to
#'             take multiple forms ('JCAMPDX' vs. 'JCAMP-DX') are converted
#'             to a single fixed form (full list to come).
#' @param entries TRUE to convert entries using a system of ad-hoc rules.
#'                All entries designated with '<>' are converted to strings.
#'                All other entries are converted to numeric if possible.
#'                Lists and tables are converted to matrices if possible.
#' @param ... Extra arguments passed into process_jcamp_tag().
#
#' @return ...
#'
process_jcamp <- function(jcamp_list, tags = TRUE, entries = TRUE, ...) {

  # Basic input check
  if (!all(c('blocks', 'header') %in% names(jcamp_list))) {
    msg <- 'Input must be a JCAMP list similar to the output of read_jcamp()'
    stop(msg)
  }

  # Copying over output
  out <- jcamp_list

  #------------------------------------------------------------------------
  # Entries first
  if (entries) {

    # Header data
    out[['header']] <- lapply(out[['header']], process_jcamp_entry)

    # Looping through blocks
    for (i in 1:length(out[['blocks']])) {

      # Processing all items other than ntuples directly
      logic <-  names(out[['blocks']][[i]]) != 'NTUPLES'
      out[['blocks']][[i]][logic] <- lapply(out[['blocks']][[i]][logic],
                                            process_jcamp_entry)

      # If there are no ntuples then skip the rest of the iteration
      if (!'NTUPLES' %in% names(out[['blocks']][[i]])) next

      # Processing all ntuples descriptors with a separator and combining
      # them into a data.frame
      for (j in 1:length(out[['blocks']][[i]][['NTUPLES']]) ) {
        out[['blocks']][[i]][['NTUPLES']][[j]] <- with(out[['blocks']][[i]], {
          logic <- names(NTUPLES[[j]]) != 'PAGES'
          NTUPLES[[j]][logic] <- lapply(NTUPLES[[j]][logic],
                                        process_jcamp_entry,
                                        sep = '[ \t]*,[ \t]*')

          # Checking on descriptor length (they should all be the same)
          all.lengths <- sapply(NTUPLES[[j]][logic], length)

          if (length(unique(all.lengths)) > 1) {
            msg <- paste('Not all NTUPLE descriptors are the same length,',
                         'there may have been a parsing error')
            warning(msg)
          }

          longest <- all.lengths[all.lengths == max(all.lengths)]

          if (longest[1] < 3) {
            msg <- paste('NTUPLE content does not seem to contain more than',
                         'one data set')
            warning(msg)
          }

          # Combining descriptors into data frame
          descriptors <- NTUPLES[[j]][logic][names(longest)]
          d <- do.call(cbind,
                       lapply(descriptors, function (x) as.data.frame(x)))
          colnames(d) <- names(descriptors)

          if ('SYMBOL' %in% names(descriptors)) {
            rownames(d) <- descriptors[['SYMBOL']]
          }

          # Replacing individual descriptor tags with a single descriptor
          # data frame
          logic <- names(NTUPLES[[j]]) %in% names(descriptors)
          NTUPLES[[j]][logic] <- NULL
          NTUPLES[[j]][['DESCRIPTORS']] <- d

          # Looping through the pages
          NTUPLES[[j]][['PAGES']] <- lapply(NTUPLES[[j]][['PAGES']],
                                            process_jcamp_entry)

          # Outputting formatted NTUPLES[[j]]
          NTUPLES[[j]]
        })
      }
    }

  }

  #------------------------------------------------------------------------
  # Validation

  # May introduce later

  #------------------------------------------------------------------------
  # Tags

  if (tags) {

    # Header data
    names(out[['header']]) <- sapply(names(out[['header']]),
                                     process_jcamp_tag, ...)

    # Looping through blocks
    for (i in 1:length(out[['blocks']])) {

      # Processing all non ntuple items
      logic <- names(out[['blocks']][[i]]) != 'NTUPLES'

      names(out[['blocks']][[i]])[logic] <-
        sapply(names(out[['blocks']][[i]])[logic], process_jcamp_tag, ...)

      # If there are no ntuples then skip the rest of the iteration
      if (!'NTUPLES' %in% names(out[['blocks']][[i]])) next

      # Processing ntuple internals first
      for (j in 1:length(out[['blocks']][[i]][['NTUPLES']]) ) {

        out[['blocks']][[i]][['NTUPLES']][[j]] <- with(out[['blocks']][[i]], {

          # Modifying dataframe names
          NTUPLES[[j]][['DESCRIPTORS']] <- with(NTUPLES[[j]], {
            colnames(DESCRIPTORS) <- sapply(colnames(DESCRIPTORS),
                                            process_jcamp_tag, ...)
            rownames(DESCRIPTORS) <- sapply(rownames(DESCRIPTORS),
                                            process_jcamp_tag, ...)
            DESCRIPTORS
          })

          # Looping through the PAGES
          for (k in 1:length(NTUPLES[[j]][['PAGES']])) {
            NTUPLES[[j]][['PAGES']][[k]] <- with(NTUPLES[[j]], {
              colnames(PAGES[[k]]) <- sapply(colnames(PAGES[[k]]),
                                             process_jcamp_tag, ...)
              PAGES[[k]]
            })
          }

          names(NTUPLES[[j]]) <- sapply(names(NTUPLES[[j]]),
                                        process_jcamp_tag, ...)
          NTUPLES[[j]]
        })
      }

      # Processing the ntuple list itself
      out[['blocks']][[i]][['NTUPLES']] <- with(out[['blocks']][[i]], {
        names(NTUPLES) <- sapply(names(NTUPLES), process_jcamp_tag, ...)
        NTUPLES
      })

      logic <- names(out[['blocks']][[i]]) == 'NTUPLES'
      names(out[['blocks']][[i]])[logic] <- process_jcamp_tag('NTUPLES')
    }
  }

  return(out)
}

#------------------------------------------------------------------------
#' Process JCAMP entry
#'
#' Converts single JCAMP character string entry into numeric vector or
#' data.frame types where possible. If no conversion options found, returns
#' character string unchanged.
#'
#' @param jcamp.entry Single character string from the output of read_jcamp().
#' @param sep A regex string used to split entries into vectors. By default,
#'            no splitting is performed.
#'
#' @return Character string, numeric vector, or data.frame depending on
#'         the form of the entry
#'
process_jcamp_entry <- function(jcamp.entry, sep = NULL) {

  by.line <- str_split(jcamp.entry, '[\r\n]+')[[1]]

  # Processing for single line items
  if (length(by.line) == 1) {

    # If there is a descriptor tag of the form (0..10), the entry is
    # a vector
    d.vector <- '^[ \t]*\\(\\d+\\.\\.\\d+\\)'
    if (str_detect(by.line[1], d.vector)) {
      descriptor <- str_extract(by.line[1], d.vector)
      extent <- str_extract(descriptor, c('(?<=\\()\\d+', '\\d+(?=\\))'))
      extent <- as.numeric(extent)

      content <- str_replace(by.line[1], d.vector, '')
      formatted <- as.numeric(str_split(content, '[ ,\n]+')[[1]])

      if (length(formatted) != (extent[2] - extent[1] + 1)) {
        msg <- sprintf('Unexpected vector length processing: \n%s', jcamp.entry)
        stop(msg)
      }
    }
    # If surrounded by '<>', character, character
    else if (str_detect(jcamp.entry, '<.*>')) {
      formatted <- as.character(str_replace(jcamp.entry, '<(.*)>', '\\1'))
      return(formatted)
    }
    # Otherwise, attempt spitting is a separator is provided
    else if (!is.null(sep)) {
      formatted <- str_split(jcamp.entry, sep)[[1]]
    }
    else {
      formatted <- jcamp.entry
    }

    # Whether splitting required or not, try conversion to numeric
    formatted <- tryCatch({as.numeric(formatted)},
                          warning = function(w) {formatted},
                          error = function(e) {formatted})
  }
  # Multi-line items
  else {
    # Parse based on the nature of the format descriptor tag,
    # e.g. (0..10) vector, (XY..XY) pairs, or (X++(Y..Y)) spectrum
    d.vector <- '\\(\\d+\\.\\.\\d+\\)'
    d.pairs <- '\\(\\w{2}\\.\\.\\w{2}\\)'
    d.spectrum <- '\\(X\\+\\+\\(\\w\\.\\.\\w\\)\\)'

    if (str_detect(by.line[1], d.vector)) {
      descriptor <- str_extract(by.line[1], d.vector)
      extent <- str_extract(descriptor, c('(?<=\\()\\d+', '\\d+(?=\\))'))
      extent <- as.numeric(extent)

      content <- paste(by.line[2:length(by.line)], collapse = '\n')
      formatted <- str_split(content, '[ \n]+')[[1]]

      # Checking to see if it's a vector of strings
      n.strings <- sum(grepl('<.*>', formatted))

      if ((n.strings > 0) && (n.strings < length(formatted))) {
        msg <- sprintf('Error processing mixed type vector: \n%s', jcamp.entry)
        stop(msg)
      }

      if (n.strings > 0) {
        formatted <- as.character(gsub('<(.*)>', '\\1', formatted))
      } else {
        formatted <- as.numeric(formatted)
      }

      if (length(formatted) != (extent[2] - extent[1] + 1)) {
        msg <- sprintf('Unexpected vector length processing: \n%s', jcamp.entry)
        stop(msg)
      }
    }
    else if (str_detect(by.line[1], d.pairs)) {
      # Picking off xy names
      descriptor <- str_extract(by.line[1], d.pairs)
      xy.names <- str_extract(descriptor, c('(?<=\\()\\w', '(?<=\\(\\w)\\w'))

      content <- paste(by.line[2:length(by.line)], collapse = '\n')

      # Splitting multiple pairs into single pairs
      xy.pairs <- str_trim(str_split(content, '[;\r\n]+')[[1]], 'both')

      # Spliting single pairs into values
      values <- str_split_fixed(xy.pairs, '[ ,]+', n = 2)

      # Converting to numeric and labelling columns
      mode(values) <- 'numeric'
      colnames(values) <- xy.names

      formatted <- values
    }
    else if (str_detect(by.line[1], d.spectrum)) {
      # Picking off xy names
      descriptor <- str_extract(by.line[1], d.spectrum)
      x.name <- str_extract(descriptor, '(?<=\\()\\w(?=\\+\\+)')
      y.name <- str_extract(descriptor, '(?<=\\(\\w\\+\\+\\()\\w')

      # Dropping format string
      by.line <- by.line[-1]

      # First, convert all spaces to a single character width
      by.line <- str_replace_all(by.line, '\\s+', ' ')

      # Then substitute a space before all leading characters
      by.line <- str_replace_all(by.line, '(?<![ ])([a-zA-Z@%+-])', ' \\1')

      # Check compression mode
      if ( any(str_detect(by.line, '[sS-Z]')) ) {
        mode <- 3 # DUP
      } else if ( any(str_detect(by.line, '[j-rJ-R%]')) ) {
        mode <- 2 # DIF
      } else if ( any(str_detect(by.line, '[a-iA-I@]')) ) {
        mode <- 1 # SQZ
      } else {
        mode <- 0
      }

      # If the data has encoded repeats, expand them
      if ( mode == 3 ) {

        # Pick off repeated elements
        pattern <- '[0-9a-zA-Z@%+-]+ [sS-Z][0-9]*'
        repeats <- unique(unlist(str_match_all(by.line, pattern)))

        # Translate repeats
        patterns <- c('S'='1', 'T'='2', 'U'='3', 'V'='4', 'W'='5',
                      'X'='6', 'Y'='7', 'Z'='8', 's'='9')
        replacements <- str_replace_all(repeats, patterns)

        # Split repeat pairs and apply the repeat
        replacements <- str_split(replacements, ' ')
        f <- function(x) paste(rep(x[1], as.numeric(x[2])), collapse = ' ')
        replacements <- unlist(lapply(replacements, f))

        # Finally, re-insert the repeated values
        names(replacements) <- repeats
        by.line <- str_replace_all(by.line, replacements)
      }

      # Convert squeezed characters
      if ( mode >= 1) {

        # Translate characters
        patterns <- c('A'='1', 'B'='2', 'C'='3', 'D'='4', 'E'='5',
                      'F'='6', 'G'='7', 'H'='8', 'I'='9',
                      'a'='-1', 'b'='-2', 'c'='-3', 'd'='-4', 'e'='-5',
                      'f'='-6', 'g'='-7', 'h'='-8', 'i'='-9', '@'=0)

        by.line <- str_replace_all(by.line, patterns)
      }

      # At this point, split by spaces
      by.line <- lapply(by.line, str_trim)
      by.line <- str_split(by.line, ' ')

      # If the data is stored as differences, apply them now
      if ( mode >= 2 ) {

        # Translate characters
        patterns <- c('J'='1', 'K'='2', 'L'='3', 'M'='4', 'N'='5',
                      'O'='6', 'P'='7', 'Q'='8', 'R'='9',
                      'j'='-1', 'k'='-2', 'l'='-3', 'm'='-4',
                      'n'='-5', 'o'='-6', 'p'='-7', 'q'='-8', 'r'='-9', '%'=0)

        # Double check to make sure there are no non-differences present
        f <- function(x) all(str_detect(x[-(1:2)], '[j-rJ-R%]'))
        all.differences <- lapply(by.line, f)
        msg <- 'Absolute values found among DIF compressed data, aborting.'
        if ( any(! unlist(all.differences)) ) stop(msg)

        by.line <- lapply(by.line, str_replace_all, patterns)

        # Convert to numerical format and apply cumulative sum
        f <- function(x) c(as.numeric(x[1]), cumsum(as.numeric(x)[-1]))
        by.line <- lapply(by.line, f)

        # With the difference comparison mode, the last value of each line
        # should match the second value of following line
        last.values <- unlist(lapply(by.line, function(x) `[`(x, length(x))))
        second.values <- unlist(lapply(by.line, function(x) `[`(x, 2)))

        n <- length(by.line)
        msg <- 'y value data check failed in decompressing DIF data, aborting.'
        differences <- abs(last.values[1:(n-1)] - second.values[2:n])
        if ( any(differences > 1e-6) ) stop(msg)

        # Stripping these last values
        by.line[-n] <- lapply(by.line[-n], function(x) x[-length(x)])

      } else {
        # Otherwise, just convert to numeric
        by.line <- lapply(by.line, as.numeric)
      }

      # Picking off the x values
      x.values <- unlist(lapply(by.line, `[`, 1))

      # Counting observations per line
      counts <- unlist(lapply(by.line, function(x) length(x) - 1))

      # Ensuring that delta x remains consistent throughout
      n <- length(by.line)
      delta.x <- (x.values[2:n] - x.values[1:(n-1)])/counts[1:(n-1)]

      msg <- 'x value data check failed in data import, aborting.'
      if ( diff(range(delta.x))/mean(delta.x) > 1e-3 ) stop(msg)

      delta.x <- median(delta.x)

      n.out <- sum(counts)
      x.out <- seq(x.values[1], by = delta.x, length.out = n.out)
      y.out <- unlist(lapply(by.line, `[`, -1))

      formatted <- data.frame(x.out, y.out)
      colnames(formatted) <- c(x.name, y.name)
    }
    else {
      # If no descriptive marker is found, then it might be a
      # multi-line string that should be stitched back together.
      formatted <- paste(by.line, collapse = '\n')
    }
  }

  return(formatted)
}

#------------------------------------------------------------------------
#' Process JCAMP tag
#'
#' Processes a JCAMP tag into standard form. This consists of selecting
#' a case conversion function, such as tolower(), as well as a map of
#' specfic conversions to make. Note, case conversions are applied after
#' the mapped conversions.
#'
#' @param jcamp.tag Single character string corresponding to list names
#'                  from the output of read_jcamp().
#' @param f_case The function used to perform the case conversion, e.g.
#'               tolower().
#' @param tag.space A string used to replace empty spaces and underscores.
#' @param tag.map A vector of strings corresponding to new tags,
#'                where the vector names correspond to the tags that will
#'                be renamed. NULL avoids map conversions while NA uses
#'                the default conversions below:
#'
#' @return A renamed character string.
#'
#'
process_jcamp_tag <- function(jcamp.tag, f_case = tolower,
                              tag.space = '.', tag.map = NA) {

  # Converting case
  jcamp.tag <- f_case(jcamp.tag)

  # Replacing spaces
  jcamp.tag <- str_replace_all(jcamp.tag, '[ \t_]+', tag.space)

  # Specifying default tag.map
  default.map <- c('rev'='reverse')

  # If no tag.map specified, setting default case
  if ((length(tag.map) == 1) && is.na(tag.map)) {
    tag.map <- default.map
  }

  jcamp.tag <- ifelse(jcamp.tag %in% names(tag.map),
                      tag.map[jcamp.tag], jcamp.tag)
}

#------------------------------------------------------------------------
#' Flatten JCAMP list file
#'
#' As most JCAMP files contain only one block/ntuple, it's not necessary
#' to use deep nesting for all the entries. This function checks the
#' the number of blocks/ntuples and then stores header, block, and ntuple
#' entries in a single flat list. Only the ntuple element remains a
#' list, containing the elements 'descriptors' and 'pages'. Note, tag names
#' are switched between processed and non-processed forms automatically.
#' If there are name conlicts between header and block elements, only
#' the block elements are kept.
#'
#' @param jcamp_list The list output of read_jcamp() or process_jcamp().
#'
#' @return A flattened list containing all jcamp entries.
#'
#'
flatten_jcamp <- function(jcamp_list) {

  # Initializing output
  out <- jcamp_list

  # Moving headers if they exist
  if ('header' %in% names(out)) {
    out[names(out$header)] <- out$header
    out$header <- NULL
  }

  # Moving blocks if there aren't too many
  if ('blocks' %in% names(out)) {
    if (length(out$blocks) > 1) {
      msg <- "JCAMP lists with more than one block can't be flattened"
      stop(msg)
    }

    out[names(out$blocks[[1]])] <- out$blocks[[1]]
    out$blocks <- NULL
  }

  # Choosing appropraite ntuples tag
  if ('NTUPLES' %in% names(out)) {
    ntuples.string <- 'NTUPLES'
    ntuples <- out[['NTUPLES']]
  } else if (process_jcamp_tag('NTUPLES') %in% names(out)) {
    ntuples.string <- process_jcamp_tag('NTUPLES')
    ntuples <- out[[ntuples.string]]
  } else {
    ntuples.string = NA
  }

  # Moving ntuples if there aren't too many
  if (!is.na(ntuples.string)) {
    if (length(ntuples) > 1) {
      msg <- "JCAMP lists with more than one ntuple entry can't be flattened"
      stop(msg)
    }

    out[[ntuples.string]] <- ntuples[[1]]
  }

  # Returning
  out
}

