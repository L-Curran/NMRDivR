# Definition of a super class structure for NMR data.

#------------------------------------------------------------------------
#' Super class combining NMR data and scan parameters.
#'
#' See NMRData1D or NMRData2D.
#'
#' @slot processed A data.frame containing chemical shift and intensity data.
#' @slot parameters A list of generic parameters.
#' @slot acqus A list of parameters specific to acqus file.
#' @slot procs A list of parameters sepcific to procs file.
#'
#' @name NMRData-class
#' @export
NMRData <- setClass("NMRData",
                    slots = c(processed = "data.frame",
                              parameters = "list",
                              acqus = "list",
                              procs = "list"))



#========================================================================>
# Basic setter and getter functions
#========================================================================>


#------------------------------------------------------------------------
#' @templateVar slot processed 
#' @template NMRData_access
#' @name processed
#' @export
setGeneric("processed", 
           function(object, ...) standardGeneric("processed"))

#' @rdname processed
#' @export
setMethod("processed", "NMRData", 
          function(object) object@processed)

#' @templateVar slot processed 
#' @template NMRData_replacement
#' @name processed-set
#' @export
setGeneric("processed<-", 
           function(object, value) standardGeneric("processed<-"))

#' @rdname processed-set
#' @export
setReplaceMethod("processed", "NMRData",
                 function(object, value) {
                   object@processed <- value
                   validObject(object)
                   object 
                 })

#------------------------------------------------------------------------
#' @templateVar slot parameters 
#' @template NMRData_access
#' @name parameters
#' @export
setGeneric("parameters", 
           function(object, ...) standardGeneric("parameters"))

#' @rdname parameters
#' @export
setMethod("parameters", "NMRData", 
          function(object) object@parameters)

#' @templateVar slot parameters 
#' @template NMRData_replacement
#' @name parameters-set
#' @export
setGeneric("parameters<-", 
           function(object, value) standardGeneric("parameters<-"))

#' @rdname parameters-set
#' @export
setReplaceMethod("parameters", "NMRData",
                 function(object, value) {
                   object@parameters <- value
                   validObject(object)
                   object 
                 })

#------------------------------------------------------------------------
#' @templateVar slot procs 
#' @template NMRData_access
#' @name procs
#' @export
setGeneric("procs", 
           function(object, ...) standardGeneric("procs"))

#' @rdname procs
#' @export
setMethod("procs", "NMRData", 
          function(object) object@procs)

#' @templateVar slot procs 
#' @template NMRData_replacement
#' @name procs-set
#' @export
setGeneric("procs<-", 
           function(object, value) standardGeneric("procs<-"))

#' @rdname procs-set
#' @export
setReplaceMethod("procs", "NMRData",
                 function(object, value) {
                   object@procs <- value
                   validObject(object)
                   object 
                 })

#------------------------------------------------------------------------
#' @templateVar slot acqus 
#' @template NMRData_access
#' @name acqus
#' @export
setGeneric("acqus", 
           function(object, ...) standardGeneric("acqus"))

#' @rdname acqus
#' @export
setMethod("acqus", "NMRData", 
          function(object) object@acqus)

#' @templateVar slot acqus 
#' @template NMRData_replacement
#' @name acqus-set
#' @export
setGeneric("acqus<-", 
           function(object, value) standardGeneric("acqus<-"))

#' @rdname acqus-set
#' @export
setReplaceMethod("acqus", "NMRData",
                 function(object, value) {
                   object@acqus <- value
                   validObject(object)
                   object 
                 })



#========================================================================>
# Defining list and data.frame like behaviour
#========================================================================>

#------------------------------------------------------------------------
#' Convert NMRData object to list
#'
#' Outputs all slot values as entries in a list.
#'
#' @param x NMRData object.
#'
#' @export
as.list.NMRData <- function(x) {
  list(processed = x@processed, parameters = x@parameters, 
       procs = x@procs, acqus = x@acqus)
}

setMethod("as.list", "NMRData", as.list.NMRData)

#------------------------------------------------------------------------
#' Convert NMRData object to data.frame
#'
#' Outputs the processed data.frame, dropping acqus and procs lists.
#'
#' @param x NMRData object.
#'
#' @export
as.data.frame.NMRData <- function(x) {
  x@processed
}

setMethod("as.data.frame", "NMRData", as.data.frame.NMRData)


