#' Display select information about hub object
#'
#' @param hub A hub object, either AnnotationHub or ExperimentHub.
#' @param rdataclass For the hub resources of interest, the class of the R 
#'     object used to represent the object when imported into R, e.g., 'GRanges'.
#'
#' @importFrom AnnotationHub subset
#' 
#' @return A DataFrame.
#'
#' @examples
#' ah = AnnotationHub()
#' dataInfo(ah, "BigWigFile")
dataInfo <- function(hub, rdataclass) {
    stopifnot(
        is(hub, "AnnotationHub") || is(hub, "ExperimentHub"),
        is.character(rdataclass),
        rdataclass %in% unique(hub$rdataclass)
    )

    sub_hub <- subset(hub, rdataclass == rdataclass)

    mcols(sub_hub) ## can display only the columns we think are of interest...
}
