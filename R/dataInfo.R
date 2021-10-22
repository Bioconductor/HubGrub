#' Display select information about hub object
#'
#' @param hub A hub object, either AnnotationHub or ExperimentHub.
#' @param sourceType For the hub resources of interest, the format of the 
#'     original resource, e.g., BED file.
#'
#' @importFrom AnnotationHub subset
#' 
#' @return A DataFrame.
#'
#' @examples
#' ah = AnnotationHub()
#' dataInfo(ah, "BigWig")
dataInfo <- function(hub, sourceType) {
    stopifnot(
        is(hub, "AnnotationHub") || is(hub, "ExperimentHub"),
        is.character(sourceType),
        sourceType %in% unique(hub$sourcetype)
    )

    sub_hub <- subset(hub, sourcetype == sourceType)

    mcols(sub_hub)[,c("dataprovider", "species", "sourcetype")]
}
