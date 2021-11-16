#' Discover what data source types are available in the specific Hub
#'
#' @param hub A hub object, either AnnotationHub or ExperimentHub. 
#' @param fileType The type of file the user is intersted in exploring, e.g.
#'     'BigWigFile'. 
#'
#' @importFrom AnnotationHub mcols
#' @importFrom tibble tibble
#' @importFrom dplyr filter count arrange desc
#'
#' @return A tibble sorted in descending order.
#' @export
#'  
#' @examples
#' library(AnnotationHub)
#' ah = AnnotationHub()
#' discoverData(ah)
discoverData <- function(hub, fileType = types()) {
    stopifnot(
        is(hub, "AnnotationHub") || is(hub, "ExperimentHub"),
        all(fileType %in% mcols(hub)$rdataclass)
    )

    DataClass <- n <- NULL

    res <- tibble(DataClass = mcols(hub)$rdataclass) |>
        filter(DataClass %in% fileType) |>
        count(DataClass) |>
        arrange(desc(n))
    res
}

types <- function() {
    c(
        "BigWigFile",
        "VcfFile"
    )
}
