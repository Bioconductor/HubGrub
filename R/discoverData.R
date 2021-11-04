#' Discover what data source types are available in the specific Hub
#'
#' @param hub A hub object, either AnnotationHub or ExperimentHub. 
#'
#' @importFrom AnnotationHub mcols
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter count arrange desc
#'
#' @return A tibble sorted in descending order.
#' 
#' @examples
#' ah = AnnotationHub()
#' discoverData(ah)
discoverData <- function(hub, fileType = types()) {
    stopifnot(
        is(hub, "AnnotationHub") || is(hub, "ExperimentHub"),
        all(fileType %in% mcols(hub)$rdataclass)
    )

    res <- as_tibble(DataClass = mcols(hub)$rdataclass) |>
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
