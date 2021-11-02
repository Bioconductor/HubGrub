#' Discover what data source types are available in the specific Hub
#'
#' @param hub A hub object, either AnnotationHub or ExperimentHub. 
#'
#' @importFrom AnnotationHub mcols
#' @importFrom dplyr arrange desc %>%
#'
#' @return A tibble sorted in descending order.
#' 
#' @examples
#' ah = AnnotationHub()
#' discoverData(ah)

discoverData <- function(hub) {
    stopifnot(is(hub, "AnnotationHub") || is(hub, "ExperimentHub"))

    table(mcols(hub)$rdataclass) %>% sort(decreasing = TRUE)
}
