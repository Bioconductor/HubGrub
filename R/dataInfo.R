#' Display select information about hub object
#'
#' @param hub A hub object, either AnnotationHub or ExperimentHub.
#' @param fileType The type of file the user in interested in exploring, e.g. 
#'     'BigWigFile'.
#'
#' @importFrom AnnotationHub subset
#' @importFrom tibble as_tibble
#' @importFrom dplyr select
#' 
#' @return A DataFrame.
#'
#' @examples
#' ah = AnnotationHub()
#' dataInfo(ah, "BigWigFile")
dataInfo <- function(hub, fileType) {
    stopifnot(
        is(hub, "AnnotationHub") || is(hub, "ExperimentHub"),
        .is_scalar_character(fileType),
        fileType %in% unique(hub$rdataclass)
    )

    sub_hub <- AnnotationHub::subset(hub, rdataclass == fileType)

    datInfo <- as_tibble(mcols(sub_hub), rownames = "ID") |>
        dplyr::select(
            -c(
                "coordinate_1_based",
                "rdatadateadded",
                "preparerclass",
                "rdataclass",
                "rdatapath",
                "sourceurl",
                "sourcetype"
            )
       )
    datInfo
}
