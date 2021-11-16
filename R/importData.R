#' @importFrom AnnotationDbi dbFileConnect dbfile
#' @importFrom DBI dbGetQuery dbDisconnect
#' @noRd
.resourceLocation <- function (x) {
    query <-
        'SELECT DISTINCT resources.ah_id, location_prefixes.location_prefix
        FROM resources, location_prefixes
        WHERE location_prefixes.id == resources.location_prefix_id'

    con <- dbFileConnect(dbfile(x))
    mat <- dbGetQuery(con, query)
    dbDisconnect(con)
    mat
}

.idsInfo <- function (hub, ids) {
    dat <- .resourceLocation(hub)
    sub_dat <- match(ids, dat$ah_id)

    res <- dat[sub_dat,]
    res
}

#' Import data of hub object
#' 
#' @param hub A hub object, either AnnotationHub or ExperimentHub
#' @param ids The hub id(s) of interest.
#' @param which A range data structure, like a 'GRanges'.
#' 
#' @importFrom rtracklayer import
#' 
#' @return A 'GRanges'.
#' @export
#'  
#' @examples
#' library(AnnotationHub)
#' library(GenomicRanges)
#' ah = AnnotationHub()
#' which <- GRanges(c("chr2", "chr2"), IRanges(c(1, 300), c(400, 1000000)))
#' importData(ah, "AH49544", which)
importData <- function (hub, ids, which) {
    tbl <- .idsInfo(hub, ids)

    datPath <- mcols(hub[ids])$rdatapath

    path <- paste(tbl$location_prefix, datPath, sep = "")    

    import(path, which = which)
}
