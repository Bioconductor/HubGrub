test_that("'importData()' works", {
    library(AnnotationHub)
    library(GenomicRanges)
    ah = AnnotationHub()
    which <- GRanges(c("chr2", "chr2"), IRanges(c(1, 300), c(400, 1000000)))

    dat1 <- importData(ah, "AH49544", which)
    expect_gte(length(dat1), 30384L)

    expect_error(importData(ah))
    expect_error(importData(ah, AH49544))
    expect_error(importData(ah, "AH49544"))
})
