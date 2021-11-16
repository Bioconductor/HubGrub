test_that("'discoverData()' works", {
    library(AnnotationHub)
    ah = AnnotationHub()

    dat1 <- discoverData(ah)
    expect_identical(class(dat1), c("tbl_df", "tbl", "data.frame"))
    expect_identical(dim(dat1), c(2L, 2L))
    expect_identical(colnames(dat1), c("DataClass", "n"))

    dat2 <- discoverData(ah, "BigWigFile")
    expect_identical(class(dat2), c("tbl_df", "tbl", "data.frame"))
    expect_identical(dim(dat2), c(1L, 2L))
    expect_identical(colnames(dat2), c("DataClass", "n"))

    expect_error(discoverData(ah, "test"))
})
