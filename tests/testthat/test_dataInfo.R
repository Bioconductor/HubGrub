test_that("'dataInfo()' works", {
    library(AnnotationHub)
    ah = AnnotationHub()
    

    dat1 <- dataInfo(ah, "BigWigFile")
    expect_identical(class(dat1), c("tbl_df", "tbl", "data.frame"))
    expect_gte(dim(dat1)[1], 10247L)
    expect_identical(dim(dat1)[2], 9L)
    expect_identical(colnames(dat1), c("ID", "title", "dataprovider", "species", 
        "taxonomyid", "genome", "description", "maintainer", "tags"))

    dat2 <- dataInfo(ah, "VcfFile")
    expect_gte(dim(dat2)[1], 8L)
    expect_identical(dim(dat2)[2], 9L)
    expect_identical(colnames(dat2), c("ID", "title", "dataprovider", "species", 
        "taxonomyid", "genome", "description", "maintainer", "tags"))

    expect_error(dataInfo(ah))
    expect_error(dataInfo(ah, c("BigWigFile", "VcfFile")))
})
