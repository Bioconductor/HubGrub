# HubGrub package

### Used to discover information about resources in the hubs

Goals:

* Create a table of file type information
	* What types of file are in the hubs?
	* How many of each type?
* Figure out how to get information for various file types
* Generate a coherent set of files based on certain characteristics
	* species
	* source
	* project

## 10/22/21 work

I've created an initial function `discoverData()` that generates a table for
source types. For example, in AnnotationHub there is:

```{r}
library(AnnotationHub)
ah = AnnotationHub()
> discoverData(ah)
# A tibble: 31 × 2
   name         value
   <chr>        <table>
 1 FASTA        14061
 2 BED          13018
 3 BigWig       12384
 4 GTF          10354
 5 ensembl       2719
 6 UCSC track    2217
 7 NCBI/UniProt  1725
 8 Chain         1115
 9 GFF            933
10 CSV            678
# … with 21 more rows
>
```

What I've discovered is there is also rdataclass for files. In some cases these
match the sourcetype, but in others they are different. The rdataclass is what
determines how the resource is loaded in R. So this might be of more interest to
us. This information look like this:

```{r}
> enframe(table(mcols(ah)$rdataclass)) %>% arrange(desc(value))
# A tibble: 29 × 2
   name          value
   <chr>         <table>
 1 GRanges       26599
 2 TwoBitFile    14056
 3 BigWigFile    10247
 4 EnsDb          2719
 5 Rle            2185
 6 OrgDb          1744
 7 ChainFile      1115
 8 TxDb            480
 9 SQLiteFile      421
10 Inparanoid8Db   268
# … with 19 more rows
>
```

For the 14061 FASTA files the rdataclass is:

```{r}
> table(mcols(fasta)$rdataclass)

AAStringSet      FaFile  TwoBitFile        TxDb
          1           3       13912         145
>
```

For the 13018 BED files the rdataclass is:

```{r}
> table(mcols(bed)$rdataclass)

   GRanges SQLiteFile     String
     13002          4         12
>
```

For the 12384 BigWig files the rdataclass is:

```{r}
> table(mcols(bigwig)$rdataclass)

BigWigFile        Rle     String
     10247       2136          1
>
```

For the 10354 GTF files the rdataclass is:

```{r}
> table(mcols(gtf)$rdataclass)

GRanges    list  String    TxDb
  10313      20       1      20
>
```

For the 2719 ensembl files the rdataclass is:

```{r}
> table(mcols(ensembl)$rdataclass)

EnsDb
 2719
>
```

For the 2217 UCSC track files the rdataclass is:

```{r}
> table(mcols(ucsc)$rdataclass)

GRanges
   2217
>
```

For the 1725 NCBI/UniProt files the rdataclass is:

```{r}
> table(mcols(ncbiUni)$rdataclass)

OrgDb
 1725
>
```

For the 1115 Chain files the rdataclass is:

```{r}
> table(mcols(chain)$rdataclass)

ChainFile
     1115
>
```

For the 933 GFF files the rdataclass is:

```{r}
> table(mcols(gff)$rdataclass)

GRanges    list    TxDb
    615       3     315
>
```

For the 678 CSV files the rdataclass is:

```{r}
> table(mcols(csv)$rdataclass)

   GRanges SQLiteFile     Tibble
       406        252         20
>
```

For the 268 Inparanoid files the rdataclass is:

```{r}
> sort(table(mcols(inparanoid)$rdataclass), decreasing = TRUE)
Inparanoid8Db
          268
>
```

For the 180 TSV files the rdataclass is:

```{r}
> sort(table(mcols(tsv)$rdataclass), decreasing = TRUE)

SQLiteFile       list  character data.frame     Tibble     SQLite
       148         17          6          4          3          2
>
```

For the 144 TwoBit files the rdataclass is:

```{r}
> sort(table(mcols(twobit)$rdataclass), decreasing = TRUE)
TwoBitFile
       144
>
```

For the 72 RData files the rdataclass is:

```{r}
> sort(table(mcols(rdata)$rdataclass), decreasing = TRUE)

data.frame       list    GRanges     igraph
        36         30          4          2
>
```

For the 70 XML files the rdataclass is:

```{r}
> sort(table(mcols(xml)$rdataclass), decreasing = TRUE)

    Tibble SQLiteFile data.table
        39         17         14
>
```

For the 57 VCF files the rdataclass is:

```{r}
> sort(table(mcols(vcf)$rdataclass), decreasing = TRUE)

    Rle VcfFile
     49       8
>
```

For the 42 RDS files the rdataclass is:

```{r}
> sort(table(mcols(rds)$rdataclass), decreasing = TRUE)
GRanges
     42
>
```

For the 21 TXT files the rdataclass is:

```{r}
sort(table(mcols(txt)$rdataclass), decreasing = TRUE)

      Rda character    sqlite
       19         1         1
>
```

For the 19 NCBI/ensembl files the rdataclass is:

```{r}
> sort(table(mcols(ncbiEnsembl)$rdataclass), decreasing = TRUE)
OrgDb
   19
```

For the 19 Zip files the rdataclass is:

```{r}
> sort(table(mcols(zip)$rdataclass), decreasing = TRUE)

                       data.frame data.frame, DNAStringSet, GRanges
                               15                                 3
                        character
                                1
>
```

For the 6 BioPaxLevel2 files the rdataclass is:

```{r}
> sort(table(mcols(biopax2)$rdataclass), decreasing = TRUE)
biopax
     6
>
```

For the 4 JSON files the rdataclass is:

```{r}
> sort(table(mcols(json)$rdataclass), decreasing = TRUE)
data.table
         4
>
```

For the 3 BioPax files the rdataclass is:

```{r}
> sort(table(mcols(biopax)$rdataclass), decreasing = TRUE)
biopax
     3
>
```

For the 2 HDF5 files the rdataclass is:

```{r}
> sort(table(mcols(hdf5)$rdataclass), decreasing = TRUE)
String
     2
>
```

For the 1 GRASP file the rdataclass is:

```{r}
> sort(table(mcols(grasp)$rdataclass), decreasing = TRUE)
SQLiteConnection
               1
>
```

For the 1 mzid file the rdataclass is:

```{r}
> sort(table(mcols(mzid)$rdataclass), decreasing = TRUE)
mzRident
       1
>
```

For the 1 mzML file the rdataclass is:

```{r}
> sort(table(mcols(mzml)$rdataclass), decreasing = TRUE)
mzRpwiz
      1
>
```

For the 1 mzTab file the rdataclass is:

```{r}
> sort(table(mcols(mztab)$rdataclass), decreasing = TRUE)
MSnSet
     1
>
```

For the 1 RDA file the rdataclass is:

```{r}
> sort(table(mcols(rda)$rdataclass), decreasing = TRUE)
list
   1
>
```

For the 1 tab file the rdataclass is:

```{r}
> sort(table(mcols(tab)$rdataclass), decreasing = TRUE)
data.frame
         1
>
```

For the 1 XLS/XLSX file the rdataclass is:

```{r}
> sort(table(mcols(xls)$rdataclass), decreasing = TRUE)
character
        1
>
```
