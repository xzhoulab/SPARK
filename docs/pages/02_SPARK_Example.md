---
layout: page
title: Example Analaysis with SPARK
permalink: /02_SPARK_Example/
---

## Example Analysis with SPARK: Breast Cancer Data

Load the `SPARK` package and Breast cancer data set, which can be downloaded [here](https://github.com/xzhoulab/SPARK/blob/master/data/Layer2_BC_Count.rds).
```R
    library('SPARK')
    load("./Layer2_BC_Count.rds")
     
```
View the expression count matrix `rawcount`, each row denotes a gene and each column represents a cell/spot.
```R
rawcount[1:5,1:5]

    17.907x4.967   18.965x5.003   18.954x5.995    17.846x5.993 20.016x6.019
GAPDH   1   7   5   1   2
USP4    1   0   0   0   0
MAPKAPK2    1   1   0   0   1
CPEB1   0   0   0   0   0
LANCL2  0   0   0   0   0
```

Extract the annotation information for each sample, i.e., location or coordinates
```R   
    ## extract the coordinates from the rawdata
    info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",1)),
                             y=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",2)),
                             total_counts=apply(rawcount,2,sum))
    rownames(info) <- colnames(rawcount)
```
Create a SPARK object for analysis. This step excludes the gene that are lowly expressed.
```R 
    ## filter genes and cells/spots and 
    spark <- CreateSPARKObject(counts=rawcount, 
                                 location=info[,1:2],
                                 percentage = 0.1, 
                                 min_total_counts = 10)

    ## total counts for each cell/spot
    spark@lib_size <- apply(spark@counts, 2, sum)

    ## Take the first ten genes as an example
    #spark@counts   <- spark@counts[1:10,]
```

Fit the statistical model under the null hypothesis.
```R 
    ## Estimating Parameter Under Null
    spark <- spark.vc(spark, 
                       covariates = NULL, 
                       lib_size = spark@lib_size, 
                       num_core = 5,
                       verbose = F)
```

Test the spatially expressed pattern genes. By default, the kernel matrices are computed automatically by coordinates, and check the positive definition of the kernel matrices. There is also an option to provide a kernel matrix by user.
```R 
    ## Calculating pval
    spark <- spark.test(spark, 
                         check_positive = T, 
                         verbose = F)
    
```

Output the final results, i.e., combined p-values, adjusted p-values, etc. 
```R 
head(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")])

    combined_pvalue   adjusted_pvalue
GAPDH   7.477171e-09    4.233403e-06
MAPKAPK2    1.016078e-01    1.000000e+00
MCL1    1.149519e-08    6.078909e-06
TMEM109 4.303998e-01    1.000000e+00
TMEM189 6.189064e-01    1.000000e+00
ITPK1   7.213287e-01    1.000000e+00
```
