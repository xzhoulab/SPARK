---
#
# By default, content added below the "---" mark will appear in the home page
# between the top bar and the list of recent posts.
# To change the home page layout, edit the _layouts/home.html file.
# See: https://jekyllrb.com/docs/themes/#overriding-theme-defaults
#
layout: home
author: "Shiquan Sun, Jiaqiang Zhu and Xiang Zhou"
date: '2019-06-08'
fontsize: 10pt
linkcolor: blue
mainfont: Times
output:
pdf_document: default
html_document: default
md_document:
variant: markdown_github
---

![SPARK\_pipeline](pipline.png)

## SPARK

**SPARK** is an efficient method to identify genes with spatial expression pattern. The intended applications are spatially resolved RNA-sequencing from e.g., Spatial Transcriptomics, or *in situ* gene expression measurements from
e.g., SeqFISH, or Merfish.


## Example: Breast Cancer Data

Load the `SPARK` package and Breast cancer data set.
```R
    library('SPARK')
    load("~/data/Layer2_BC_Count.rds")
     
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
Create a SPARK object for analysis. This step excludes the gene that are lowly expressed in each gene.
```R 
    ## filter genes and cells/spots and 
    spark <- CreateSPARKObject(counts=rawcount, 
                                 location=info[,1:2],
                                 prectage = 0.1, 
                                 min_total_counts = 10)

    ## total counts for each cell/spot
    spark@lib_size <- apply(spark@counts, 2, sum)

    ## Take the first ten genes as an example
    spark@counts   <- spark@counts[1:10,]
```

Fit the statistical model under the null hypothesis.
```R 
    ## Estimating Parameter Under Null
    spark <- spark.vc(spark, 
                       covariates = NULL, 
                       lib_size = spark@lib_size, 
                       num_core = 1,
                       verbose = F)
```

Test the spatially expressed pattern genes. By default, the kernel matrices are computed automatically, and check the positive definition of the kernel matrices. There is also an option to provide a kernel matrix by user.
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
GAPDH   7.477171e-09    1.683453e-07
MAPKAPK2    1.016078e-01    5.952118e-01
MCL1    1.149519e-08    1.683453e-07
TMEM109 4.303998e-01    1.000000e+00
TMEM189 6.189064e-01    1.000000e+00
ITPK1   7.213287e-01    1.000000e+00
```
