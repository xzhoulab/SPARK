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
output:
pdf_document: default
html_document: default
md_document:
variant: markdown_github
---

![SPARK\_pipeline](MethodOverview.png)

## SPRINT

**SPRINT** is an efficient method to identify genes with spatial expression pattern. 
The intended applications are spatially resolved RNA-sequencing from e.g.
Spatial Transcriptomics, or *in situ* gene expression measurements from
e.g. SeqFISH, Merfish.

## System requirements
SPRINT has been tested on R 3.3.1 and is platform independent (tested on Linux, OS X and Windows)
Installation can then be done via the devtools package:

```R
library('devtools')
devtools::install_github('xzhoulab/SPRINT')
```
Alternatively, installation can then be done from a local binary package from the shell:
```bash
R CMD INSTALL SPRINT_1.0.0.tar.gz
```



## Sample Code: Analysis of Breast Cancer Data
```R
    library('SPRINT')
    load("~/data/Layer2_BC_Count.rds")
     
    ## rawcount matrix of genes by cells/spots
    rawcount[1:5,1:5]
```

|               | 17.907x4.967  | 18.965x5.003 | 18.954x5.995|17.846x5.993  |20.016x6.019  |
| ------------- | ------------- |------------- |-------------|------------- |------------- |
|GAPDH          |     1         |   7          |  5          |  1           |      2       |
|USP4           |     1         |   0          |  0          |  0           |      0       |
|MAPKAPK2       |     1         |   1          |  0          |  0           |      1       |
|CPEB1          |     0         |   0          |  0          |  0           |      0       |
|LANCL2         |     0         |   0          |  0          |  0           |      0       |
```R    
    ## extract the coordinates from the rawdata
    info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",1)),
                            y=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",2)),
                            total_counts=apply(rawcount,2,sum))
    rownames(info) <- colnames(rawcount)

    ## filter genes and cells/spots and create the SPRINT object for following analysis
    sprint <- CreateSPRINTObject(counts=rawcount, location=info[,1:2],
                                 prectage = 0.1, 
                                 min_total_counts = 10)

    ## total counts for each cell/spot
    sprint@lib_size <- apply(sprint@counts, 2, sum)

    ## Take the first ten genes as an example
    sprint@counts   <- sprint@counts[1:10,]

    ## Estimating Parameter Under Null
    sprint <- sprint.vc(sprint,covariates = NULL, lib_size=sprint@lib_size, num_core=1,verbose=F)

    ## Calculating pval
    sprint <- sprint.test(sprint, check_positive = T, verbose=F)
    
    ## Check pvals 
    head(sprint@res_mtest[,c("combined_pvalue","adjusted_pvalue")])
```

|               | combined_pvalue | adjusted_pvalue |
| ------------- | ------------- |------------- |
|GAPDH       |7.477171e-09   | 1.683453e-07 |
|MAPKAPK2    |1.016078e-01   | 5.952118e-01 |
|MCL1        |1.149519e-08   | 1.683453e-07 |
|TMEM109     |4.303998e-01   | 1.000000e+00 |
|TMEM189     |6.189064e-01   | 1.000000e+00 |
|ITPK1       |7.213287e-01   | 1.000000e+00 |
