---
layout: page
title: Example Analaysis
permalink: /02_SPARK_Example/
---

### Example Analysis with SPARK: Breast Cancer Data

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



### Example Analysis with SPARK-X: HDST Data

Load the `SPARK` package and HDST data set, which can be downloaded [here](https://github.com/xzhoulab/SPARK-X-Analysis/blob/main/processed_data/CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds).
```R
    library('SPARK')
    load("./CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds")
     
```
View the expression count matrix `sp_count`, each row denotes a gene and each column represents a spot.
```R
sp_count[26:30,1:5]

        1000x100 1000x103 1000x113 1000x114 1000x116
Gm42418        1        .        2        .        .
Gm10925        .        .        .        .        .
Gm7135         .        .        .        .        .
Atrx           .        .        .        .        .
Celf2          .        .        .        .        .
```

View the coordinates information `location`, each row denotes a gene and each column represents a spot.
```R   
## extract the coordinates from the rawdata
info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(sp_count),split="x"),"[",1)),
                         y=as.numeric(sapply(strsplit(colnames(sp_count),split="x"),"[",2)))
rownames(info)  <- colnames(sp_count)
location        <- as.matrix(info)
```

Removing mitochondrial genes
```R 
mt_idx      <- grep("mt-",rownames(sp_count))
if(length(mt_idx)!=0){
    sp_count    <- sp_count[-mt_idx,]
}
```

Analyze the data with SPARK-X
```R 
sparkX <- sparkx(sp_count,location,numCores=1,option="mixture")
## ===== SPARK-X INPUT INFORMATION ==== 
## number of total samples: 177455 
## number of total genes: 19913 
## Running with single core, may take some time 
## Testing With Projection Kernel
## Testing With Gaussian Kernel 1
## Testing With Gaussian Kernel 2
## Testing With Gaussian Kernel 3
## Testing With Gaussian Kernel 4
## Testing With Gaussian Kernel 5
## Testing With Cosine Kernel 1
## Testing With Cosine Kernel 2
## Testing With Cosine Kernel 3
## Testing With Cosine Kernel 4
## Testing With Cosine Kernel 5
```

Output the final results, i.e., combined p-values, adjusted p-values, etc. 
```R 
head(sparkX$res_mtest)
        combinedPval adjustedPval
Rcn2      0.11886766            1
Mycbp2    0.03704549            1
Mprip     0.10150035            1
Mroh1     0.92949857            1
Zfp560    0.28341540            1
Cacna1e   0.08444910            1
```

