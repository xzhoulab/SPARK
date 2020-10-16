---
layout: page
title: Experiments
permalink: /03_experiments/
---
### SPARK analysis
#### Data sets for experiments
  * [All raw data](https://github.com/xzhoulab/SPARK-Analysis/tree/master/raw_data)
  * [All processed data](https://github.com/xzhoulab/SPARK-Analysis/tree/master/processed_data)
  
#### Code for experiments
  * [Simulations](https://github.com/xzhoulab/SPARK-Analysis/tree/master/simulation)
  * [Real data analysis](https://github.com/xzhoulab/SPARK-Analysis/tree/master/analysis)
  
#### Example
```R
rm(list = ls())
source("./funcs/funcs.R")
load("./output/MOB_Pattern_SpatialDE.rds")

tmplist <- datlist
for (i in 1:5) {
    tmplist[[i]][, 2] <- relative_func(datlist[[i]][, 2])
}

patterns = c("I", "II", "III")
## three major pattern were used for simulation
df <- setNames(cbind.data.frame(tmplist[[1]][, 1], 
                do.call(cbind, sapply(tmplist[c(5, 1, 4)], "[", 2))), 
                c("xy", paste0("Pattern ", patterns)))
pp <- lapply(1:3, function(x) {
    pattern_plot2(df, x, xy = F, main = T, titlesize = 1.5)
})


grid.arrange(grobs = pp, ncol = 3)

```
![summarized patterns from mouse olfactory bulb](mouseOB_pattern.png)


### SPARK-X analysis
#### Data sets for experiments
* Raw data can be downloaded through the link listed in the manuscript
* [All processed data](https://github.com/xzhoulab/SPARK-X-Analysis/tree/main/processed_data)

#### Code for experiments
* [Simulations](https://github.com/xzhoulab/SPARK-X-Analysis/tree/main/simulation) with [example](http://htmlpreview.github.io/?https://github.com/xzhoulab/SPARK-X-Analysis/blob/main/simulation/hotspot_simulation_example.html)
* Real data analysis
  * [HDST Analysis with SPARK-X](http://htmlpreview.github.io/?https://github.com/xzhoulab/SPARK-X-Analysis/blob/main/analysis/HDST_SPARKX.html)
  * [Slide-seq Analysis with SPARK-X](http://htmlpreview.github.io/?https://github.com/xzhoulab/SPARK-X-Analysis/blob/main/analysis/SV1_SPARKX.html)
  * [Slide-seqV2 Analysis with SPARK-X](http://htmlpreview.github.io/?https://github.com/xzhoulab/SPARK-X-Analysis/blob/main/analysis/SV2_SPARKX.html)

  