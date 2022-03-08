# bryozoan-data-paper

This repository contains data and code for the paper "Regression, transformations, and mixed-effects with marine bryozoans" (Evans, 2022), which discusses how to use data from Pettersen *et al.* (2015) to teach statistics students about data wrangling, regression, transformations, and mixed effects models.

## Files

### Data

* bryozoan_raw.csv : The original data set made available by Pettersen *et al.* It contains an error in the recorded mass for late-stage *Bugula neritina*
* bryozoan_data.csv : A tidy version of bryozoan_raw.csv, with the untransformed mass and metabolic rate. The error for late-stage *Bugula* is still present
* bryozoan_data_fixed.csv : The fully cleaned data, with the error fixed

### Analysis

* bryozoan_paper_analysis.R : An R script containing all analysis used for data processing and analysis
* `*`.pdf : images for data visualization and modeling in the paper

### Activity

* bryozoan_paper_activity.html : An activity that walks students through an initial reading of Pettersen *et al.* (2015). This activity is also hosted at [https://ciaran-evans.github.io/files/bryozoan_paper_activity.html](https://ciaran-evans.github.io/files/bryozoan_paper_activity.html)
* bryozoan_paper_activity.Rmd : The source file for the activity


## References

Pettersen, A. K., White, C. R., & Marshall, D. J. (2015). Why does offspring size affect performance? Integrating metabolic scaling with life-history theory. *Proceedings of the Royal Society B: Biological Sciences*, 282(1819), 20151946. URL: [https://royalsocietypublishing.org/doi/10.1098/rspb.2015.1946](https://royalsocietypublishing.org/doi/10.1098/rspb.2015.1946)