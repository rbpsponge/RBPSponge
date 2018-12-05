# RBPSponge

RBPSponge is a web tool to identify potential sponge lncRNAs for RBPs (www.rbpsponge.com). Here we provide the scripts and some example data sets to perform the statistical analyses available on the RBPSponge web tool.

## Setup

### Requirements
1. R (3.3.3 )
* ggplot2(2.2.1)
* data.table (1.10.4-3)
* DAAG (1.22)
* lmtest (0.9-35)

### Usage

There are four  R scripts, please make sure EMTAB2706_mapping_frame.RData, EMTAB2706_tpm.RData and EMTAB2706_median_tpm.RData are in the same directory as the scripts.  

* plot_expression.R plots the expression of a given lncRNA and a RBP. 	This script can be run as follows:

```	
Rscript plot_expression.R <lncRNA_id>  <RBP_id> <dataset_id> <job_id>
```

Here job_id is used to specify the output filenames. 

Here is an example run:

``` 
Rscript plot_expression.R ENSG00000260032 PUM2 EMTAB2706 1
```

* corr_analysis.R calculates the Spearman correlation between the expression of target/background genes of an RBP and the expression of the lncRNA. Then a box plot is generated where the distribution of correlation values are compared between target genes and background genes. The plot also displays the p-value obtained with Wilcoxon ranksum test. This script can be run as follows:
``` 
Rscript corr_analysis.R <lncRNA_id>  <RBP_id> <dataset_id> <job_id>
``` 
Here is an example run:
``` 
Rscript corr_analysis.R ENSG00000260032 PUM2 EMTAB2706 1
``` 
	
* regress_analysis.R assesses whether  lncRNA expression has added predictive value in addition to RBP expression in determining the expression levels of target genes. A simple linear regression analysis is performed to predict target gene expression where in one case only RBP expression is used and in the other case both RBP and lncRNA expression are used as features. The Spearman correlation coefficient values between actual and predicted expression values of target genes are calculated on held-out datasets using 10-fold cross-validation. If the lncRNA of interest acts as a sponge for the RBP, we expect to see an improved predictive performance when lncRNA expression is included as an additional feature. The significance of change is evaluated with Wilcoxon rank-sum test and likelihood ratio test. The script can be run as follows:

```
Rscript regress_analysis.R <lncRNA_id>  <RBP_id> <dataset_id> <job_id>
```

Here is an example run:

```
Rscript regress_analysis.R ENSG00000260032 PUM2 EMTAB2706 1
```

* knockdown_analysis. R investigates the expression changes upon the  knockdown of the lncRNA of interest, when available. Because lncRNA activity is minimized we expect to see an increased RBP activity and a more pronounced effect (either stabilizing or de-stabilizing) on the expression of target genes. To this end, we compare the distribution of expression changes of target and background genes with a cumulative distribution frequency (CDF) plot.  Wilcoxon rank-sum test is used to assess the significance between the two distributions. The script can be run as follows:

```
Rscript knockdown_analysis.R <lncRNA_id>  <RBP_id>  <job_id>
```

Here is an example run:
```
Rscript knockdown_analysis.R ENSG00000260032 PUM2 1
```

--------------

[EMTAB-2706](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2706/) RNA-seq of 675 commonly used human cancer cell lines
