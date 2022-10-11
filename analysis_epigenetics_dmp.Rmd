---
title: "Example - Epigenetic Differential Expression, DMP"
author: "Brett Maroni-Rana"
output:
  html_document:
    df_print: paged
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(
	echo = FALSE,
	warning = FALSE,
	collapse = TRUE,
	comment = "#>",
	out.width = "100%"
)
str(knitr::opts_chunk$get())
```

```{r library, include = FALSE}
# remove all objects
#rm(list = ls())
# load libraries
library("tidyverse")
library('minfi')
library('sva')
library('limma')
library('ggplot2')
library('GenomicRanges')
library('TxDb.Hsapiens.UCSC.hg19.knownGene') # confirm hg19
library('DMRcate')
library('RPMM')
library('ComplexHeatmap')
library('circlize')
library('RColorBrewer')
library('GenomicFeatures')
library('Gviz')
```

```{r wd, echo=TRUE}
# set working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
```
### Load Data
Raw data has been processed previously and saved as a quantile normalized, probe filtered, genomic annotated ratio set. Here, the object loaded is __gr_filt.rda__.
```{r load gr_filt}
# load normalized and filtered genomic ratio set
load('/Volumes/CRDShared/GradyLab/Brett Maroni-Rana/methylation/radiogenomics/preprocessing/20211215_init/gr_filt.rda')
```

## Data Exploration
### QC check : PCA most variable probes
Prior to batch corrections and further analysis, the most variable probes according to m-value standard deviation can visualized with PCA. This quick QC check assess distribution of the two most important aspects of the data, principal component 1 and principal component 2. This looks at covariation and can determine if covariates need to be accounted for if the distribution is not random. Batch correction and SVA can account for covariation.
```{r pca code, include=FALSE}
# PCA based on 1000 most variable (sd) cpg probes
beta <- getBeta(gr_filt)
m_value <- logit2(beta)

# find sd, sort sd, subset most variable probes
row_sd <- apply(m_value, 1, sd)
row_sd_sort <- sort(row_sd, decreasing = T)
mvp <- names(row_sd_sort)[1:1000] #1000

# run pca on mvps
pca_mvp <- prcomp(t(m_value[mvp,]))$x
```

```{r pca plot code, include=FALSE}
# plot batch
p_batch <- ggplot(data.frame(PC1=pca_mvp[,1], PC2=pca_mvp[,2], Batch=as.factor(gr_filt$batch)), aes(x = PC1, y = PC2, color=Batch)) +
  geom_point(size=3, alpha=0.8) +
  labs(title = "PCA of 1,000 MVPs: Batch")

# plot location
p_location <- ggplot(data.frame(PC1=pca_mvp[,1], PC2=pca_mvp[,2], Location=as.factor(gr_filt$location)), aes(x = PC1, y = PC2, color = Location)) +
  geom_point(size=3, alpha=0.8) +
  labs(title = "PCA of 1,000 MVPs: Location")

# plot sex
p_gender <- ggplot(data.frame(PC1=pca_mvp[,1], PC2=pca_mvp[,2], Gender = as.factor(gr_filt$sex)), aes(PC1, PC2, color = Gender)) +
  geom_point(size=3, alpha=0.8) +
  labs(title = "PCA of 1,000 MVPs: Sex")

# potential progressor
p_prog <- ggplot(data.frame(PC1=pca_mvp[,1], PC2=pca_mvp[,2], Prog = gr_filt$progression), aes(PC1, PC2, color = Prog)) +
  geom_point(size=3, alpha=0.8) +
  labs(title = "PCA of 1000 MVPs: Progression")
```

```{r pca plot}
p_batch
p_location
p_gender
p_prog
```
```{r, include=FALSE}
# pca pdf
#pdf("pca_phenotypes.pdf", height=11, width=8.5)
#cowplot::plot_grid(p_path, p_batch, p_location, p_gender, p_prog, p_vol_year, ncol=2)
#dev.off()
```

PCA analysis on the top 1000 most variable probes by standard deviation shows a few interesting patterns.
  
- __Batch__: Significant batch effect. Batch 1 and 2 seem random. Batch 3 concentrates toward negative PC2. There is a need for batch correction. 
- __Location__: Distinct distribution of PCs between left and right sided colon samples. This is a biological affect and should be preserved, although analyzing these samples separately according to location may aid in analysis. 
- __Gender__: Seemingly randomly distributed
- __Progression__: High progressor is randomly distributed. Progressor and Static have unique patterns, P more than S. This indicates the methylation signitures are different between the three groups and warrents further investigation

### Goals of Study
Following PCA, how can the samples best be characterized? and what are the goals of this study?
  
1. What are the differentially expressed CpGs between potential progressor phenotypes?
2. Of those DMPs, are they found close together in differentially expressed regions?
3. What are the characterisitics DMPs?
  
### Combat : Batch Effect Removal  
  
Accounting for and removing known and unknown variability in an experiment can reduce error rates. To do so, two methods are used and can be found within the __sva__ R package. The first of these functions is _ComBat()_ in which an empirical Bayesian framework is applied to the dataset. Once batch effects are removed, unknown variability can be removed with _sva()_ or further processed with downstream analysis. 
  
comBat code and visualization is below:
```{r, echo=TRUE}
# set model matrix
mod = model.matrix(~ gr_filt$progression +
                     gr_filt$location)
# run batch correction
m_value_combat <- ComBat(dat = m_value, 
                         mod = mod,
                         batch = gr_filt$batch)

# get beta
beta_combat <- ilogit2(m_value_combat)
```

```{r, include=TRUE}
# PCA based on 1000 most variable (sd) cpg probes
# find sd, sort sd, subset most variable probes
row_sd <- apply(m_value_combat, 1, sd)
row_sd_sort <- sort(row_sd, decreasing = T)
mvp <- names(row_sd_sort)[1:1000] #1000

# run pca on mvps
pca_mvp_combat <- prcomp(t(m_value_combat[mvp,]))$x
```

```{r, include= FALSE}
# create PCA df for plotting
df1 <- data.frame(PC1 = pca_mvp_combat[,1], PC2 = pca_mvp_combat[,2], Batch = gr_filt$batch)
df2 <- data.frame(PC1 = pca_mvp_combat[,1], PC2 = pca_mvp_combat[,2], Sex = gr_filt$sex)
df3 <- data.frame(PC1 = pca_mvp_combat[,1], PC2 = pca_mvp_combat[,2], Location = gr_filt$location)
df4 <- data.frame(PC1 = pca_mvp_combat[,1], PC2 = pca_mvp_combat[,2], Prog = gr_filt$progression)
```

```{r, include = FALSE}
# plot pca combat batch
p_batch_combat <- ggplot(df1, aes(x = PC1, y = PC2, color=Batch)) +
  geom_point(size=3, alpha=0.8) +
  labs(title = "PCA of 1,000 MVPs: Batch after ComBat")

# plot pca combat path
p_gender_combat <- ggplot(df2, aes(x = PC1, y = PC2, color=Sex)) +
  geom_point(size=3, alpha=0.8) +
  labs(title = "PCA of 1,000 MVPs: Sex after ComBat")

# plot pca combat location
p_location_combat <- ggplot(df3, aes(x = PC1, y = PC2, color=Location)) +
  geom_point(size=3, alpha=0.8) +
  labs(title = "PCA of 1,000 MVPs: Location after ComBat")

# plot pca combat location
p_prog_combat <- ggplot(df4, aes(x = PC1, y = PC2, color=Prog)) +
  geom_point(size=3, alpha=0.8) +
  labs(title = "PCA of 1,000 MVPs: Progression after ComBat")
```


```{r}
p_batch_combat
p_gender_combat
p_location_combat
p_prog_combat
```
```{r, include=FALSE}
# plot combat adjust
#pdf("pca_phenotypes_combat.pdf", height=11, width=11)
#cowplot::plot_grid(p_batch_combat, p_path_combat, p_location_combat, p_prog_combat, nrow = 2)
#dev.off()
```

Following batch correction and the PCA plotting of 1,000 most variable probes, batch 3 looks to be more randomly distributed. Before batch correction, batch 3 is slightly more distributed toward a negative PC1 and PC2. After batch correction, batch 3 is slightly more distributed toward a positive PC2 resulting in a slight more evenly distributed batch.  
  
Importantly, Progression and Location are relatively unchanged with their distinct patterns remaining. 

### Surrogate Variable Analysis (SVA) 
  
Surrogate Variable Analysis removes unknown variability in a dataset by first estimating surrogate variables, then by including them into a model. In this study, the surrogate variables are used in the limma package along with the variables from the data.   
  
Importantly, this study focuses on methylation expression vs potential progression, without interest in demographics. If demographics were integral to this analysis, sva may not be appropriate in an effort to preserve the integrity of the demographic data. Polyp location is exception of this and is not included in SVA model.  
  
Covariates taken into consideration and correcting for include final age and sex. Additional phenotypic data may be appropriate.
  
sva code is below:  
```{r, echo=TRUE}
# sva for surrogate variables
mod <- model.matrix(~ gr_filt$progression +
                      gr_filt$age_t0 +
                      gr_filt$sex)

mod0 <- model.matrix(~ gr_filt$age_t0 +
                       gr_filt$sex)

svobj <- sva(dat = m_value_combat, mod = mod, mod0 = mod0) # 17 variables
#save(svobj, file='svobj.rda')
```
### EWAS : Linear Models for Microarray Data (limma) 
  
Next, a linear model can be fit to our data using the limma package. An experiment is designed with model.matrix and used to distinguish phenotypic differences between samples. This experimental design and methylation data matrix is used by lmFit to fit linear models to each row of the data. Each row in our methylation matrix corresponds to a CpG probe. Finally, a moderated t-statistic and log-odds of differential expression is calculated with eBayes and summarized using topTable. Differntially methylated probes can be annotated downstream or visualized together in a heatmap. See below:
```{r, echo=T}
# EWAS DMP analysis with limma
# single comparison
fit <- eBayes(lmFit(m_value_combat, 
                    model.matrix(~ gr_filt$progression +
                                   gr_filt$location +
                                   svobj$sv)), robust = T)

res <- topTable(fit, coef = 2, number = nrow(m_value_combat))
res <- res[,c(1,4,5)]
colnames(res) <- c("logFC", "pval","adj_pval")
```

__Top CpGs__  
Criteria: adjusted p-value < 0.05 and absolute LogFC >= 1.5
```{r, echo=T}
# limma summary
top_cpgs <- sum(res$adj_pval < 0.05 & abs(res$logFC) > 1.5)
cat("Top CpGs:", top_cpgs)
```
Following limma, there are __62__ CpG probes with an adjusted P-Value of < 0.01 and an absolute logFC greater than 2.

__Annotation__  
Each CpG probe has accompanying annotatin data from Illumina. Minfi package is used to retrieve the annotation for downstream anlaysis.
```{r, include=FALSE}
# add probe annotation
probeInfo <- as.data.frame(minfi::getAnnotation(gr_filt))[,c(1,2,3,18,19,22:24)]
tmp.name <- strsplit(as.character(probeInfo$UCSC_RefGene_Name), ';')
tmp.name <- lapply(tmp.name, FUN = unique)
tmp.name <- lapply(tmp.name, paste, collapse=';')
probeInfo$UCSC_RefGene_Name <- unlist(tmp.name)
tmp.group <- strsplit(as.character(probeInfo$UCSC_RefGene_Group), ';')
tmp.group <- lapply(tmp.group, FUN = unique)
tmp.group <- lapply(tmp.group, paste, collapse=';')
probeInfo$UCSC_RefGene_Group <- unlist(tmp.group)
res <- cbind(probeInfo[rownames(res),], res)
```

DMPs are annotated and written as a tab-delimited text file. __See "result_dmp.txt"__
_example RMD does not write to file_
```{r, echo=T}
# example - dmp annotation
head(res)

# write to .csv
#write.table(res, 'result_dmp.txt', sep='\t', quote = F)

# save to .rda
#save(res, file='res.rda')
```

## Visualization
### Heatmap with Top CpGs
```{r, include=FALSE}
# remove objects and load packages
library(gplots)
```

__Get Top CpGs__  
Following DMP analysis with Limma, we can visualize the results in a heatmap. We can get the most statistically relevant results by filtering our results with an FDR adjusted p-value < 0.05 and log2FC > +-1.5
```{r, include=TRUE}
## top cpgs tubvil vs tub
res_sort <- res %>%
  filter(adj_pval < 0.05 & abs(logFC) > 1.5) %>%
  arrange(adj_pval & abs(logFC))
head(res_sort)
top_probes <- rownames(res_sort)
```

__Generate Heatmap__  
The heatmap package _pheatmap_ is used to generate a heatmap and includes a hierarchical clustering method. Here, default settings of _hclust()_ are used. Euclidean distance is measured between samples, based on thier methylation value, while complete linkage clustering option is used to generate a dendogram.
```{r, include=FALSE}
# heatmap - pheatmap
library(pheatmap)
# create df for categorical variables
df <- data.frame(sex = gr_filt$sex,
                 group = gr_filt$progression,
                 location = gr_filt$location)

# match names to array_id
rownames(df) <- colnames(ilogit2(m_value_combat)[top_probes, colnames(gr_filt)])
# make heatmap
# large growing vs small not growing
h1 <- pheatmap(ilogit2(m_value_combat)[top_probes, colnames(gr_filt)],
               show_colnames = F,
               show_rownames = F,
               cluster_cols = T,
               cluster_rows = T,
               cutree_rows = 1,
               cutree_cols = 5,
               annotation_col = df,
               scale = 'none',
               main = 'EWAS DMP Analysis, Potential Progressor Characterization\n62 CpGs Adj P-Val < 0.05 & Log2FC > +-1.5')

#pdf('heatmap_limma_lgsg.pdf')
#h1
#dev.off()
```

```{r}
# draw heat map
h1

# save to pdf
#pdf("heat_map_probes.pdf")
#draw(ht1)
#dev.off()
```

Following dmp analysis, 62 differentially methylated probes are found. Above, they are mapped to a heatmap according to demographic data. Visually, the 'H group' does populate to the left side of the heatmap but the probes only show marginal differences between group. This visual pattern can be tested statistically, but a deeper look at the potential progressor categorical variable should also be undertaken. Some samples may be misclassified.
