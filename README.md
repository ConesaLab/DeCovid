# DeCovid: A tool for conducting a Wide-spread transcriptional analysis for COVID-19 Disease Map genes

- Tianyuan Liu
- Leandro Balzano Nogueira
- Ana Lleo
- Ana Conesa

### About DeCovid
This app provides information on gene expression differences between man and women and old versus young individuals for genes in the COVID-19 Disease Map. Gene expression data was obtained from the GTEx project, that collects expression data for multiple human tissues.
- You can select tissue of interest, companions type (man vs female or young vs. old), a p.value and logFC of differential expression to identify regulated COVID-19 genes.
- Also, you can simply introduce a gene name to see the expressions data across these population groups.
- Finally, you can implement a GO term enrichment analysis for the DE genes.

## Video Tutorial
[![Watch the video](https://github.com/ConesaLab/DeCovid/blob/master/www/AA65B516-2B6A-463E-AEA5-2C9D7FD4C2D2.jpeg)](https://youtu.be/aBwrSgVLSqQ)

## Introduction
Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) can cause severe medical complication and lead patients to die. Based on current data, males have a higher risk than females, which means that males are more likely to die when they are diagnosed as cases. Based on the expression data from the Genotype-Tissue Expression (GTEx) project and genes related to COVID-19  in COVID-19 Disease Map, we wrote DeCovid,  a shiny app using edgeR R package to explore baseline gene expression difference and understand the sex bias.
<img src="https://github.com/ConesaLab/DeCovid/blob/master/www/idea.png">

## Method
There are two main **Datasets** integrated with DeCovid:
- RNA-seq data: Genotype-Tissue Expression project (GTEx) (Lonsdale *et al.,* 2013).
- COVID-19 Disease Map: A list of genes mined 9996 PMC papers in COVID-19 Open Research Dataset (CORD-19) using machine learning approaches(Lu Wang *et al.,* 2020).

The **DeCovid software** is a Shiny app written in R with a user-friendly interface. DeCovid already has integrated all essential data from GTEx, and no additional downloads are necessary. Differential gene expression is calculated using edgeR (Robinson *et al.,* 2010), and multiple testing correction is applied following the Benjamini and Hochberg（BH）method (Benjamini and Hochberg, 1995). Results are provided either as lists of differentially expressed genes, heatmaps of sex and age mean expression values, and gene-specific boxplots showing the distribution of expression values across demographic groups. GO enrichment analysis of significant gene sets is provided through the cluster profile R package and uses the COVID-19 disease map list as a reference set (Yu *et al.,* 2012).


## How to Install and Run DeCovid

There are two main options to run DeCovid:

1.  Run DeCovid as *Docker image* 
2.  Directly download and run DeCovid as *shiny app*.

### Run  DeCovid as Docker Image

- install the [Docker engine](https://docs.docker.com/engine/install/).
- Run DeCovid with the following command in terminal (Mac/Linux) or PowerShell (Win):
```
docker run --rm -p 3838:3838 conesalab/decovid:decovid
```
- Open the URL shows in the terminal (typically [http://[::]:3838](http://[::]:3838)) in any **web browser**.

### Run  DeCovid as shiny app
- Download the folder
- Unzip the file in data folder
- Install all the packages listed in the script
- Run the app

## Reference
- Lonsdale,J. et al. (2013) The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45, 580–585.
- Lu Wang,L. et al. (2020) CORD-19: The Covid-19 Open Research Dataset. ArXiv.
- Robinson,M.D. et al. (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26, 139–140.
- Benjamini,Y. and Hochberg,Y. (1995) Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society. Series B (Methodological), 57, 289–300.
- Yu,G. et al. (2012) clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. OMICS, 16, 284–287.

