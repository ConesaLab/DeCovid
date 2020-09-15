# DeCovid: A tool for conducting a Wide-spread transcriptional analysis for COVID-19 Disease Map genes between males and females.

- Tianyuan Liu
- Ana Lleo
- Leandro Balzano Nogueira
- Ana Conesa

### About DeCovid
This app provides information on gene expression differences between man and women and old versus young individuals for genes in the COVID19 Disease Map. Gene expression data was obtained from the GTEx project, that collects expression data for multiple human tissues.
- You can select tissue of interest, companions type (man vs female or young vs old), a p.value and logFC of differential expression to identify regulated COVID-19 genes.
- Also, you can simply introduce a gene name to see the expressions data across these population groups.
- Finally, you can implement a GO term enrichment analysis for the DE genes.

## Video Tutorial
[![Watch the video](https://github.com/TianYuan-Liu/Basic-app-COVID-19-/blob/master/www/video.png)](https://www.youtube.com/watch?v=vwtjcw7cJG4)

## Introduction
Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) can cause severe medical complication and lead patients to die. Based on current data, males have a higher risk than females, which means that males are more likely to die when they are diagnosed as cases. Based on the expression data from the Genotype-Tissue Expression (GTEx) project and genes related to COVID-19  in COVID-19 Disease Map, we wrote DeCovid,  a shiny app using edgeR R package to explore baseline gene expression difference and understand the sex bias.


## How to Install and Run DeCovid

There are two main options to run MirCure:

1.  Run DeCovid as *Docker image* 
2.  Directly download and run DeCovid as *shiny app*.

### Run  DeCovid as Docker Image

- install the [Docker engine](https://docs.docker.com/engine/install/).
- Run DeCovid with the following command:
```
docker run --rm -p 3838:3838 tianyuanliu/explore_covid-19_disease_genes
```
- Open the URL shows in the terminal (typically [http://[::]:3838](http://[::]:3838)) in any **web browser**.

### Run  DeCovid as shiny app
- Download the folder
- Unzip the file in data folder
- Install all the packages listed in the script
- Run the app
