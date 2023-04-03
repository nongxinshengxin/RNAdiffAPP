# RNAdiffAPP
[![docker](https://img.shields.io/badge/Docker_image-ready-red.svg)](https://hub.docker.com/r/mouther/rnadiff) [![source](https://img.shields.io/badge/Source_code-support-blue.svg)](https://github.com/nongxinshengxin/RNAdiffAPP/tree/main/R)
## Introduction
### What is RNAdiff App ?
This app is used to perform **RNA-Seq downstream analysis** and is available to the user through an interactive interface.
### How is this RNAdiff App built ?
The app is built on the R language and shiny packages, calling DESeq2, edgeR, ggplot2 and other R packages.
### What can the RNAdiff App do ?
The app is divided into five main sections.<br/>
- *__Function 1:__ Analysis of differentially expressed genes, which can be performed using both DESeq2 and edgeR methods.*<br/>
- *__Function 2:__ Plotting volcano maps, based on the results of differentially expressed gene analysis.*<br/>
- *__Function 3:__ Matrix normalization and visualization.Draw heatmap or corrplot, or only calculate TPM, CPM and FPKM.*<br/>
- *__Function 4:__ GO or KEGG enrichment analysis. Enrichment analysis based on the clusterProfiler package.*<br/>
- *__Function 5:__ Plotting bubble maps, based on the results of GO or KEGG enrichment analysis.*<br/>
## Installation
### Install RNAdiff App in Rstudio
Before installing this App, you will need to install some **dependent R packages** on your R.

```{r}
if (!require("BiocManager"))
  install.packages("BiocManager")
library(BiocManager)
if (!require("DESeq2"))
  BiocManager::install("DESeq2")
if (!require("edgeR"))
  BiocManager::install("edgeR")
if (!require("clusterProfiler"))
  BiocManager::install("clusterProfiler")
if (!require("devtools"))
  install.packages("devtools")
```
Once you have completed the installation of the dependencies, start downloading and installing the RNAdiffAPP.
```{r}
devtools::install_github("nongxinshengxin/RNAdiffAPP")
```
Once the installation is complete, run `RNAdiffAPP::run_app()` to open the APP page.
### Install RNAdiff App using docker
```
docker pull mouther/rnadiff
docker run --rm -ti --name h2 -v /home/path/to/shinylog/:/var/log/shiny-server/ -p 3838:3838 mouther/rnadiff
```
Access http://127.0.0.1:3838/RNAdiff/ with a browser
## APP Interface
Analysis of differentially expressed genes

![Alt1](/image/img1.png)

Plotting volcano maps

![Alt2](/image/img2.png)

Calculating TPM and plotting heatmap

![Alt1](/image/img3.png)

GO or KEGG enrichment analysis

![Alt1](/image/img4.png)

Plotting bubble maps

![Alt1](/image/img5.png)
## Documentation
The English documentation is available in - <https://github.com/nongxinshengxin/RNAdiffAPP>

The Chinese documentation is available in - 微信公众号**农心生信工作室**

## Contact us
- Email: nongxinshengxin@163.com
- Wechat Official Account：

![Alt1](/image/wx.png)
