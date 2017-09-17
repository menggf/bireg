---
title: 'bireg: A tool to predict the transciptional regulation associated with part of samples'
author: "Guofeng Meng"
date: '2017-09-02'
output:
  pdf_document:
    keep_tex: yes
    number_sections: yes
    toc: yes
    toc_depth: 2
  html_document:
    toc: yes
    toc_depth: '2'
vignette: "\\VignetteIndexEntry{bireg: A tool to predict the
transciptional regulation associated with part of samples} \n\\VignetteEncoding{UTF-8} \n\\VignetteEngine{knitr::rmarkdown}\n"
---

<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{bireg: A tool to predict the transciptional regulation associated with part of samples}
-->

## Introduction

bireg: A tool to predict the transciptional regulation associated with part of samples. This package impletement a bi-clustering analysis to predicted the target genes of input regulators. In this work, regulation is measured by Spearman's correlation and uses a strict correlation threshold to declare the regulation relationship. Different from the traditional co-expression network, this package can discover the regulated only observed in part of patients. 


To use bireg:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("bireg")
```

## Usage
