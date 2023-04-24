# CytoAnalyze

## Description

The CytoAnalyze package is an analysis pipeline for single-cell cytometry data. The package includes high-dimensional analysis, statistical discovery of differential cell abundances and bio-marker expressions, and a variety of visualization tools. This pipeline is based on open source R packages and is flexible to accommodate different study designs. It aims primarily to bring in one pipeline the tools to explore single-cell cytometry data to facilitate the work of biology researchers, and increase time-efficiency and reproducibility of single-cell cytometry data mining.

## Summary

The package includes several functions to perform different analysis tasks, and one wrapper function to streamline all the required analysis steps. For the workflow demonstration, please see [the workflow vignette](doc).

## Installation

``` {r }
if(!require(devtools)) {
  install.packages("devtools")
}

devtools::install_github("RamyGadalla/CytoAnalyze")
```

## Usage

CytoAnalyze input is cleaned-up single cell cytometry data as .FCS files or SummarizedExperiment-class object. CytoAnalyze does not include data cleaning feature. The output of CytoAnalyze are a variety of visualization plots and tables with all the necessary values to explore and provide insight about the data set.
