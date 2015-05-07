##---
##title: "QPAD estimation"
##author: "Peter Solymos"
##date: "May 7, 2015"
##output: 
##  pdf_document: 
##    toc: true
##    toc_depth: 2
##---

### Preliminaries

## Define root folder where data are stored
ROOT <- "c:/bam/May2015"

## Load required packages
library(mefa4)
library(pbapply)
library(detect)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

## Load preprocesses data
load(file.path(ROOT, "out", "new_offset_data_package_2015-05-07.Rdata"))


