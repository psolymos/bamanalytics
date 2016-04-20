# Settings

1. number of nodes to use (12 cores per node), use 5 or 10 as set in the PSB file
2. species acronym (e.g. `"CAWA"`)

## Loading git on jasper

`module load application/git/1.7.10.1`

## Loading R on jasper

`module load application/R/3.1.2`

## script to update the runs

```
rm ~/bam/*
module load application/git/1.7.10.1
cd ~/repos/bamanalytics/
git pull
cd ~/repos/bragging/
git pull
cd ~/bam/
cp ~/repos/bamanalytics/wg/* ~/bam/
```
