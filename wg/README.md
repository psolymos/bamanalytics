# Settings

1. number of nodes to use (12 cores per node),
2. `fid` that is file ID to load the data: 1, 2, 3, 4
3. species acronym (e.g. CAWA)

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
