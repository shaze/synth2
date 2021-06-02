# synth2
Workflow for creating artificial GWAS data sets based on real data


Several separate projects

#gensimul#

This was our first version -- good thing is that it keeps structure very well but the resulting data is relatively identifiable and
also we land up with smaller data sets than our original

This is C program that takes a source plink file, anonymises
  and does some randomisation. Here is an example. 

```
   ./gensimul src   200 src.clt mutes  new  
```
   Here
   * `src` is the base name for the source PLINK files
   * 200 is the number of people you want in the new data set
   * `src.clt` is a cluster file (see below)
   * `mutes` is a list of mutations (see below)
   * `new`: the base name of the new PLINK files

## Cluster file ##

The cluster file contains a list of clusters. Each cluster is on a line by itself with space separated integers.  e.g.

```
1 15 39
2 7
3
4 8 124
20 11 5
```

There should be no overlap in clusters. In the current version, there should be a maximum of 3 members of each cluster. We hope/plan/sort of promise to generalise. The numbers are the 0-indexed entries in the FAM file.

The idea is that each cluster consists of closely related individuals. In our new data set we have a randomised amalgam from thse individuals. See below for how to create this.

## Mutation file

The file consists of a list of mutations, each mutation on a line by itself, withree entries per line:
* The 0-indexed entry in the BIM file of the mutation
* The desired MAF of cases (between 0 and 1)
* The desired MAF of controls

###  Creating a cluster file

The cluster file can be created in two steps.

Use PLINK to make a PLINK-style cluster file

```
plink --bfile data --cluster --ppc 0.01  --mc 3  --out data
```

Then call our auxiliary script to produce in the right format

```
python3 munge_cluster.py data data.clt
```



# hapsimul

This is version 2 -- takes a largish input data set and produces a larger output randomised data set that preserves some LD but collapses
population structure


# simulstruct

This is a Nextflow script that uses hapsimul to get some population structure back in the data
