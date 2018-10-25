# synth2
Workflow for creating artificial GWAS data sets based on real data


# gensimul 

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

## Other scripts ##

* `refam.py`: Python program that takes a FAM file, anonymises, removes sex and phenotype data

   The options are:
   * `--desex`  Puts 0 in each entry in the sex column (else leaves alone)
   * `--depheno` Puts a 0 in each entry in the pheno colum (else leaves alone)
   * `--batch`  If used creates a pheno file with fake batches. Give as arguments comma separated lists of ranges (1-indexed). The individuals in each range are given a distinct batch number. e.g. `--batch 30-50,9000-10000` would put all individuals from line 30 in the fam file to and including line 49 in batch 1, and individuals from 9000 up to and including 9999 in batch 2. All others would be in batch 0.
```python3 /projects/scott/shaze/synth2/scripts/refam.py --batch 4149-6548 ```

   * `--miss-batch x.imiss num_under num_above cut`  This is used to batch according to missingness. Using missingness specified by the imiss file, choose num_under below the cut-off given and num_above above the cut-off given. 

# Sample run

```
PDIR=...
DATA=..
MUTES=..
RESDIR=
$PDIR/gensimul/gensimul $DATA   2018 ${DATA}.clt $MUTES /tmp/step0
cp ${DATA}.bim /tmp/step0.bim
cd /tmp
plink --bfile step0 --impute-sex --make-bed --out step1
plink --bfile step0  --missing --out step1
python3 ${PDIR}/scripts/refam.py  --miss-batch step1.imiss 640 60 0.014 KT step1.fam  step2.fam
cp step1.{bed,bim} $RESDIR/
cp step2.fam $RESDIR/step1.fam
cp step2.phe $RESDIR/step1.phe
```
