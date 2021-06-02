
#simulstruct

A workflow to create simulated GWAS data that has some population structure in it

NB: THis is dependant on the `hapsimul` code  being made


## Parameters

The parameters and their default values are shown below

* `source`  : no default -- name of input file
* `out_dir="output"  : which directory output should go
* `out` ("sim)       : name of the PLINK output files
* `expansion` (2)       : factory by which the output is bigger than the inputs
* `K` (10)           : Number of groups we should look for in input
