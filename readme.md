# AccuVIR: Accurate viral genome polisher using long reads
=======================================================================

**AccuVIR**--a **Acc**urate **VIR**al genome polisher-- utilizes path searching and sampling in sequence alignment graphs to polish long reads assembly of viral genomes. 

AccuVIR requires the following as input:
+ Reads file for graph construction (Optimized for third generateion sequencing data).
+ Backbone sequence file for graph construction (Can be either draft assembly from aseembly tools including Canu/Shasta/Flye  etc, or an available genome of another viral subtype/haplotype).
 
## Installation

### Option_1: Conda Package Installation
The convenient way of installation is using the conda package as follows:
```console
conda install -c bioconda AccuVIR
conda create --name AccuVIR
conda activate AccuVIR 
conda install -c bioconda AccuVIR
```

### Option_2: Installation from the source
CmakeLists is provided in the project root folder. 

#### Dependencies
- Conda
- python
- Required python package: Biopython>=1.70, networkx >= 2.5.1, pandas >= 1.1.3

#### Installation
```console
git clone --recursive https://github.com/rainyrubyzhou/AccuVIR AccuVIR
cd AccuVIR/src
python AccuVIR_main.py -h
```
Successfull installation will end with usage information using above commands.


## Usage of the AccuVIR: 

>**Usage:**
```console
python AccuVIR_main.py <args>
```

>**Mandatory args:**
```console
-r | <str, e.g. "raw_HIV.fasta">
Reads file for graph construction (in fasta format). 

-b | <str, e.g. "canu_contig.fasta">
Backbone sequence file for graph construction (in fasta format). 
```
>**Optional args:**
```console
--beamwidth | <int, e.g. 100>
Beamwidth for diverse beam search (default: 500).

-h | Print the usage information. 
```

>**Output Results:** 

Outputs will be created in the input file directory.
+ `X_ON_Y_filtered.fa` is the intermediate output before MRR. It containing outputs from both DBS and path sampling modules.  

`X` is the name of the reads file and `Y` is the backbone sequence. 

## Contact
Other than raising issues on Github, you can also contact YU Runzhou (runzhouyu2-c@my.cityu.edu.hk) for help in installation/usage or any other related query.




