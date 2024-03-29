# AccuVIR: Accurate viral genome assembler and polisher using long reads
=======================================================================

**AccuVIR**--an **Acc**urate **VIR**al genome assembler and polisher -- utilizes path searching and sampling in sequence alignment graphs to assemble or polish draft assembly of viral genomes. Users are welcome to read our published paper at Bioinformatics [here](https://doi.org/10.1093/bioinformatics/btac827) and try our tool as guided below.

AccuVIR requires the following as input:
+ Reads file for graph construction (Optimized for third-generation sequencing data).
+ Backbone sequence file for graph construction (Can be either draft assembly from assembly tools including Canu/Shasta/Flye etc, or an available genome of another viral subtype/haplotype).
 
## Installation

#### Dependencies
- Conda
- python
- [blasr](https://anaconda.org/bioconda/blasr)
- Required python package: Biopython>=1.70, networkx >= 2.5.1, pandas >= 1.1.3, seaborn >= 0.11.1

#### Installation
```console
git clone --recursive https://github.com/rainyrubyzhou/AccuVIR AccuVIR
cd AccuVIR/src
python AccuVIR_main.py -h
```
Successful installation will end with usage information using the above commands.


## Usage of the AccuVIR: 
### Step 1: Paths generation using two modules.

>**Command Usage:**
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
+ `X_ON_Y_filtered.fa` is the intermediate output before MRR. It contains multiple sequences from both DBS and path sampling modules.  

    `X` is the name of the reads file and `Y` is the backbone sequence. Outputs will be created in the input file directory. 

### Step 2: Apply gene prediction tool (Genemark recommended)

Due to the license requirement of [Genemark](http://exon.gatech.edu/GeneMark/) tools, users need to preinstall the tool or run it [online](http://exon.gatech.edu/GeneMark/gmhmmp.cgi) in this step.   
For online running, please 
1. upload the file `X_ON_Y_filtered.fa` as input
2. tick `GFF` as output format 
3. save the output as `X_ON_Y_filtered.fa.gtf` for next step

For offline running, we use the version `GeneMark.hmm for prokaryotes`. Users need to 
1. download `GeneMarkS` at [this page](http://exon.gatech.edu/GeneMark/license_download.cgi)
2. download `gm_key` and put it at required directory
3. run Genemark using the command below
>**Example usage of Genemark (GeneMark.hmm for prokaryotes in this example):**
```console
gmhmmp -m heu_11.mod -f G -o X_ON_Y_filtered.fa.gtf X_ON_Y_filtered.fa
```

`.gtf` output file is required for next step. (e.g. `X_ON_Y_filtered.fa.gtf`)


### Step 3: Call ranking module for final output.
Pass in the `X_ON_Y_filtered.fa` and the sequence of greatest MRR value will be output.
>**Command Usage:**
```console
python AccuVIR_MRR.py -r X_ON_Y_filtered.fa
```
>**Output Results:** 

 + `**_final.fa` is final output of AccuVIR. It contains the single sequence that ranks best using MRR. 

    Users can also pass this sequence as the backbone to step 1 to iteratively refine the output. 

## Data availability
Users can test our tool using the test data [here](https://drive.google.com/drive/folders/1iCNVjkw_LEhd8pYfS4QDXEAmVAHZW2N9). `ref.fa` is the ground truth for this dataset.
```console
python AccuVIR_main.py -r test_reads.fa -b backbone.fa
```


Simulated datasets used in our experimentd are available [here](https://drive.google.com/drive/folders/1jIIBaANO5Gi0EeECuxq_7IHYScds4dDB).
## Contact
Other than raising issues on Github, you can also contact YU Runzhou (runzhouyu2-c@my.cityu.edu.hk) for help in installation/usage or any other related query.




