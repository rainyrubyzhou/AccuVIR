# AccuVIR: Accurate viral enome long reads polisher
=======================================================================

**AccuVIR**--a **Acc**urate **VIR**al genome polisher-- utilises long reads within a single run to polish a long reads assembly of small and large genomes. 

AccuVIR requires the following as input:
+ Reads file for graph construction (Optimized for third generateion sequencing data).
+ Backbone sequence file for graph construction.
 
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
  cd AccuVIR
  python python AccuVIR_main.py -v 
```
Successfull installation will  end with version number using above commands.


## Usage of the AccuVIR: 
```console
 Usage: python AccuVIR_main.py <args>

 ** Mandatory args: 
	-r | <str, e.g. "raw_HIV.fasta">
	Reads file for graph construction(in fasta format). 

	-b | <str, e.g. "canu_contig.fasta">
	Backbone sequence file for graph construction (in fasta format). 

** Optional args:
	--beamwidth | <int, e.g. 100>
	Beamwidth for diverse beam search (default: 500).

	-o | <str>
	Output file name. 
	[Default] AccuVIR_<draft_file_name>.fasta in the working directory.

	-v | Print the version number.
 	-h | Print the usage. 
```

### Output Results
If no `--output` (or `-o`) is provided, `AccuVIR_X.fasta` will be created in the working directory where <X> is the name of the draft file. 

## Contact
Other than raising issues on Github, you can also contact YU Runzhou (runzhouyu2-c@my.cityu.edu.hk) for help in installation/usage or any other related query.




