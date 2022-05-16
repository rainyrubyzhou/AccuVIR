# AccuVIR: Accurate viral enome long reads polisher
=======================================================================

AccuVIR--a **Acc**urate **VIR**al genome polisher-- utilises long reads within a single run to polish a long reads assembly of small and large genomes. 

AccuVIR requires the following as input:
+ Short reads/HiFi reads (in FASTA/FASTQ format; can be compressed)
+ Draft contigs (in FASTA/FASTQ format; can be compressed)
 
### E-mail: runzhouyu2-c@my.cityu.edu.hk

## Installation
AccuVIR is only available for Unix-like platforms (Linux and MAC OS). We recommend using *Option_1* for a convenient installation and in the case where target machine is different from the one on which AccuVIR is installed. On the other hand, *Option_2* is more suitable if the binary of AccuVIR is to be run on the same machine on which it is compiled because then a machine-specific optimised (and thus slightly faster) binary can be produced using the flag `-Doptimise_for_native=ON`.


### Option_1: Conda Package Installation
The convenient way of installation is using the conda package as follows:
```console
conda install -c bioconda AccuVIR
```
Htslib 1.10 may cause conflicts with some of the already installed packages. In such a case, AccuVIR can be installed in a new environment as follows:
```console
conda create --name AccuVIR_env
conda activate AccuVIR_env  
conda install -c bioconda AccuVIR
```

### Option_2: Installation from the source
CmakeLists is provided in the project root folder. 

#### Dependencies

- bioconda
- Biopython>=1.70
- Required python package: Biopython>=1.70, networkx >= 2.5.1, pandas >= 1.1.3
- [HTSLIB](https://github.com/samtools/htslib) (version >=1.10)
  + If htslib version 1.10 or higher is not installed, we recommend using `install_deps.sh` in the project folder to install it locally.

- [KMC3](https://github.com/refresh-bio/KMC)
  + Required for suk
  + Can also be installed using conda as follows: 
  ```console
  conda install -c bioconda kmc
  ``` 

#### Building Executable
Run the following commands to build a binary (executable) `AccuVIR` in `build/bin` :
If htslib version 1.10 or higher is installed:
```console
  git clone --recursive https://github.com/rainyrubyzhou/AccuVIR AccuVIR
  cd AccuVIR
```
If htslib is not installed or the version is smaller than 1.10:
```console
  git clone --recursive https://github.com/kensung-lab/AccuVIR AccuVIR
  cd AccuVIR
  chmod +x install_deps.sh
  ./install_deps.sh
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release -Doptimise_for_native=ON ..
  make -j 8
```


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

 	-h, --help
	Print the usage. 
```

### Output Results
If no `--output` (or `-o`) is provided, `AccuVIR_X.fasta` will be created in the working directory where <X> is the name of the draft file. 

### Resource Usage
The option `--processing-size` (or `-p`) controls the number of contigs processed in one batch. By default, all contigs will be processed in a single batch. More the number of contigs processed in a batch higher will be the memory used. Lower number, on the other hand, may not utilise the number of threads efficiently. As a reference, for the whole human genome we used `-p 96` on a 48 core machine which used about 380G RAM and finished its run in about 3 hours (only Illumina polishing). We recommend using `-p` to be a multiple of `-t`. If the genome size is not too large for the available RAM, we recommend processing all the contigs in a single batch (i.e. avoid specifying `-p`).

### Example 1 (using Illumina paired-end reads)

#### Mapping the short reads to contigs:
Assuming `$R1`, `$R2`, `$Draft` contain the names of the files containing short reads (paired-end) and draft contigs, respectively. Let `$NUMTH` represents the number of threads  to be used.
```console
minimap2 --secondary=no --MD -ax sr -t $NUMTH $DRAFT $R1 $R2 | samtools view -Sb - > mapped-sr.bam
samtools sort -@$NUMTH -o mapped-sr.sorted.bam mapped-sr.bam
samtools index mapped-sr.sorted.bam
rm mapped-sr.bam
```

#### Mapping the long reads to contigs:
Assuming `$LONGR` and `$Draft` contain the names of the files containing the long reads (PacBio or ONT) and draft contigs, respectively. Let `$NUMTH` represents the number of threads to be used. `$RTYPE` will be `map-pb` for PacBio and `map-ont` for ONT reads.
```console
minimap2 --secondary=no --MD -ax $RTYPE -t $NUMTH $DRAFT $LONGR | samtools view -Sb - > mapped-lg.bam
samtools sort -@$NUMTH -o mapped-lg.sorted.bam mapped-lg.bam
samtools index mapped-lg.sorted.bam
rm mapped-lg.bam
```

#### Running AccuVIR:
Create a text file containing the names of the short reads files.
```console
echo -e "$R1\n$R2" > il_names.txt
```

Let genome size be around 3Gbp and average coverage of Illumina reads is around 55. Say, we would like to use 48 threads and process 96 contigs in a single batch.

A sample run of AccuVIR (for only short reads polishing) can be:
```console
./bin/AccuVIR -d $DRAFT -r @il_names.txt -s 3g -c 55 -b mapped-sr.sorted.bam -p 96 -t 48 -o whole_genome.h.fa
```

A sample run of AccuVIR (for short reads as well as long reads polishing) can be:
```console
./bin/AccuVIR -d $DRAFT -r @il_names.txt -s 3g -c 55 -b mapped-sr.sorted.bam -B mapped-lg.sorted.bam -p 96 -t 48 -o whole_genome.h2.fa
```
### Example 2 (using CCS reads)

#### Mapping the CCS reads to contigs:
Assuming `$READS` and `$Draft` contain the names of the files containing CCS reads and draft contigs, respectively. Let `$NUMTH` represents the number of threads  to be used.
```console
minimap2 --secondary=no --MD -ax asm20 -t $NUMTH $DRAFT $READS | samtools view -Sb - > mapped-ccs.bam
samtools sort -@$NUMTH -o mapped-ccs.sorted.bam mapped-ccs.bam
samtools index mapped-ccs.sorted.bam
rm mapped-ccs.bam
```

#### Running AccuVIR:
Let genome size be around 3Gbp and average coverage of CCS reads is around 30. Say, we would like to use 48 threads and process all the contigs in a single batch.
A sample run of AccuVIR (for only CCS polishing) can be:
```console
./bin/AccuVIR -d $DRAFT -r $READS -s 3g -c 30 -b mapped-ccs.sorted.bam -t 48 -o whole_genome.h.fa
```

## Method and Results
For the whole human genome (HG002) assembly, AccuVIR took about 3 hours and about 380G RAM to polish (on a 48 cores machine with 48 threads) using only Illumina reads. For polishing using Illumina as well as PacBio reads, time taken was about 4 hours and 15 minutes using about 410G RAM. The method and partial results can be found at [BioRxiv](https://www.biorxiv.org/content/10.1101/2019.12.19.882506v1).

## External Libraries

 * [sdsl-lite](https://github.com/simongog/sdsl-lite) has been used for rank-select and bit-vectors data-structures.
 * [slog](https://github.com/Ritu-Kundu/slog) has been used to print time and memory usage.
 * [suk](https://github.com/Ritu-Kundu/suk) has been used as the module to compute the solid (unique) kmers.
 * An adapted version of [spoa](https://github.com/rvaser/spoa.git) library has been used for POA.

## Contact
Other than raising issues on Github, you can contact Ritu Kundu (dcsritu@nus.edu.sg) or Joshua Casey (joshuac@comp.nus.edu.sg) for getting help in installation/usage or any other related query.




