![alt text](rebar.png "Ralstonia eutropha BarSeq analysis")

# _Ralstonia eutropha_ BarSeq analysis

Scripts for analysis of _Ralstonia eutropha_ Tn-BarSeq data.

## Usage

The pipeline currently calculates gene fitness values by using a TnSeq knockout genome mapping "poolfile", a "metadata" file, and gzipped FASTQ files organized in "project" directories under `data/`, `intermediate/`, and `results/`.

### Input data

Before starting, a poolfile, a metadata file, and gzipped FASTQ files must be stored under `data/projects/`:

```
data/projects/rebar.poolfile.tab
data/projects/rebar.metadata.tab
data/projects/rebar/*.fastq.gz
```

The poolfile describes TnSeq mappings and looks like this:

```
barcode	rcbarcode	nTot	n	scaffold	strand	pos	n2	scaffold2	strand2	pos2	nPastEnd
TCTTGGTCCAGGAAGCCGAT	ATCGGCTTCCTGGACCAAGA	83	80	AY305378	+	325252	1	NC_008313	-	2508065	0
TGGGTGATCTTGGGTGTAGG	CCTACACCCAAGATCACCCA	32	31	NC_008313	-	100293	1	NC_008314	+	1850688	0
CCGGGTAAGCGTTTGTGTTG	CAACACAAACGCTTACCCGG	10	9	NC_008314	-	649441	1	NC_008313	-	342601	0
GCTATGTATACATAGTGACT	AGTCACTATGTATACATAGC	956	925	NC_008313	+	701767	6	NC_008313	+	701768	0
...
```

The metadata file contains sample descriptions for the FASTQ files:

```
ID	Date	Condition	Sample
IT001	8/14	Fructose	1a
IT002	8/14	Succinate	2a
IT003	8/28	Formate	3a
IT004	8/14	Fructose	1b
...
```

The gzipped FASTQ files (one per sample) need to be placed in a project-unique directory under `data/`:
```
data/projects/rebar/
```

...and they should be named as such with IDs matching the metadata file:

```
Fructose_IT001_1.fastq.gz
Succinate_IT002_1.fastq.gz
Formate_IT003_1.fastq.gz
Fructose_IT004_1.fastq.gz
...
```

### Running the analysis

Note that the paths to the Perl scripts are hardcoded and need to be changed for each installation.

#### Step 1: Extract per-sample barcode counts

Run the `MultiCodes.pl` script from FEBA on the project:

```
source/run_MultiCodes.sh rebar
```

Creates "codes" and "counts" files in `ìntermediate/rebar/`.

#### Step 2: Combine BarSeq data and genome mappings

Run the `combineBarSeq.pl` script on the project:

```
source/run_combineBarSeq.sh rebar
```

Creates "colsum" and "poolcount" files in `results/rebar/`.

#### Step 3: Calculate gene fitness

Calculate gene fitness for the project samples:

```
source/calculate_gene_fitness.R rebar
```

Creates the fitness table `results/projects/rebar/rebar.fitness.tab.gz`.

#### Step 4: Check fitness results with PCA

Generate a PCA plot based on the gene fitness values:

```
source/PCA.R rebar
```

Creates the PCA plot `results/projects/rebar/rebar.PCA.pdf"`.

## Authors

Johannes Asplund-Samuelsson, KTH (johannes.asplund.samuelsson@scilifelab.se)

Qi Chen, KTH (qiche@kth.se)
