![alt text](rebar.png "Refactored BarSeq analysis")

# Pipeline for BarSeq analysis

Pipeline for analysis of Tn-BarSeq data. The pipeline is based on the scripts of [Morgan Price's Feba repository](https://bitbucket.org/berkeleylab/feba/src/master/). The user is referred to the Feba repository or the original Tn-BarSeq publications for further information.

## Usage

The pipeline currently calculates gene fitness values by using a TnSeq knockout genome mapping "poolfile", a "metadata" file, and gzipped FASTQ files organized in "project" directories under `data/`, `intermediate/`, and `results/`.


### Retrieving data from Illumina basespace *via* command line

Data in form of `*.fastq` files can be manually downloaded from the basespace website on MacOS or Windows.
For Linux systems, only the command line option is available via Illumina's basespace client `bs-cp`. Files are in Illumina's proprietary format. Execute the following line in a terminal and replace `<your-run-ID>` with the number you will find in the URL of your browser. For example, log in to basespace, navigate to `runs`, select a sequencing run and copy the ID you find in the URL: `https://basespace.illumina.com/run/200872678/details`.

```
bs-cp -v https://basespace.illumina.com/Run/<your-run-ID> /your/target/directory/
```

The data must then be converted to `*.fastq` (plain text) files using Illumina's `bcl2fastq` tool. If it complains about indices being too similar to demultiplex, the command has to be executed with option `--barcode-mismatches 0`.

```
cd /your/target/directory/
bcl2fastq
```

The gzipped `*.fastq.gz` files will be stored in `data/projects/example/Data/Intensities/BaseCalls/`. To merge several lanes of the same sample into a new `*.fastq.gz` file, run the following script (pointer has to be in the `rebar/` directory). Output files will be saved to `data/projects/example/`.

```
source source/merge_fastq_files.sh example
```

### Input files

Before starting, a poolfile, a metadata file, and gzipped FASTQ files must be stored under `data/projects/`:

```
data/projects/example.poolfile.tab
data/projects/example.metadata.tab
data/projects/example/*.fastq.gz
```

The poolfile describes TnSeq mappings, and is obtained as output from the [TnSeq pipeline](https://github.com/m-jahn/TnSeq-pipe).
It has the following structure (note: the column for gene names is currently hard-coded as `old_locus_tag`).

```
barcode  rcbarcode  nTot     n scaffold strand    pos   begin     end gene_strand desc  old_locus_tag new_locus_tag gene_length pos_relative central
CAGAAGA… CCCCGCCC…     1     1 NC_0083… +      1.02e6 1014705 1015328 -           gene  H16_B0896     H16_RS23215           623        0.475 TRUE   
CTGTTGG… ACCAACCC…     1     1 NC_0083… +      3.12e6 3118049 3119266 +           gene  H16_A2889     H16_RS14395          1217        0.411 TRUE   
AGCCGCG… GTCCCCCT…     1     1 NC_0083… -      3.44e6 3442926 3443798 -           gene  H16_A3183     H16_RS15875           872        0.861 TRUE   
CGTCATG… CCACCGCT…     1     1 NC_0083… +      3.46e6 3464096 3464620 -           gene  H16_A3206     H16_RS34325           524        0.960 FALSE  
...
```

The tab-separated metadata file contains sample descriptions for the FASTQ files. The `Filename`s must match the names in `/data/projects/example/`

```
Filename          ID Condition Replicate Date       Time          Reference
01_cond1gen0_1     1 generic           1 2020-12-05 0 generations TRUE     
02_cond1gen0_2     2 generic           2 2020-12-05 0 generations TRUE     
03_cond1gen8_1     3 generic           1 2020-12-08 8 generations FALSE    
04_cond1gen8_2     4 generic           2 2020-12-08 8 generations FALSE  
```

The gzipped FASTQ files (one per sample) need to be placed in a project-unique directory under `data/`:

```
data/projects/example/samplefile.fastq.gz
```

### Running the analysis

#### Step 1: Extract per-sample barcode counts

Run the wrapper for the `MultiCodes.pl` script from [FEBA](https://bitbucket.org/berkeleylab/feba/src/master/) on your project:

```
source/run_MultiCodes.sh example
```

Creates `.codes` and `.counts` files in `ìntermediate/example/`.

#### Step 2: Combine BarSeq data and genome mappings

Run the `combineBarSeq.pl` script from [FEBA](https://bitbucket.org/berkeleylab/feba/src/master/) on the project:

```
source/run_combineBarSeq.sh example
```

Creates `.colsum` and `.poolcount` files in `results/example/`.

#### Step 3: Calculate gene fitness

Calculate gene fitness for the project samples using the method described in [Wetmore 2015](https://mbio.asm.org/content/6/3/e00306-15.full):
Note that the column for mapped gene names is currently hard-coded. This can be customized in the header of the script.

```
source/calculate_gene_fitness.R example
```

Creates the fitness tables `results/projects/example/example.fitness.tab.gz` for all strains (barcodes), including data per gene (columns `Strains_per_gene`, `Norm_fg`, `t`, `Significant`), and `results/projects/example/example.gene_fitness.tab.gz` for gene fitness data only (columns `Counts` and `n0` are summed over all strains [barcodes] for each gene, column `log2FC` is log2(`Counts`/`n0`)).

The columns have the following contents:

| Column | Description |
| ------ | ----------- |
| barcode | Barcode sequence |
| locusId | Locus ID (gene name) |
| scaffoldId | Name of DNA molecule |
| Date | Sample date batch; Variable connecting samples to t0 samples |
| Sample | Sample name |
| Condition | Growth condition |
| Counts | Read count for strain (barcode) in sample (summed per gene in gene_fitness table) |
| n0 | Read count in corresponding t0 samples |
| Strains_per_gene | Number of strains (barcodes) for the current locusId |
| Strain_fitness | Strain (barcode) fitness; f_s on p.12 in Wetmore 2015 |
| Norm_fg | Normalized gene fitness; (iii) on p.13 in Wetmore 2015 |
| t | t-like test statistic; calculated on p.13 in Wetmore 2015 |
| Significant | Significant gene if \|t\| > 4; stated on p.3 in Wetmore 2015 |
| log2FC | log2(Counts/n0), only in gene_fitness table |

#### Step 4: Check fitness results with PCA

Generate a PCA plot based on the gene fitness values:

```
source/PCA.R example
```

Creates the PCA plot `results/projects/example/example.PCA.pdf`.

## Authors

Johannes Asplund-Samuelsson, KTH (johannes.asplund.samuelsson@scilifelab.se)

Qi Chen, KTH (qiche@kth.se)

Michael Jahn, KTH (michael.jahn@scilifelab.se)
