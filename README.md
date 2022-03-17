# Simulate Outbreak Dataset
This pipeline is designed to perform a **rudimentary** simulation of an infectious disease outbreak, and to generate simulated sequence datasets
that resemble the data that would be obtained by isolating and sequencing the microbes responsible for the outbreak. It is intended to be used
to produce datasets that can be used as part of a validation experiment for other microbial genomics pipelines, in particular those that involve
variant calling and core SNP phylogenetics. As such, no great effort is currently made to make the disease outbreak process particularly realistic.
The emphasis of the design is to quickly and easily create relatively large datasets where mutations accumulate over generations along a set of lineages.

As the mutated genomes are created, we also simulate sequencing reads from those genomes. We collect records of all mutations that are introduced, so
we have a 'ground truth' of mutations, variants and phylogenetic relationships that are often difficult (or impossible) to obtain for real datasets. 

Genomic mutations are simulated using [`simuG`](https://github.com/yjx1217/simuG), and illumina paired-end reads are simulated using
[`ART`](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm).

**Note**: this pipeline uses Nextflow features that are currently experimental to support recursive execution of processes. In order to run this
pipeline, you must be using a recent [edge release](https://www.nextflow.io/docs/latest/getstarted.html#stable-edge-releases) of nextflow.

## Dependencies
Dependencies are listed in [environments.yml](environments/environment.yml). When using the `conda` profile, they will be installed automatically.
Other runtime environments such as docker, podman or singularity are not currently supported.

## Usage
The pipeline takes a single reference genome as input. That genome will be treated as the infectious organism that causes the index case of infection.

To run the pipeline with default parameters:

```
nextflow run dfornika/simulate-outbreak-dataset \
  -profile conda \
  --cache ~/.conda/envs \
  --ref </path/to/reference.fa> \
  --outdir </path/to/outdir>
```

Several parameters are available to control the outbreak simulation.

| Parameter                | Description                                                        | Default |
|:-------------------------|:-------------------------------------------------------------------|--------:|
| `iterations`             | Number of times to iterate through the mutation-selection process. | 8       |
| `min_snps_per_iteration` | Minimum number of SNPs to introduce during mutation phase.         | 2       |
| `max_snps_per_iteration` | Maximum number of SNPs to introduce during mutation phase.         | 20      |

## Outputs
An output directory will be created using the filename of the initial reference genome (excluding the file extension).
For each derived genome, a separate directory will be created, using a 4-character randomly generated identifier.

```
output
├── 381D
│   ├── 381D.fa
│   ├── 381D_mutation_info.csv
│   ├── 381D_R1.fastq.gz
│   ├── 381D_R2.fastq.gz
│   ├── 381D.variants.tsv
│   └── 381D.vcf
├── 4798
│   ├── 4798.fa
│   ├── 4798_mutation_info.csv
│   ├── 4798_R1.fastq.gz
│   ├── 4798_R2.fastq.gz
│   ├── 4798.variants.tsv
│   └── 4798.vcf
├── 5A68
│   ├── 5A68.fa
│   ├── 5A68_mutation_info.csv
│   ├── 5A68_R1.fastq.gz
│   ├── 5A68_R2.fastq.gz
│   ├── 5A68.variants.tsv
│   └── 5A68.vcf
├── F272
│   ├── F272.fa
│   ├── F272_mutation_info.csv
│   ├── F272_R1.fastq.gz
│   ├── F272_R2.fastq.gz
│   ├── F272.variants.tsv
│   └── F272.vcf
└── NC_000962.3
    ├── NC_000962.3_R1.fastq.gz
    └── NC_000962.3_R2.fastq.gz
```

The `_mutation_info.csv` file includes two headers that make it possible to re-construct the simulated chain of transmission.

```
parent_id
num_snps
```

For each simulated isolate, a FASTA-formatted 'true' genome is produced, along with simulated illumina paired-end sequence data.

A `.variants.tsv` file is produced with headers:

```
ref_chr
ref_start
ref_end
ref_strand
ref_allele
sim_chr
sim_start
sim_end
ref_strand
sim_allele
variant_type
variant_id
donor_chr_in_ref
donor_start_in_ref
donor_end_in_ref
donor_strand_in_ref
duplication_type
inserted_copy_number
total_copy_number
```

A `.vcf` formatted file of the simulated mutations is also provided.
