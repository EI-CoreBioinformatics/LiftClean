# LiftClean
Gemy George Kaithakottil, David Swarbreck  
[![DOI](https://zenodo.org/badge/1063856495.svg)](https://doi.org/10.5281/zenodo.18166872)

## Abstract

LiftClean detects potential problems in lifted-over (projected) annotations that may be incompatible with coding gene models. These issues include internal stop codons, as well as structural anomalies such as book-ended exons (zero-length introns), overlapping CDS segments, or overlapping exon features.

LiftClean uses [Mikado](https://github.com/EI-CoreBioinformatics/mikado?tab=readme-ov-file#mikado---pick-your-transcript-a-pipeline-to-determine-and-select-the-best-rna-seq-prediction) Prepare to identify and categorise these issues, including producing summary plots. It then filters the problematic annotations (with user-configurable options) and generates a cleaned GFF file suitable for downstream genome annotation workflows.

The LiftClean utility can therefore refine [LiftOn](https://github.com/Kuanhao-Chao/LiftOn) and [Liftoff](https://github.com/agshumate/Liftoff) outputs by producing:
- a cleaned, ready-to-use GFF file,
- a detailed report, and
- PNG plots summarising the warnings detected during processing.

## Installation

All installation methods below will install LiftClean along with its dependencies.

### Docker Installation
LiftClean can be installed with Docker. If you don't have Docker, please install [Docker](https://docs.docker.com/get-docker/) first. Then you can pull the Docker image with LiftClean installed

```console
VERSION=0.1.0
docker run gemygk/liftclean:v${VERSION} liftclean -h
```

### Singularity Installation
LiftClean can be installed with Singularity. If you don't have Singularity, please install [Singularity](https://docs.sylabs.io/guides/3.9/user-guide/quick_start.html#quick-installation-steps) first. Then you can pull the singularity image with LiftClean installed.

We can directly run LiftClean from the Singularity image hosted on DockerHub
```console
VERSION=0.1.0
singularity exec docker://gemygk/liftclean:v${VERSION} liftclean -h
```

Or, we can build and run a Singularity image, following the steps below:
```console
# Create a Singularity definition file, like below:

$ cat liftclean-0.1.0.def
bootstrap: docker
from: gemygk/liftclean:v0.1.0

# Build the Singularity image
$ sudo singularity build liftclean-0.1.0.sif liftclean-0.1.0.def

# Execute LiftClean from the Singularity image
$ singularity exec liftclean-0.1.0.sif liftclean -h
```

### Manual Installation

#### Install dependencies

* Mikado - [Installation](https://github.com/EI-CoreBioinformatics/mikado?tab=readme-ov-file#installation)  
* GffRead - [Installation](https://github.com/gpertea/gffread?tab=readme-ov-file#installation)  
* GenomeTools - [Installation](https://github.com/genometools/genometools?tab=readme-ov-file#building-and-installation)  

#### Get LiftClean
First, obtain the source code:
```console
git clone https://github.com/EI-CoreBioinformatics/LiftClean.git
cd LiftClean
```
Build and install using [UV](https://github.com/astral-sh/uv?tab=readme-ov-file#uv)
```console
version=0.1.0 \
     && uv build \
     && pip install --prefix=/path/to/software/liftclean/${version}/x86_64 -U dist/*whl
```

Also, make sure that both `PATH` and `PYTHONPATH` (below is for python3.10) environments are updated
```console
export PATH=/path/to/software/liftclean/${version}/x86_64/bin:$PATH
export PYTHONPATH=/path/to/software/liftclean/${version}/x86_64/lib/python3.10/site-packages
```

### Usage
```console
$ liftclean --help
usage: liftclean [-h] -g GENOME_FASTA [-n LIFTON_GFF] [-f LIFTOFF_GFF] [--alt_lifton_label ALT_LIFTON_LABEL] [--alt_liftoff_label ALT_LIFTOFF_LABEL] [-p PREFIX] [-o OUTPUT] [-s] [-t THREADS] [-e EXCLUDE_FROM_FILTERING] [-m MINIMUM_CDNA_LENGTH]
                 [-i MIN_INTRON_LENGTH] [--limit_filters] [--check_dup_ids] [--transcript_types TRANSCRIPT_TYPES] [--gene_types GENE_TYPES] [--gffread_params GFFREAD_PARAMS] [--gt_gff3_params GT_GFF3_PARAMS]
                 [--mikado_prepare_params MIKADO_PREPARE_PARAMS] [--force] [-d]

        Lifton/Liftoff transcript filtering and comparison pipeline


options:
  -h, --help            show this help message and exit
  -g GENOME_FASTA, --genome_fasta GENOME_FASTA
                        Provide reference genome FASTA file (default: None)
  -n LIFTON_GFF, --lifton_gff LIFTON_GFF
                        Provide Lifton output GFF file (default: None)
  -f LIFTOFF_GFF, --liftoff_gff LIFTOFF_GFF
                        Provide Liftoff output GFF file (default: None)
  --alt_lifton_label ALT_LIFTON_LABEL
                        Alternative label for Lifton in plots and outputs [default:lifton]
  --alt_liftoff_label ALT_LIFTOFF_LABEL
                        Alternative label for Liftoff in plots and outputs [default:liftoff]
  -p PREFIX, --prefix PREFIX
                        Provide a label for file name prefix and plot titles to distinguish images from multiple runs. If provided, it will be suffixed to the --alt_lifton_label and --alt_liftoff_label. For example, you can use 'Wheat_Accession1' and the output will be labeled accordingly like lifton_Wheat_Accession1 or liftoff_Wheat_Accession1 [default:None]
  -o OUTPUT, --output OUTPUT
                        Provide output directory [default:output]
  -s, --single          Only one GFF allowed. Default is to process both --lifton_gff and --liftoff_gff inputs [default:False]
  -t THREADS, --threads THREADS
                        Number of threads to use for copying (default: 1)
  -e EXCLUDE_FROM_FILTERING, --exclude_from_filtering EXCLUDE_FROM_FILTERING
                        Provide a comma-separated list of types to be excluded from being filtered, for example, 'Size under minimum,Redundant' etc [default:]
  -m MINIMUM_CDNA_LENGTH, --minimum_cdna_length MINIMUM_CDNA_LENGTH
                        Provide minimum cDNA length for filtering. Anything lower will be excluded [default: 48]
  -i MIN_INTRON_LENGTH, --min_intron_length MIN_INTRON_LENGTH
                        Provide intron size. Any models with intron size under this value will be removed. Bookended and overlapping exon features are merged with the default setting [default: 0]
  --limit_filters       Enable this to not exclude models based on categories: Size under minimum, Incorrect fusions of splice junctions, Cannot reverse strand of coding transcript, Redundant [default:False]
  --check_dup_ids       Check for duplicate IDs, but only for features that are parents; if found, remove all instances and all descendants [default:True]
  --transcript_types TRANSCRIPT_TYPES
                        Comma-separated list of feature types to treat as transcripts (default: mRNA,primary_transcript,transcript,lnc_RNA,ncRNA,miRNA,rRNA,tRNA,snoRNA,snRNA,scaRNA,pseudogenic_transcript,antisense_RNA)
  --gene_types GENE_TYPES
                        Comma-separated list of feature types to treat as gene-level (default: gene,ncRNA_gene,pseudogene)
  --gffread_params GFFREAD_PARAMS
                        Additional parameters to pass to gffread (enclosed in quotes). Only modify defaults if you know what you are doing. MUST use the assignment format --gffread_params="<params>" [default:"--keep-genes -F"]
  --gt_gff3_params GT_GFF3_PARAMS
                        Additional parameters to pass to gt (enclosed in quotes). Only modify defaults if you know what you are doing. MUST use the assignment format --gt_gff3_params="<params>" [default:"-sort -tidy -retainids yes"]
  --mikado_prepare_params MIKADO_PREPARE_PARAMS
                        Additional parameters to pass to mikado (enclosed in quotes). Only modify defaults if you know what you are doing. MUST use the assignment format --mikado_prepare_params="<params>" [default:""]
  --force               Enable this option to force removal of pre-existing output folder [default:False]
  -d, --debug           Enable this option for debugging [default:False]

Note:
  - In --single mode, you must provide exactly one of --lifton_gff or --liftoff_gff.
  - In paired mode, both --lifton_gff and --liftoff_gff are required.
  - The script will symlink input files into the output directory.
  - Filtering is based on Mikado Prepare identified warnings. Categories of warnings can be excluded i.e. not filtered via --exclude_from_filtering. Full list of categories can be found below:
	5'UTR present with a truncated ORF,
	Assertion failure,
	Assertion failure start must be less than end,
	Both UTR present with truncated ORF,
	CDS which straddles 2 different exons,
	Cannot reverse strand of coding transcript,
	Debords its exon,
	Defined UTRs but no CDS feature,
	Duplicate parent feature ID,
	General,
	Incorrect fusions of splice junctions,
	Internal stop codons found,
	Invalid CDS length,
	Invalid number of coding exons,
	Invalid start and stop of the ORF,
	Overlapping CDS,
	Overlapping exons found,
	Redundant,
	Seqid mismatch*,
	Short intron,
	Size under minimum,
	Strand conflict child*,
	Strand conflict gene*,
	Strand conflict gene-child*
    * Categories marked with * cannot be excluded from being filtered.

   - The script assumes the presence of external tools:
      - gffread
      - gt (GenomeTools)
      - mikado
```

## Running LiftClean

### Example Command
To run LiftClean, use the command line interface with the required arguments for the genome FASTA file and GFF3 annotation file (Lifton or Liftoff or both). For example:

#### Using both LiftOn and Liftoff GFFs
```console
$ cd /path/to/work_directory
$ liftclean \
    --genome_fasta input_genome.fasta \
    --lifton_gff input_lifton.gff
    --liftoff_gff input_liftoff.gff

# Output will be saved to 'output' directory by default
output
├── genome.fasta -> /path/to/work_directory/input_genome.fasta
├── liftoff.gff -> /path/to/work_directory/input_liftoff.gff
├── lifton.gff -> /path/to/work_directory/input_lifton.gff
├── genome.fasta.fai
├── output_Rejected_Transcripts_summary_plot.png
├── analysis
│   ├── lifton.corrected.gff
│   ├── lifton.corrected.gff_strand_checker.log
│   ├── lifton.short_introns.corrected.gff
│   ├── lifton.short_introns.mapping.tsv
│   ├── lifton.short_introns.tsv
│   ├── lifton.corrected.gffread.gff
│   ├── lifton.sorted.gff
│   ├── lifton.sorted.gff.log
│   ├── liftoff.corrected.gff
│   ├── liftoff.corrected.gff_strand_checker.log
│   ├── liftoff.short_introns.corrected.gff
│   ├── liftoff.short_introns.mapping.tsv
│   ├── liftoff.short_introns.tsv
│   ├── liftoff.corrected.gffread.gff
│   ├── liftoff.sorted.gff
│   ├── liftoff.sorted.gff.log
│   ├── liftoff.sorted.ids.retained.txt
│   ├── liftoff.sorted.ids.tsv
│   ├── lifton.sorted.ids.retained.txt
│   ├── lifton.sorted.ids.tsv
│   ├── mikado_prepare_liftoff
│   │   ├── list.txt
│   │   ├── mikado_prepare_liftoff.fasta
│   │   ├── mikado_prepare_liftoff.gtf
│   │   ├── mikado_prepare_liftoff.log
│   │   ├── mikado_prepare_liftoff_parsed_summary.csv
│   │   ├── mikado_prepare_liftoff_summary_stats.csv
│   │   └── mikado_prepare_liftoff_parsed_summary.ids.rejected.txt
│   └── mikado_prepare_lifton
│       ├── list.txt
│       ├── mikado_prepare_lifton.fasta
│       ├── mikado_prepare_lifton.gtf
│       ├── mikado_prepare_lifton.log
│       ├── mikado_prepare_lifton_parsed_summary.csv
│       ├── mikado_prepare_lifton_summary_stats.csv
│       └── mikado_prepare_lifton_parsed_summary.ids.rejected.txt
├── output_presence_absence_data.tsv
├── lifton.sorted.ids.retained.tsv
├── output_upset_plot.png
├── liftoff.sorted.ids.retained.tsv
├── lifton.sorted.retained.gff
├── lifton.sorted.retained.gff.mikado_stats.summary.tsv
├── lifton.sorted.retained.gff.mikado_stats.tsv
├── liftoff.sorted.retained.gff
├── liftoff.sorted.retained.gff.mikado_stats.summary.tsv
└── liftoff.sorted.retained.gff.mikado_stats.tsv
```


#### Using LiftOn GFF only
```console
$ cd /path/to/work_directory
$ liftclean \
    --genome_fasta input_genome.fasta \
    --lifton_gff input_lifton.gff

# Output will be saved to 'output' directory by default
output
├── genome.fasta -> /path/to/work_directory/input_genome.fasta
├── lifton.gff -> /path/to/work_directory/input_lifton.gff
├── genome.fasta.fai
├── output_Rejected_Transcripts_summary_plot.png
├── analysis
│   ├── lifton.corrected.gff
│   ├── lifton.corrected.gff_strand_checker.log
│   ├── lifton.short_introns.corrected.gff
│   ├── lifton.short_introns.mapping.tsv
│   ├── lifton.short_introns.tsv
│   ├── lifton.corrected.gffread.gff
│   ├── lifton.sorted.gff
│   ├── lifton.sorted.gff.log
│   ├── lifton.sorted.ids.retained.txt
│   ├── lifton.sorted.ids.tsv
│   └── mikado_prepare_lifton
│       ├── list.txt
│       ├── mikado_prepare_lifton.fasta
│       ├── mikado_prepare_lifton.gtf
│       ├── mikado_prepare_lifton.log
│       ├── mikado_prepare_lifton_parsed_summary.csv
│       ├── mikado_prepare_lifton_summary_stats.csv
│       └── mikado_prepare_lifton_parsed_summary.ids.rejected.txt
├── lifton.sorted.ids.retained.tsv
├── lifton.sorted.retained.gff
├── lifton.sorted.retained.gff.mikado_stats.summary.tsv
└── lifton.sorted.retained.gff.mikado_stats.tsv
```

#### Using Liftoff GFF only
```console
$ cd /path/to/work_directory
$ liftclean \
    --genome_fasta input_genome.fasta \
    --liftoff_gff input_liftoff.gff

# Output will be saved to 'output' directory by default
output
├── genome.fasta -> /path/to/work_directory/input_genome.fasta
├── liftoff.gff -> /path/to/work_directory/input_liftoff.gff
├── genome.fasta.fai
├── output_Rejected_Transcripts_summary_plot.png
├── analysis
│   ├── liftoff.corrected.gff
│   ├── liftoff.corrected.gff_strand_checker.log
│   ├── liftoff.short_introns.corrected.gff
│   ├── liftoff.short_introns.mapping.tsv
│   ├── liftoff.short_introns.tsv
│   ├── liftoff.corrected.gffread.gff
│   ├── liftoff.sorted.gff
│   ├── liftoff.sorted.gff.log
│   ├── liftoff.sorted.ids.retained.txt
│   ├── liftoff.sorted.ids.tsv
│   └── mikado_prepare_liftoff
│       ├── list.txt
│       ├── mikado_prepare_liftoff.fasta
│       ├── mikado_prepare_liftoff.gtf
│       ├── mikado_prepare_liftoff.log
│       ├── mikado_prepare_liftoff_parsed_summary.csv
│       ├── mikado_prepare_liftoff_summary_stats.csv
│       └── mikado_prepare_liftoff_parsed_summary.ids.rejected.txt
├── liftoff.sorted.ids.retained.tsv
├── liftoff.sorted.retained.gff
├── liftoff.sorted.retained.gff.mikado_stats.summary.tsv
└── liftoff.sorted.retained.gff.mikado_stats.tsv
```


### LiftClean Main Output Files
* `*sorted.retained.gff`: Cleaned GFF3 file containing only the transcripts that passed all filtering criteria.
* `*Rejected_Transcripts_summary_plot.png`: PNG plot summarising the types of warnings found in the rejected transcripts.
* `*sorted.retained.gff.mikado_stats.tsv`: Detailed statistics of the cleaned GFF3 file generated by Mikado.
* `*sorted.retained.gff.mikado_stats.summary.tsv`: Summary statistics of the cleaned GFF3 file generated by Mikado.
* `*upset_plot.png`: UpSet plot showing transcript IDs that are common/unique between the retained and rejected Liftoff and liftOn files. This is only generated when both Liftoff and liftOn GFFs are provided.

## License

MIT
