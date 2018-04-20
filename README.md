# crispr-mageck-nf
Statistical analysis of CRISPR and shRNA functional genetic screening data using MAGeCK.

## Installation

### Nextflow
Install `nextflow` following the instructions at https://www.nextflow.io/docs/latest/getstarted.html

### Singularity
Install `singularity` following the instructions at
https://singularity.lbl.gov/install-linux

### crispr-mageck-nf pipeline
The most convenient way is to install `crispr-mageck-nf` is to use `nextflow`'s built-in `pull` command
```bash
nextflow pull zuberlab/crispr-mageck-nf
```

## Test
Before you start, make sure `nextflow` and `singularity` are properly installed on your system.

Clone git repository from Github and run the pipeline using the provided test data.
```bash
git clone https://github.com/ZuberLab/crispr-mageck-nf.git
cd crispr-mageck-nf
./test
```

## Documentation
```bash
nextflow run zuberlab/crispr-mageck-nf --help
```

## Credits
Nextflow:  Paolo Di Tommaso - https://github.com/nextflow-io/nextflow

Singularity: Singularityware - https://singularity.lbl.gov

MAGeCK: Wei Li / Shirley Liu  - http://mageck.sourceforge.net
