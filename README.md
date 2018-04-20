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

## Documentation
```bash
nextflow run zuberlab/crispr-mageck-nf --help
```

## Credits
Nextflow:  Paolo Di Tommaso - https://github.com/nextflow-io/nextflow

Singularity: Singularityware - https://singularity.lbl.gov

MAGeCK: Wei Li / Shirley Liu  - http://mageck.sourceforge.net
