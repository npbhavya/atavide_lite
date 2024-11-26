# Atavide

metagenome processing pipeline

## Install

### Source install 

Conda environment
```bash
conda create -y -n atavide python==3.12
conda activate atavide
```
Git clone 
```bash 
git clone https://github.com/npbhavya/atavide_lite.git
git branch nextflow 
cd atavide_lite
pip install . 
```

#### Install databases 

```bash
nextflow run atavide/workflow/install.nf
```

This downloads 
- human reference genome (Host Filtering step)
- mmseqs formatted UniRef50 database

#### Running atavide

```bash
#For running read based annotations without host filtering
nextflow run atavide/workflow/main.nf

#For running read based annotation with host filtering (default to search against human GRCh38)
nextflow run atavide/workflow/main.nf --hostFilter true

#For running the MAG worklow
nextflow run atavide/workflow/mag.nf

#To resume a workflow, just add -resume
```

Notes: Right now input and output directories are hardcoded in main.nf (fix later)
