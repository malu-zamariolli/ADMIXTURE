# Global Ancestry (ADMIXTURE)

# Usage

This Github contains the workflow pipeline that was used to run global ancestry using [Admixture.v1.3](https://dalexander.github.io/admixture/admixture-manual.pdf). The workflow is organized as a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. 

# Input data
For this pipeline, we need:
- Plink binary files for the sample of interest (Target Cohort) after QC
- 1000 Genomes reference samples with known populations (.vcf files)
- Configuration file used by the snakemake pipeline 

## Target Cohort: plink binaries
Plink binaries after quality control for the admixed population of interest

Files should be named in order to be recognized by defined wildcards in snakemake:

*{dataset}_QCed_final.bed -> testData_QCed_final.bed*

*{dataset}_QCed_final.bim -> testData_QCed_final.bim*

*{dataset}_QCed_final.fam -> testData_QCed_final.fam*

## 1000 Genomes Reference Sample
Admixture will estimate population proportions (K) for each individual. Adding 1000 Genome reference samples, allows us to distinguish the populations.

Information about 1000 Genomes data: https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html

- File with IDs and populations
*../1KGenomes/1000GP_Phase3/1000GP_Phase3.sample*

```console
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz 
tar -xf 1000GP_Phase3.tgz
```
- VCF files per chromosome in *../1KGenomes/Reference_vcf/*:
```console
nohup wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{1..5}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz{,.tbi} &

nohup wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{6..10}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz{,.tbi} &

nohup wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{11..16}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz{,.tbi} &

nohup wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{17..22}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz{,.tbi} &
```

## Configuration file (.yaml)
- File with modifiable parameters for the pipeline. 
- The file needs to be in the ADMIXTURE directory with snakefile. But it'll get copied (snakemake rule config) to the directory with output so that specific parameters of the run can be saved

*config_admixture.yaml*

```
# Date: # Fill with data (not needed in the pipeline)
dir: panel1_testData # Name of the directory that will be created to store output files (should identify different runs) (used as wildcard)
panel: panel1 # Given name of the reference population panel from 1000 Genomes (used as wildcard)
ref_samples: "GBR|CEU|TSI|IBS|PEL|MXL|PUR|CLM|ESN|GWD|LWK|MSL|YRI" #Popolutations to be selected from 1000 Genomes data
maf_ref: 0.01 #MAF applied for 1000 Genomes data filtering
# QC thr for reference panel #QC thresholds applied for 1000 Genomes reference panel
qcmaf: 0.01
qcgeno: 0.01
qcmind: 0.01
qchwe: 0.000001
# Target sample 
dataset: testData # Name of the target data (used as wildcard)
prune_window: 50 # Prunning parameters for target data (plink)
prune_step: 10
prune_r2: 0.1
# Admixture 
Kstart: 3 #First K to be tested in admixture
Kend: 14 #Last K to be testes in admixture - Plus one from the K of interest, in this case 13
# In the case of this config file, admixture will run from K 3 to 13
```

# Tools
All tools are created by the pipeline as conda environments
- [Plink.v1.9](https://www.cog-genomics.org/plink/):  *../envs/plink.yaml*
- [BCFtools.v.1.18](https://samtools.github.io/bcftools/bcftools.html): *../envs/bcftools.yaml*
- [Admixture.v1.3](https://dalexander.github.io/admixture/admixture-manual.pdf):  *../envs/admixture.yaml*
- [R.v.4.2](https://www.r-project.org):  *../envs/R.yaml* and *../envs/Rplotting.yaml*

# References
Admixture plot was adapted from https://github.com/mkanai/cancritls/blob/master/admixture/admixture.plot.r
