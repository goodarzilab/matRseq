# matRseq
mature tRNA sequencing

### Setting up the environment
Use the following command to create a conda environment with the essential packages used by matRseq.
```bash
conda env create --file=matRseq.yaml
conda activate matRseq
```
### Run matRseq.
#### Create the metadata file
The metadata file resides in the working directory and lists the required information for each sample. 
For example:

| Sample.name  | Sample.prefix  | R1  | R2  | sample.type | cell.line | sample.rep |
|---|---|---|---|---|---|---|
SW480_LVM2_ribo_r1 | SW480_LVM2_ribo_6_S15 | SW480_LVM2_ribo_6_S15_R1_001.fastq.gz | SW480_LVM2_ribo_6_S15_R2_001.fastq.gz | ribo |	SW480-Par |	1 |
SW480_LVM2_ribo_r2 | SW480_LVM2_ribo_7_S16 | SW480_LVM2_ribo_7_S16_R1_001.fastq.gz | SW480_LVM2_ribo_7_S16_R2_001.fastq.gz | ribo | SW480-Par | 2 | 
SW480_P_total_r1 | SW480_P_total_6_S9 | SW480_P_total_6_S9_R1_001.fastq.gz | SW480_P_total_6_S9_R2_001.fastq.gz	total |	SW480-LvM2 | 1 |
SW480_P_total_r2 | SW480_P_total_7_S10 | SW480_P_total_7_S10_R1_001.fastq.gz | SW480_P_total_7_S10_R2_001.fastq.gz | total | SW480-LvM2 | 2 |

##### human mature tRNA reference included in the package

#### Run the analysis
For bi-variate analysis the following command will run the analysis:
```bash
python matRseq.py --runMode metadata.txt 'sample.type~cell.line'  tRNA_LvM2-v-Par.txt
```
and for uni-variate analysis the following command will run the analysis:
```bash
python matRseq.py --runMode --ref=SW480Par metadata_univariate.txt '~cell.line'  tRNA_LvM2-v-Par.txt
```

#### Options
Run `python matRseq.py` for usage.
The following are the options:
1. `--runMode` or `-r` vs. `--printMode` or `-p` : `--printMode` prints all the commands that are run
2. `-a` or `--aligner`: Choose aligner package (default is BWA)
3. `-l` or `--read` : R1, R2, both; Use both R1 and R2 (miSeq) as opposed to R2 (10X). Default is R2
4. `--hasUMI` or `--noUMI`: Whether reads contain UMI or not (Default is hasUMI)
5. `--paired` or `--single` : if the reads are paired end. (Default is paired)
6. `--ref`: The sample that is the reference in univariate analysis.
