# Freyria_et_al_in_revision

Scripts for "**Salinity tolerance and stress responses in an Arctic Pelagophyte revealed by comparative transcriptomic and gene expression analysis**" 
by **Nastasia J. Freyria, Alan Kuo, Mansi Chovatia, Jenifer Johnson, Anna Lipzen, Kerrie W. Barry, Igor V. Grigoriev and Connie Lovejoy**.

All transcriptome samples are in the DOE JGI Genome Portal under Sequencing Project ID 1253386 and Analysis Project ID 123385, from the Sequence Read Archive (SRP284677-SRP284681). CCMP 2097 reference genome and the annotated genome are available at JGI Genome Portal and PhycoCosm Portal under JGI Project ID 1020062.

## Transcriptomic analysis

### Pipeline steps
1. Quality control of clean sequences (available on DOE JGI Genome Portal) using FastQC/0.11.9
```
zcat *fastq.gz | fastqc stdin
```
2. Trimmomatic/0.36
```
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 8 -phred33 ../R1.fastq.gz ../R2.fastq.gz R1.paired.fastq.gz R1.unpaired.fastq.gz R2.paired.fastq.gz R2.unpaired.fastq.gz ILLUMINACLIP:trimmomatic/0.36/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:100
```
3. Mapping against genome reference using STAR/2.6.1a
```
# STAR index
STAR --runMode genomeGenerate --genomeDir /Genome_directory/ --genomeFastaFiles /Genome_directory/Pelago2097_1_AssemblyScaffolds_Repeatmasked.fasta --runThreadN 24 --sjdbGTFfile /Genome_directory/Pelago2097_1_GeneCatalog_20160408.gff3 --sjdbOverhang 99 --sjdbGTFtagExonParentTranscript Parent

# Mapping
STAR --genomeDir /Genome_directory/ --readFilesIn /Trimmomatic/R1_paired.fq.gz /Trimmomatic/R2_paired.fq.gz --runThreadN 8 --readFilesCommand zcat --outFileNamePrefix Star_results/Star_mapp_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --quantMode GeneCounts --outReadsUnmapped Fastx --sjdbGTFfile /Genome_directory/Pelago2097_1_GeneCatalog_20160408.gtf
```
4. FeatureCounts
```
featureCounts -T 4  -a /Genome_directory/Pelago2097_1_GeneCatalog_20160408.gtf -o Results_featureCounts_matrix/featurecounts.txt *.out.bam >  Results_featureCounts_matrix/featurecounts_screen_output.log
```
5. Differential gene expression analysis
- R analyses and R figures were run using the master script "**Script_R_analysis_figures_data_Freyria_et_al.R**" in the folder.
- Figures 4c, 4d, 5a and 6a were made using EXCEL/16.52
- Figures 5b and 5c were made using POWERPOINT/16.52
- Figure 6b was made using MUSCLE/3.8.1551 and RAXML/8.2.11 with the following command:
````
# Alignment using MUSCLE
muscle -in file.fa -out file_muscle.fa

# Tree was constructed using RAXML
raxmlHPC-HYBRID-AVX2 -f a -#1000 -m PROTGAMMAGTR -p 12345 -x 12345 -s file_muscle.fa -n file_muscle.fa_PROTGAMMAGTR.tree
```
