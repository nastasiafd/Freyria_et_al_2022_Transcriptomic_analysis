# Freyria_et_al_in_revision

Scripts for "**Salinity tolerance and stress responses in an Arctic Pelagophyte revealed by comparative transcriptomic and gene expression analysis**" 
by **Nastasia J. Freyria, Alan Kuo, Mansi Chovatia, Jenifer Johnson, Anna Lipzen, Kerrie W. Barry, Igor V. Grigoriev and Connie Lovejoy**.

All transcriptome samples are in the DOE JGI Genome Portal under Sequencing Project ID 1253386 and Analysis Project ID 123385, from the Sequence Read Archive (SRP284677-SRP284681). CCMP 2097 reference genome and the annotated genome are available at JGI Genome Portal and PhycoCosm Portal under JGI Project ID 1020062.

## Transcriptomic analysis

### Pipeline steps
1. Quality control of raw sequences
2. Trimmomatic
3. Mapping against genome reference
4. FeatureCounts
5. Differential gene expression analysis

- R analyses and R figures were run using the master script "**Script_R_analysis_figures_data_Freyria_et_al.R**" in the Scripts folder.
- Figures 4c, 4d, 5a and 6a were made using EXCEL/16.52
- Figures 5b and 5c were made using POWERPOINT/16.52
- Figure 6b was made using MUSCLE/3.8.1551 and RAXML/8.2.11 with the following command:
  - muscle -in file.fa -out file_muscle.fa
  - raxmlHPC-HYBRID-AVX2 -f a -#1000 -m PROTGAMMAGTR -p 12345 -x 12345 -s file_muscle.fa -n file_muscle.fa_PROTGAMMAGTR.tree
