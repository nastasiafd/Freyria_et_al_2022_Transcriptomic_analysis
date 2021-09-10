# Transcriptomic_analysis

Freyria_et_al_in_revision

Scripts for "Salinity tolerance and stress responses in an Arctic Pelagophyte revealed by comparative transcriptomic and gene expression analysis" by Nastasia J. Freyria, Alan Kuo, Mansi Chovatia, Jenifer Johnson, Anna Lipzen, Kerrie W. Barry, Igor V. Grigoriev and Connie Lovejoy.

R analyses and R figures were run using the master script "Script_R_analysis_figures_data_Freyria_et_al.R" in the Scripts folder.

Figures 4c, 4d, 5a and 6a were made using excel

Figures 5b and 5c were made using powerpoint

Figure 6b was made using MUSCLE and RAXML/8.2.11 with the following command:
- muscle -in file.fa -out file_muscle.fa
- raxmlHPC-HYBRID-AVX2 -f a -#1000 -m PROTGAMMAGTR -p 12345 -x 12345 -s file_muscle.fa -n file_muscle.fa_PROTGAMMAGTR.tree
