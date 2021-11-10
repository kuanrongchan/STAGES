# Static and Temporal Analysis of Gene Expression Studies (STAGES)
![Untitled (4)](https://user-images.githubusercontent.com/91276553/141084154-7d84695a-b220-43c5-bd41-08a38fd0ec70.png)
STAGES is a multi-app that integrates data visualisation and pathway analysis for static and temporal gene expression studies. STAGES is an open source and community funded web tool for creating beautiful charts from gene expression datasets. The multi-page web app built using Streamlit, which currently allows users to:
1. Plot interactive volcano plots
2. Identify differentially expressed genes (Users can apply their preferred fold-change and p-value cut-offs to identify number of upregulated and downregulated genes)
3. Build customised clustergrams based on either identified up-regulated DEGs (UP) and down-regulated DEGs (DOWN)
4. Build customised clustergrams based on selected gene list
5. Plot interactive correlation matrix comparing across different time-points or experimental conditions

## Getting started
To use the app, you will need one comparison file which should minimally contain:
1. Gene names on the first column
2. Ratio values (relative transcript expression versus control or baseline)
3. Adjusted p-value (or p-value)

For the app to be able to recognise your ratio and p-values, please label:
1. Ratio as ratio_X_vs_Y
2. Adjusted p-values (or p-values) as pval_X_vs_Y,

where X and Y are the comparison variables. 

Some examples include: ratio_virus_vs_ctrl, ratio_drugA_vs_placebo, ratio_hr6_vs_0, ratio_day1_vs_day0. 

For multiple comparisons to be made, simply insert more comparison columns (e.g. ratio_A_vs_Y, pval_A_vs_Y, ratio_B_vs_Y, pval_B_vs_Y ...), but please ensure that  "Y" is consistently present in all comparisons. Also, ensure that no icons or symbols are in "X" and "Y"





The second file, also called the comparisons file, should contain the fold-change, p-value, adjusted p-values and ratio values. These values can be easily obtained from specialised data processing software such as such as Partek Genomics Suite, Limma package, Transcriptomic Analysis Console and Nanostring analysis. Typically, for time-course studies, these values are obtained by comparing a specific time (e.g. days can be annotated as 1d, 3d and 7d for day 1, 3 and 7 respectively; hours can be annotated as 1hr, 3hr, 7hr) versus time 0. You may also choose to compare experimental groups versus control/placebo group. The comparisons file is required for plotting volcano plots, cumulative DEG plots and for pathway enrichment analysis. The first column should be gene names or gene symbols, and subsequent columns should contain the fold-change, p-value, adjusted p-values and ratio values.
