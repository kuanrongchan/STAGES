# Static and Temporal Analysis of Gene Expression Studies (STAGES)
![Untitled (4)](https://user-images.githubusercontent.com/91276553/141084154-7d84695a-b220-43c5-bd41-08a38fd0ec70.png)
STAGES is a multi-app that integrates data visualisation and pathway analysis for static and temporal gene expression studies. STAGES is an open source and community funded web tool for creating beautiful charts from gene expression datasets. The multi-page web app built using Streamlit, which currently allows users to analyse an omics data in distinct stages:
1. Plot interactive volcano plots
2. Identify differentially expressed genes (Users can apply their preferred fold-change and p-value cut-offs to identify number of upregulated and downregulated genes)
3. Build customised clustergrams based on identified up-regulated DEGs (UP) or down-regulated DEGs (DOWN)
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

Some examples of labelling "X" include: ratio_virus_vs_ctrl, ratio_drugA_vs_placebo, ratio_hr6_vs_0, ratio_day1_vs_day0. 

Some examples of labelling "Y" include: pval_virus_vs_ctrl, pval_drugA_vs_placebo, pval_hr6_vs_0, pval_day1_vs_day0. 

For multiple comparisons to be made, simply insert more comparison columns (e.g. ratio_A_vs_Y, pval_A_vs_Y, ratio_B_vs_Y, pval_B_vs_Y ...), but please ensure that  "Y" is consistently present in all comparisons. Also, ensure that no icons or symbols used for labelling "X" and "Y." If you have other column statistics, it is not necessary to remove them.

Demo examples and descriptions of data formats are provided. You can try out the demo examples to familiarise yourself with the apps before uploading your dataset

## Data safety and security
The data you upload is safe and is never stored anywhere.

## Contributors
These apps are jointly made by myself (Kuan Rong Chan), Clara Koh, Justin Ooi and Gabrielle Lee from Duke-NUS, Department of Emerging Infectious Diseases. I am also thankful for Eugenia Ong and Ayesa Syenina from VIREMICS for their constructive feedback. These apps are now free for everyone to use, but for a limited period of time as we are constantly upgrading the apps. For more details on what we do, feel free to visit us at omicsdiary.com.
