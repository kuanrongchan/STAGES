# Static and Temporal Analysis of Gene Expression Studies (STAGES)
![Untitled (4)](https://user-images.githubusercontent.com/91276553/141084154-7d84695a-b220-43c5-bd41-08a38fd0ec70.png)

STAGES is an easy-to-use web tool that integrates data visualisation and pathway enrichment analysis for both static and temporal gene expression studies. STAGES is free and open to all users and there is no login requirement. The web tool works by running the Python programming language at backend to perform the data analysis and graph plotting, while the Streamlit framework is used to display the output data tables and graphs at frontend. Overall, STAGEs allow users to perform the following:
    1.	Plot interactive volcano plots
    2.	Filter data to identify and quantify number of differentially expressed genes based on users’ pre-defined fold change and p-value cut-offs
    3.	Pathway analysis by Enrichr against Gene Ontology (GO) Biological Processes, GO Molecular Function. GO Cellular Component, Reactome databases. Also allows pathway enrichment analysis against customised gene sets such as the Blood Transcriptomic Modules (BTM) and our in-house curated database (Vaccinomics)
    4.	GSEA pre-ranked analysis against the Reactome database, BTM and Vaccinomics
    5.	Plot clustergrams based on DEGs, genes from Enrichr pathway enrichment analysis or leading edge genes from GSEA
		6.  STRING query based on DEGs or user-input gene list.
    7.	Correlation matrix comparing transcriptomics responses between different experimental conditions.

## Getting started
To use the web tool, you will need at least one comparison file which contain:

    1.	Gene names on the first column
    2.	Ratio values (relative transcript expression versus control or baseline)
    3.	Adjusted p-values (or p-values)
    
    For the tool to be able to recognise your ratio and p-value columns, please format the ratio and p-value header columns to be parsed by underscores ( _ ):
    
    1.	Ratio as ratio_X_vs_Y
    2.	Adjusted p-values (or p-values) as pval_X_vs_Y,
    
    where X and Y are comparison variables parsed by underscores ( _ ). The X and Y variables can be time-point comparisons (e.g.ratio_hr6_vs_0) or experimental-control comparisons (e.g. ratio_drugA_vs_placebo, ratio_virus_vs_ctrl).
    
    For multiple comparisons to be analysed simultaneously, users can add more ratio and p-value columns (e.g. ratio_A_vs_Y, pval_A_vs_Y, ratio_B_vs_Y, pval_B_vs_Y). Users do not need to manually remove any other column labels and values that are present within the file.

Demo examples and descriptions of data formats are provided. You can try out the demo examples to familiarise yourself with the apps before uploading your dataset. Further documentation breaking down each feature of the webtool is available within the website's 'Read the Docs'.

## Data safety and security
The data you upload is safe and is never stored anywhere.

## Contributors
These apps are jointly made by myself (Kuan Rong Chan), Clara Koh, Justin Ooi and Gabrielle Lee from Duke-NUS, Department of Emerging Infectious Diseases. I am also thankful for Eugenia Ong and Ayesa Syenina from VIREMICS for their constructive feedback. These apps are now free for everyone to use, but for a limited period of time as we are constantly upgrading the apps. For more details on what we do, feel free to visit us at omicsdiary.com.
