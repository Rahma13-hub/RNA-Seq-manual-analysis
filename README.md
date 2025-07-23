# Functional Analysis of Top Expressed Genes in Lung Squamous Cell Carcinoma (LUSC)
This project offers a functional enrichment analysis for the top expressed genes from a training RNA-seq sample labeled lusc rsem, pertaining to Lung Squamous Cell Carcinoma (LUSC). This enrichment analysis was conducted using manual coding in R and subsequently with DAVID Functional Annotation Tool, aimed at biological meaning and cancer context.
#Project Overview
Data: Normalized RNA-seq expression data from a LUSC training sample
Gene Selection: Top 50 most highly expressed genes were selected for downstream analysis.
Tools: R for visualization, DAVID for functional annotation and enrichment.
Goal: To explore key biological processes and pathways linked to cancer behavior, especially in LUSC.
# Analyses Performed
# Differential Expression Visualization
- Volcano Plot: Highlights the most top significant genes.
- Heatmap: Shows expression patterns across selected genes.
  - ggplot Plots: Focused on the top enriched terms
-All visualizations were generated manually in R.
## Functional Enrichment
Using DAVID platform, we analyzed the top expressed genes to identify overrepresented biological terms and pathways.
Term                     | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| Cytoprotection           | Suggests cellular resistance mechanisms common in therapy-resistant tumors. |
| Signal                   | Includes signaling molecules important in tumor growth and communication.  |
| Cell-Matrix Adhesion     | Implicated in invasion and metastasis in solid tumors.            |
| Extracellular Region     | Involves secreted proteins shaping the tumor microenvironment.              |
| Angiogenesis             | Supports formation of new blood vessels — essential for tumor expansion.    |
These results are consistent with what is already known about LUSC progression mechanisms including evasion of apoptosis, extracellular remodeling, and angiogenesis associated.
##Files Included
volcano plot.png - Volcano plot of selected genes
Heatmap. png - Heatmap of gene expression
15 enriched terms. png - for the eriched terms
gene.list.txt - List of top expressed genes used in DAVID
DAVID_results.txt - Enrichment output from DAVID 
# R Script
The analysis pipeline including data loading, gene selection, heatmap and volcano plotting — is included in the file:
RNA-Seq manual script.R
# Biological Revelance
The LUSC gene signature also highlighted in this project indicates:
- Stress survival and anti-apoptotic phenotype (cytoprotection)
- Disassociation of adhesion and communication
- Pro-angiogenic signaling
- Factors secreted that modify the lung tumor microenvironment

These characteristics are recognized contributors to lung squamous carcinoma pathobiology, and provide evidence that using this training sample is a valid basis for pathway-based understanding of tumor biology.
# About the Project
This project was created as part of a personal initiative to develop skills in:
- RNA-seq interpretation
- Functional annotation using DAVID
- Manual visualization with R
  # Connect
 This project is part of my bioinformatics learning journey and academic portfolio.  
Feel free to explore, share feedback, or reach out for collaboration or discussion.

 [LinkedIn – Rahma Reda](http://www.linkedin.com/in/rahma-reda-2269a236b)
 
