# Pareto-Front-analysis-for-scRNA-seq-data

I worked on this computational biology project using single cell RNA expression data from metastatic lung adenocarcinoma and performed Pareto front analysis, and differential expression analysis on the genes expressed near the archetypes for each tumor sample to explore the intra tumor heterogeneity (ITH). I was provided with a TPM normalized single cell RNA sequencing dataset to perform the Pareto front analysis using the ParTI MATLAB package developed by Uri Alon group. I then performed differential expression analysis on the genes expressed near the archetypes from Pareto analysis, using various Bioconductor packages in R.

Tumors show a wide range of phenotypic and genetic characteristics at both intertumor and intratumor levels. Intratumor heterogeneity, also called as intralesion heterogeneity, describes the presence of several tumor cell populations within the same tumour tissue (Jamal-Hanjani et al., 2015).  Recent evidences have shown that heterogeneity in tumor is not only because of genetic factors, but also encompasses epigenetic, transcriptional, phenotypic, metabolic, and secretory components (Marusyk, Almendro & Polyak,  2012; Ramón Y Cajal et al., 2020; Teixeira et al., 2019; Sharma et al., 2019; van Galen et al., 2019; Grosselin et al., 2019). Most tumors are now understood to be complex ecosystems that originate and evolve in response to strong selective pressure from their microenvironment, which includes trophic, metabolic, immunological, and therapeutic components and these pressure promote the diversification malignant and non-malignant tumors in the tumor microenvironment (Vitale et al., 2019).

Intratumor heterogeneity (ITH) is believed to be one of the most important factors of therapeutic resistance and treatment failure, as well as one of the primary causes of poor overall survival in cancer patients with metastatic illness (Jamal-Hanjani et al., 2015; Jamal-Hanjani et al., 2017). Tumors are made up of a mosaic of cancer cells with varied features and susceptibilities to anticancer drugs, and this heterogeneity is very important for tumor progression, malignancy, metastasis and drug resistance. Research in multiregion genome-sequencing has found considerable variations in the genetic makeup of malignant cells in distinct regions of the same tumor, called as spatial ITH. Similarly, results from longitudinal studies have shown that genetic features of the tumor can vary over time and are called temporal ITH (Wang et al., 2016). The tumor mutational landscape and evolutionary trajectory that have been defined by single or multisample sequencing investigations have revealed extensive genetic ITH in both spatial and temporal dimensions (Figure 1). There are different layers of complexity in this tumor. Since tumor harbors such heterogeneity, individual patients, lesions, and cell populations should be characterized before treatment to achieve the goals of precision medicine (McGranahan et al., 2015; McGranahan & Swanton, 2017).
<p align="center">
  <img width = "390" alt="image" src="https://user-images.githubusercontent.com/121320501/236418184-a8f136b9-8a51-45cf-b0ec-21889c41c9f1.png"> 
  
  Figure 1: Clonal heterogeneity and cellular consortium (Source: Ramón Y Cajal et al., 2020)
</p>

All biological systems that perform multiple evolutionary tasks cannot be optimal in all tasks and face a trade-off. The implementation of Pareto front concept in biological datasets can be used for finding the best-trade-offs of phenotypes using the weighted averages of archetypes. Here the archetypes are phenotypes optimal for each task. The Pareto front describes the optima for fitness functions that are increasing functions of each task. This Pareto optimality analysis predicts that the multiple tasks performed by cells or organisms have phenotypes that fall on the low dimensional polytopes with the optima for the task being the archetypes. The two assumptions for the Pareto front analysis are: (1) A single phenotype called archetype has the maximized (optimal)  performance function, and (2) the performance of the function decreases with the distance from the archetype. The tasks being performed can be identified and inferred from the phenotypes expressed near the archetypes (Shovel et al., 2012).

The Pareto task inference method (ParTI) was developed by the Uri Alon group (http://www.weizmann.ac.il/mcb/UriAlon/download/ParTI) as a Matlab package for inferring biological tasks from high-dimensional biological data. For analyzing the high dimensionality data, dimensionality reduction techniques like principal-component analysis (PCA) (Markus Ringnér, 2008), t-distributed stochastic neighbor embedding (t-SNE) (van der Maaten & Hinton, 2008) and methods that split data points into groups, such as clustering and Gaussian mixture models (GMMs) (Hastie, Tibshirani & Friedman, 2009). The Pareto approach suggests that if cells have to perform multiple tasks, then no gene expression can be optimal for all the tasks faced by the cell. A biological data can be represented by a Pareto-optimal situation when (1) the data falls inside s polytope (2) the data near the vertex of the polytope (archetype) correspond to specific biological tasks. To implement the ParTI tool, a set of N datapoint, described as vector of K traits, annotated by a vector of M additional features is provided as input. The two stages employed in the method are computing (1) the minimal polytope enclosing the data and its statistical significance, and (2) enriching of each feature as a function of the distance from each archetype. The number of archetypes are determined by fitting the polytopes with n vertices to the data using principal convex hull analysis (PCHA) (Mørup & Hansen, 2012). The data dimensionality is first reduced by principal comment analysis (PCA), and then for each number of archetypes n, the best fit polytope is found using PCHA algorithm. Based on the ‘elbow’ test, the n value is chosen beyond which there is an improvement in explained variance (EV). The archetype positions are then determined by hyper spectral unmoving algorithms.  <img width="216" alt="Screen Shot 2023-05-05 at 2 41 43 PM" src="https://user-images.githubusercontent.com/121320501/236419637-b9b93fbf-d841-4972-9aab-ae7fc2aeb650.png">

Where pi is the with datapoint and si is the closest point to pi in the polytope. And for points inside the polytope, <img width="95" alt="Screen Shot 2023-05-05 at 2 43 23 PM" src="https://user-images.githubusercontent.com/121320501/236419994-71461f61-0769-4253-9597-02397cf83188.png">

The data near the archetypes (closest to the vertices of the polytope) denote the features maximally enriched and allow the identification of the tasks the vertices represent. Hart et al, 2015, have depicted human breast tumors and mouse tissues by tetrahedrons in gene expression space, with specific tumor types and biological functions enriched at the vertices, suggesting four key tasks (Hart et al., 2015).

<p align="center">
  <img width = "500" alt="image" src="https://user-images.githubusercontent.com/121320501/236420484-d8535d94-5642-4c43-af10-e952c9a75b6e.png"> 
  
  Figure 2: Maximally enriched cancer feature near the archetypes(Source: Hart et al., 2015)
</p>

The TPM data from single-cell RNA sequencing lung adenocarcinoma (LUAD) was used for the Pareto front analysis (accession code GSE131907). The data contained single-cell mRNA expression profiles acquired from 58 samples from forty-four patients diagnosed with lung adenocarcinoma (LUAD). The TPM expression dataset was separated into smaller datasets based on the sample IDs using Python code. 

# splitting.py or multiprocess code.py can be used for this step. (The multiprocess code.py code uses 8 forks and can be modified depending on the users processor)

The large TPM dataset, separated based on the patient sample ID for matched cancer samples were analyzed by using the ParTI_lite function in the ParTI Matlab package. 

The ParTI_lite function executed is as follows:

<img width="404" alt="Screen Shot 2023-05-05 at 2 54 21 PM" src="https://user-images.githubusercontent.com/121320501/236422166-9ddc0bbb-3d98-4b6f-bedf-b0d3708fc874.png">

Where, 
1. LUNGN06 is the TPM cancer expression dataset 
2. 5 is the PCHA algorithm selected for finding the simplex
3. 10 is the dimension for Explained Variance (EV) to be calculated
4. DiscFeatName is a 1x2 string array with values "Index" &  “Cell_type.refined"
5. LUNGN06META is the summary data provided as a string cell array in which columns represent discrete feature and rows represent datapoints
6. 0 is set as value for the index of the columns
7. An empty matrix [] is provided for Continuous Feature labels
8. An empty matrix [] is provided for EnMatCont, a real matrix with attributes
9. An empty matrix [] is provided for GOcat2Genes
10. 'LUNGN06' is the output name of the files to be saved after the Pareto front analysis

The code executed is provided in the ParTI_litesri.m file. The default number of archetypes generated by the elbow method was used for constructing the Pareto front plot. Additional codes were included in the function for finding archetypes in the ParTI_lite code - findArchetypesnew2 provided in findArchetypesnew2.m . The additional lines introduced were for representing the datapoint in the plot with different colors based on the cell type the each datapoint belongs to, and to store the dimensionality reduced PCA points as an additional CSV file for further analysis.

The barcodes of the points for the samples near the archetypes were extracted from the dimensionality reduced PCA plots obtained after the Pareto front analysis by using Euclidean distance formula for the points from a fixed distance from the archetypes using a Matlab code and saving the points as a CSV file, and then comparing the obtained CSV with the original expression dataset to extract only those points expressed near the archetypes using a Python code. The CSV files stored were further subjected to differential gene expression analysis using R/bioconductor packages.

The expression dataset for points near each archetype for each of the cancer sample in the matched cancer samples, stored as a CSV, were processed for pairwise differential expression analysis using various packages in R/bioconductor like MAST (McDavid, Finak, & Yajima, 2021), DESeq2 (Love, Huber, & Anders, 2014) and DESingle (Miao et al., 2018). These packages mentioned are the most commonly employed packages for differential expression analysis with reliable results, with MAST being employed for single-cell RNA sequencing TPM datasets (Wang et al., 2019). The other packages were used in order to compare the results and look for correlation. The packages were executed by considering the CSV files containing the datapoints extracted near the archetypes for each sample type and run as a nested loop with the preceding archetypes being considered as the reference archetype against the succeeding archetypes (i.e. for archetype1 vs 2, archetype 1 was considered as reference condition). 

The Pearson and Spearman correlation for the results using different R packages were obtained by using R code Correlation.R

The results for the entire analysis can be availble on request. Please contact srisruthi02.ss@gmail.com






References:
Grosselin, K., Durand, A., Marsolier, J., Poitou, A., Marangoni, E., Nemati, F., Dahmani, A., Lameiras, S., Reyal, F., Frenoy, O., Pousse, Y., Reichen, M., Woolfe, A., Brenan, C., Griffiths, A. D., Vallot, C., & Gérard, A. (2019). High-throughput single-cell ChIP-seq identifies heterogeneity of chromatin states in breast cancer. Nature genetics, 51(6), 1060–1066. https://doi.org/10.1038/s41588-019-0424-9
Hart, Y., Sheftel, H., Hausser, J., Szekely, P., Ben-Moshe, N. B., Korem, Y., Tendler, A., Mayo, A. E., & Alon, U. (2015). Inferring biological tasks using Pareto analysis of high-dimensional data. Nature methods, 12(3), 233–235. https://doi.org/10.1038/nmeth.3254
Hastie, T., Tibshirani, R., & Friedman, J. H.(2009). The elements of statistical learning: data mining, inference, and prediction. 2nd ed. New York: Springer.520–528.
Jamal-Hanjani, M., Quezada, S. A., Larkin, J., & Swanton, C. (2015). Translational implications of tumor heterogeneity. Clinical cancer research : an official journal of the American Association for Cancer Research, 21(6), 1258–1266. https://doi.org/10.1158/1078-0432.CCR-14-1429
Jamal-Hanjani, M., Wilson, G. A., McGranahan, N., Birkbak, N. J., Watkins, T., Veeriah, S., Shafi, S., Johnson, D. H., Mitter, R., Rosenthal, R., Salm, M., Horswell, S., Escudero, M., Matthews, N., Rowan, A., Chambers, T., Moore, D. A., Turajlic, S., Xu, H., Lee, S. M., … TRACERx Consortium (2017). Tracking the Evolution of Non-Small-Cell Lung Cancer. The New England journal of medicine, 376(22), 2109–2121. https://doi.org/10.1056/NEJMoa1616288
Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
Marusyk, A., Almendro, V., & Polyak, K. (2012). Intra-tumour heterogeneity: a looking glass for cancer?. Nature reviews. Cancer, 12(5), 323–334. https://doi.org/10.1038/nrc3261
McDavid, A., Finak, G., Yajima, M. (2021). MAST: Model-based Analysis of Single Cell Transcriptomics. R package version 1.18.0, https://github.com/RGLab/MAST/.
McGranahan, N., & Swanton, C. (2017). Clonal Heterogeneity and Tumor Evolution: Past, Present, and the Future. Cell, 168(4), 613–628. https://doi.org/10.1016/j.cell.2017.01.018
McGranahan, N., Favero, F., de Bruin, E. C., Birkbak, N. J., Szallasi, Z., & Swanton, C. (2015). Clonal status of actionable driver events and the timing of mutational processes in cancer evolution. Science translational medicine, 7(283), 283ra54. https://doi.org/10.1126/scitranslmed.aaa1408
Miao, Z., Deng, K., Wang, X., & Zhang, X. (2018). DEsingle for detecting three types of differential expression in single-cell RNA-seq data. Bioinformatics (Oxford, England), 34(18), 3223–3224. https://doi.org/10.1093/bioinformatics/bty332
Mørup, M., & Hansen, L.K. (2012). Archetypal analysis for machine learning and data mining. Neurocomputing, 80, 54-63.
Ramón Y Cajal, S., Sesé, M., Capdevila, C., Aasen, T., De Mattos-Arruda, L., Diaz-Cano, S. J., Hernández-Losa, J., & Castellví, J. (2020). Clinical implications of intratumor heterogeneity: challenges and opportunities. Journal of molecular medicine (Berlin, Germany), 98(2), 161–177. https://doi.org/10.1007/s00109-020-01874-2
Ringnér M. (2008). What is principal component analysis?. Nature biotechnology, 26(3), 303–304. https://doi.org/10.1038/nbt0308-303
Sharma, A., Merritt, E., Hu, X., Cruz, A., Jiang, C., Sarkodie, H., Zhou, Z., Malhotra, J., Riedlinger, G. M., & De, S. (2019). Non-Genetic Intra-Tumor Heterogeneity Is a Major Predictor of Phenotypic Heterogeneity and Ongoing Evolutionary Dynamics in Lung Tumors. Cell reports, 29(8), 2164–2174.e5. https://doi.org/10.1016/j.celrep.2019.10.045
Shoval, O., Sheftel, H., Shinar, G., Hart, Y., Ramote, O., Mayo, A., Dekel, E., Kavanagh, K., & Alon, U. (2012). Evolutionary trade-offs, Pareto optimality, and the geometry of phenotype space. Science (New York, N.Y.), 336(6085), 1157–1160. https://doi.org/10.1126/science.1217405
Teixeira, V. H., Pipinikas, C. P., Pennycuick, A., Lee-Six, H., Chandrasekharan, D., Beane, J., Morris, T. J., Karpathakis, A., Feber, A., Breeze, C. E., Ntolios, P., Hynds, R. E., Falzon, M., Capitanio, A., Carroll, B., Durrenberger, P. F., Hardavella, G., Brown, J. M., Lynch, A. G., Farmery, H., … Janes, S. M. (2019). Deciphering the genomic, epigenomic, and transcriptomic landscapes of pre-invasive lung cancer lesions. Nature medicine, 25(3), 517–525. https://doi.org/10.1038/s41591-018-0323-0
van der Maaten, L., & Hinton, G.,. Visualizing data using t-SNE. J Mach Learn Res. 2008;9:2579–2605.
van Galen, P., Hovestadt, V., Wadsworth Ii, M. H., Hughes, T. K., Griffin, G. K., Battaglia, S., Verga, J. A., Stephansky, J., Pastika, T. J., Lombardi Story, J., Pinkus, G. S., Pozdnyakova, O., Galinsky, I., Stone, R. M., Graubert, T. A., Shalek, A. K., Aster, J. C., Lane, A. A., & Bernstein, B. E. (2019). Single-Cell RNA-Seq Reveals AML Hierarchies Relevant to Disease Progression and Immunity. Cell, 176(6), 1265–1281.e24. https://doi.org/10.1016/j.cell.2019.01.031
Vitale, I., Sistigu, A., Manic, G., Rudqvist, N. P., Trajanoski, Z., & Galluzzi, L. (2019). Mutational and Antigenic Landscape in Tumor Progression and Cancer Immunotherapy. Trends in cell biology, 29(5), 396–416. https://doi.org/10.1016/j.tcb.2019.01.003
Wang, J., Cazzato, E., Ladewig, E., Frattini, V., Rosenbloom, D. I., Zairis, S., Abate, F., Liu, Z., Elliott, O., Shin, Y. J., Lee, J. K., Lee, I. H., Park, W. Y., Eoli, M., Blumberg, A. J., Lasorella, A., Nam, D. H., Finocchiaro, G., Iavarone, A., & Rabadan, R. (2016). Clonal evolution of glioblastoma under therapy. Nature genetics, 48(7), 768–776. https://doi.org/10.1038/ng.3590
Wang, T., Li, B., Nelson, C. E., & Nabavi, S. (2019). Comparative analysis of differential gene expression analysis tools for single-cell RNA sequencing data. BMC bioinformatics, 20(1), 40. https://doi.org/10.1186/s12859-019-2599-6

