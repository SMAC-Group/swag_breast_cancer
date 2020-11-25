Chameleon microRNAs in breast cancer: their elusive role as regulatory factors in cancer progression
====================================================================================================

Replication of the findings of the research article
---------------------------------------------------

In this repository we provide the codes necessary to replicate the findings presented in the research article **Chameleon microRNAs in breast cancer: their elusive role as regulatory factors in cancer progression** authored by Cesare Miglioli, Gaetan Bakalli, Stephane Guerrier, Samuel Orso, Roberto Molinari, Mucyo Karemera and Nabil Mili. This article is currently under review on *PLOS ONE*.

The statistical analysis performed in this study is based on the data presented in the paper **Subtype-specific micro-RNA expression signatures in breast cancer progression** by Haakensen in 2016. We thank the authors for having made available the AHUS data set on the free access ArrayExpress platform ().

Get the data from ArrayExpres in R
----------------------------------

First install the **ArrayExpress** package, through bioconductor, with the following code:

``` r
### R version 3.5 and older ###

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ArrayExpress", version = "3.8")

### R version 3.6 ###

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ArrayExpress")


require(ArrayExpress) #load the new package
```

Now you can create a temporary directory to store the raw data listed as *E-MTAB-3759* on **ArrayExpress**. Then you just need to extract the data thanks to the function *ae2bioc()* of the package. To have a first visual impression of the AHuS dataset, you can go to .

``` r
temp_dir <- "C:/..." #choose your preferred path

dir.create(temp_dir)
setwd(temp_dir)


# Download the raw files
# To run only the first time!

ae_obj <- ArrayExpress::getAE('E-MTAB-3759', type = 'full', path = temp_dir)

str(ae_obj) #to explore the new object

### Extraction ###

mtab3759raw <- ae2bioc(mageFiles = ae_obj)
```

The response *y* (i.e. if the subject has breast cancer or not) can be obtained directly from the *mtab3759raw* object together with the type of breast cancer (i.e. either benign, DCIS or invasive). However, in order to get the final design matrix *X* of miRNAs used in the study, we need to apply the function *normalizeBetweenArrays()* of the *limma* package to our previously found object. This function normalizes expression intensities so that the same intensities (or log-ratios) have similar distributions across a set of arrays.

``` r
# Response variable

y <- mtab3759raw@phenoData@data$Factor.Value.disease.

y_sub <- mtab3759raw@phenoData@data$Characteristics.clinical.history. 

# Design matrix X of miRNAs

mtab3759raw <- normalizeBetweenArrays(mtab3759raw,method="quantile")

X <- mtab3759raw@.Data[[1]]

X <-  t(X) #to have subjects on the rows and miRNAs on the columns
```

Reproduce the results
---------------------

Having now the final design matrix *X* and the response *y*, we are able to reproduce the results and the graphs contained in the research article. Please refer to the file `breast_cancer_data_analysis.R` in this repository for a detailed explanation.
