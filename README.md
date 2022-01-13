# exosomePurity
A tumour purity deconvolution model to estimate tumour purity in serum exosomes of cancer patients based on miRNA signatures

The steps of exosomePurity are as follows:


# Prerequisites
exosomePurity is written in R(version 4.0.3), requiring R packages DESeq2, preprocessCore, quadprog.



Preparing expresssion files of exosome cancer


#Usage

The steps of exosomePurity are as follows:

#Step1
Performing the differential analysis between cancer cell line-derived exosomes and healthy cell-derived exosomes using miRNA-Seq data
#Input data
The microRNA-seq expression files in format of rawcount.

#applying R script 01signature_DEGS.r

#output data
Differential gene between tumour and normal
