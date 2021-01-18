# PAPAA: A galaxy enabled suite of tools for PanCancer Aberrant Pathway Activity Analysis

**Vijay Nagampalli and Daniel Blankenberg**

## Detecting aberrant PI3K activity in The Cancer Genome Atlas PanCancerAtlas

In this repository we adapted and modified a previously developed machine learning approach that utilized elastic net penalized logistic regression classification modeling to derive transcription signature or pathway alterations to measure aberrant RAS activity in The Cancer Genome Atlas (TCGA) [PanCancerAtlas project](https://gdc.cancer.gov/about-data/publications/pancanatlas). The logic and dataset processing steps for computing aberrant pathway activity analysis have been previously described in detail [Way et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29617658/). The link to original source [code](https://github.com/greenelab/pancancer) and changes that this work is based upon is [here](https://github.com/greenelab/pancancer/commit/d1b3de7fa387d0a44d0a4468b0ac30918ed66886).

We enabled this machine learning approach to measure aberrant pathway activity for any given combination of genes and diseases from TCGA. We named this collection of tools as PanCancer Aberrant Pathway Activity Analysis (PAPAA) suite. To facilitate accessibility and generalized reuse, we have made these software tools available through both command-line interfaces as well as graphical user interfaces via incorporation into the [Galaxy](https://usegalaxy.org/) platform. Galaxy is an open-source software project with large collections of bioinformatics tools that are made accessible through user-friendly web interfaces [Afgan et al., 2016](https://pubmed.ncbi.nlm.nih.gov/27137889/),[2018](https://pubmed.ncbi.nlm.nih.gov/29790989/). Users can generate, share, and easily run complex analysis pipelines. We developed a generalized pipeline using Galaxy to enable researchers to easily apply our analysis methodology to any given genes/cancer types found in TCGA, or on user-provided datasets, without the need for programming experience and computing environmental control. Additionally we created a Galaxy training Material to demonstrate the galaxy PAPAA tools and generated a PI3K- oncogenes based classifier to detect PI3K activation in TCGA cancer types.

In this repository, we provide scripts to run the tools in command line and also developed a conda package to install and run the tools for any given combination of genes and diseases listed in TCGA link is below. 



Using these tools, we designed classifiers featuring either oncogenes or tumor suppressors and reflecting either gain of function (activating) or loss of function (deactivating) mutations. We demonstrated a generic transcriptional signature difference between tumors with gain of function and loss of function mutations. Our models can predict aberrant PI3K and drug sensitivity among GDSC and CCLE common cell lines (unpublished data)

The code in this repository is flexible and can build a Pan-Cancer classifier for any combination of genes and cancer-types using gene expression, mutation, and copy number data.

## Open Access Data

All data was released by the TCGA PanCancerAtlas project. The logic and dataset processing steps for computing aberrant pathway activity analysis have been previously described in detail in [Way et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29617658/). The data files used in this analyses presented here are archived on [PAPAA-Zenodo](https://zenodo.org/record/4306639#.X9FJF-lKgZE)

## Installing from bioconda repository:
papaa can be installed from [bioconda](https://anaconda.org/bioconda/papaa)

```
$ conda install -c bioconda papaa
```

## Installing PAPAA tools in your Galaxy instance:

PAPAA tool suite (suite_papaa) can be downloaded and installed through [galaxy toolshed](https://toolshed.g2.bx.psu.edu/)  as a standard galaxy tool in your own galaxy instance.  

## Prebuild Docker Image: 

A prebuild build docker image based on the recent galaxy release can be obtained by the link below for a quick installation. 

[Docker Image](https://github.com/nvk747/galaxy_papaa/) 

## Galaxy training material:

Tutorial for running these tools and for building classifier using oncogenes ERBB2,KRAS,PIK3CA,ATK1 was created and accessible in Galaxy training materials check for the Aberrant pi3k pathway analysis is [here](https://training.galaxyproject.org/training-material/topics/statistics/).
