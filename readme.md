# PAPAA: A galaxy enabled suite of tools for PanCancer Aberrant Pathway Activity Analysis

**Vijay Nagampalli and Daniel Blankenberg**

## Detecting aberrant PI3K activity in The Cancer Genome Atlas PanCancerAtlas

In this repository we adapted and modified a previously developed approach that utilized elastic net penalized logistic regression classification modeling to derive transcription signature or pathway alterations to measure aberrant RAS activity in The Cancer Genome Atlas (TCGA) PanCancerAtlas project (link). The logic and dataset processing steps for computing aberrant pathway activity analysis have been previously described in detail (Way et al., 2018). 

We enabled this machine learning approach to measure aberrant pathway activity for any given combination of genes and diseases from TCGA. We named this collection of tools as PanCancer Aberrant Pathway Activity Analysis (PAPAA) suite. To facilitate accessibility and generalized reuse, we have made these software tools available through both command-line interfaces as well as graphical user interfaces via incorporation into the Galaxy platform. Galaxy is an open-source software project with large collections of bioinformatics tools that are made accessible through user-friendly web interfaces (Afgan et al., 2016, 2018). Users can generate, share, and easily run complex analysis pipelines. We developed a generalized pipeline using Galaxy to enable researchers to easily apply our analysis methodology to any given genes/cancer types found in TCGA, or on user-provided datasets, without the need for programming experience and computing environmental control. Additionally we created a Galaxy training Material to demonstrate the galaxy PAPAA tools and generated a PI3K- oncogenes based classifier to detect PI3K activation in TCGA cancer types.


In this repository, we provide directions for running this suite of tools in GALAXY. of each for building classifiers to detect aberration in _TP53_ and Ras signalling.We investigated gene expression states in the presence of PI3K oncogene mutations using transcriptome data from TCGA. for understanding transcriptomic -tutorial A transcriptome can describe the total state of a tumor at a snapshot in time.
In this repository, we use cancer transcriptomes from The Cancer Genome Atlas PanCancerAtlas project to interrogate gene expression states induced by deleterious mutations and copy number alterations.

The code in this repository is flexible and can build a Pan-Cancer classifier for any combination of genes and cancer-types using gene expression, mutation,
and copy number data.

## Open Access Data

All data was released by the TCGA PanCancerAtlas project. The logic and dataset processing steps for computing aberrant pathway activity analysis have been previously described in detail. The data files used in the analyses presented here are archived on[PAPAA-Zenodo](https://zenodo.org/record/4306639#.X9FJF-lKgZE)

## Usage

### Initialization:

First steps: install conda
==========================

Download the ``Miniconda`` installer 
from http://conda.pydata.org/miniconda.html and run it.   

If you already have conda, create a new environment papaa from bioconda.

```
conda create -n papaa -c bioconda papaa
source activate papaa
```

### Example Scripts
**Galaxy_training_material_link:**

**papaa conda package**

**papaa galaxy docker container** 

