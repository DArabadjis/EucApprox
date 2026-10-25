# Experimental datasets: sources, attribution, and reuse terms

This document covers third-party datasets used in the experiments accompanying *EucApprox: Euclidean Embedding of Data Using Local Distances*. Experimental data and results are distributed separately from the [EucApprox code repository](https://github.com/DArabadjis/EucApprox).

## MNIST

**Dataset:** MNIST handwritten digit database.  
**Original attribution:** Yann LeCun and Corinna Cortes; derived from original NIST datasets.  
**Source:** https://yann.lecun.com/exdb/mnist/  
**License:** Creative Commons Attribution–ShareAlike 3.0 (CC BY-SA 3.0), as documented by the Keras MNIST dataset page: https://keras.io/api/datasets/mnist/  
**License terms:** https://creativecommons.org/licenses/by-sa/3.0/

The experimental archive contains MATLAB-formatted/flattened MNIST test-set data and separately computed experimental results. The original image and label data retain their source attribution and applicable license conditions. The format conversion does not transfer ownership of the original data.

## Fashion-MNIST

**Dataset:** Fashion-MNIST.  
**Provider:** Zalando SE.  
**Source:** https://github.com/zalandoresearch/fashion-mnist  
**License:** MIT License, Copyright © 2017 Zalando SE.  
**License notice:** [THIRD_PARTY_LICENSES/FASHION-MNIST-LICENSE-MIT.txt](THIRD_PARTY_LICENSES/FASHION-MNIST-LICENSE-MIT.txt)

The experimental archive contains MATLAB-formatted/flattened Fashion-MNIST test-set data and separately computed experimental results. The original image and label data retain Zalando's copyright and the applicable MIT notice.

Suggested publication citation: H. Xiao, K. Rasul, and R. Vollgraf, “Fashion-MNIST: a Novel Image Dataset for Benchmarking Machine Learning Algorithms,” arXiv:1708.07747, 2017.

## Lung RNA-seq data from NCI Genomic Data Commons

**Source:** National Cancer Institute (NCI), Genomic Data Commons (GDC), https://portal.gdc.cancer.gov/  
**Access status of downloaded source files:** Open access.  
**Experimental dataset preparation and processing:** Dimitris Arabadjis, EucApprox project.

The EucApprox lung RNA-seq dataset was assembled and processed by Dimitris Arabadjis from open-access gene-expression and associated metadata files downloaded directly from the GDC. The archive includes the downloaded source files in its `GDC-raw_data/` folder and the assembled/processed MATLAB datasets and experiment results. The original downloaded material and the EucApprox-derived data are distinguishable by their provenance.

The processing code is `EucldnEmbdngPublic/readGDCdata.m`; it associates GDC gene-expression TSV files with GDC file/case metadata and clinical annotations, producing aligned expression matrices and diagnostic labels. The experimental outputs include the transformed input data, distance graphs, embeddings, and performance measurements. The original GDC source data were not created by the EucApprox author.

**GDC source acknowledgment:** A. P. Heath, V. Ferretti, S. Agrawal et al., “The NCI Genomic Data Commons,” *Nature Genetics*, 53, 257–262 (2021). https://doi.org/10.1038/s41588-021-00791-5  
**GDC citation guidance:** https://gdc.cancer.gov/content/how-do-i-cite-nci-gdc

This citation credits the GDC infrastructure and data source; it does **not** attribute the author's dataset selection, assembly, processing, analysis, or EucApprox results to the cited paper. Where source files originate from The Cancer Genome Atlas (TCGA), acknowledgment of the TCGA Research Network also applies: https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/citing.

GDC open-access data may be used with proper accreditation, subject to the GDC's access, sharing, and data-analysis policies. The public experimental archive contains open-access downloads and derived materials; it does not grant permission to access or redistribute GDC controlled-access data.

## Synthetic datasets

Synthetic datasets are generated for EucApprox experiments; some generators and comparative methods are drawn from DRToolbox. Any redistributed code derived from DRToolbox remains subject to DRToolbox's non-commercial permission terms, documented in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Copyright separation

The rights in third-party source datasets are distinct from rights in independently authored EucApprox algorithms, processing scripts, analyses, and original documentation. See [COPYRIGHT.md](COPYRIGHT.md) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
