# EucApprox — research data attribution and usage notices

This notice accompanies research datasets and computed results made available separately from the EucApprox source-code repository: https://github.com/DArabadjis/EucApprox.

**The included data are not all owned or licensed by the EucApprox author.** Credit and comply with each originating data provider's applicable terms. See the repository's [`DATASETS.md`](https://github.com/DArabadjis/EucApprox/blob/main/DATASETS.md) for provenance and processing details.

## MNIST

Original handwritten-digit data: Yann LeCun and Corinna Cortes (derived from NIST data). License described by Keras as [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/). Original images/labels may be supplied here in MATLAB-tabular form (flattened pixels); that format transformation and derived computations are not a claim of original ownership. See https://keras.io/api/datasets/mnist/.

## Fashion-MNIST

Copyright (c) 2017 Zalando SE. Distributed under the [MIT License](https://github.com/zalandoresearch/fashion-mnist/blob/master/LICENSE); when reproducing the dataset, include its original copyright and full license notice. The source-code repository includes a copy: [`THIRD_PARTY_LICENSES/FASHION-MNIST-LICENSE-MIT.txt`](https://github.com/DArabadjis/EucApprox/blob/main/THIRD_PARTY_LICENSES/FASHION-MNIST-LICENSE-MIT.txt). Converted files and experimental results are not a change of copyright ownership of the original images.

## Author-prepared lung RNA-seq dataset from NCI GDC open-access files

This lung RNA-seq experimental dataset was **assembled and processed by Dimitris Arabadjis for EucApprox** using open-access files downloaded directly from the [NCI Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/). The original downloaded source files are included in the archive's `GDC-raw_data/` folder, alongside processed MATLAB datasets and experimental outputs. The MATLAB processing code is provided in the EucApprox repository as `EucldnEmbdngPublic/readGDCdata.m`.

**Source acknowledgment, not attribution of downstream processing:** A. P. Heath, V. Ferretti, S. Agrawal, et al., *The NCI Genomic Data Commons*, *Nature Genetics* **53**, 257–262 (2021), https://doi.org/10.1038/s41588-021-00791-5. GDC asks researchers using its data to cite this paper. The author performed the dataset assembly, transformations and EucApprox evaluations. If the GDC downloads include TCGA source projects, also credit the [TCGA Research Network](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/citing); follow attribution policies of any other originating projects.

GDC open-access materials may be analyzed and used in publications with proper accreditation. Users must not attempt to reidentify participants; see the [GDC policies](https://gdc.cancer.gov/analyze-data/data-analysis-policies). This archive is intended to contain only the open-access source files selected by the dataset preparer and their derived results, not controlled-access records.

## Algorithmic outputs

Graph distances, embeddings, classification indices, and numeric experiment results produced by EucApprox are computational outputs. Their distribution does not change the rights in underlying third-party image, RNA-seq, or clinical data they may contain.
