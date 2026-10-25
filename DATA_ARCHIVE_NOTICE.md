# EucApprox — research data attribution and usage notices

This notice accompanies research datasets and computed results made available separately from the EucApprox source-code repository: https://github.com/DArabadjis/EucApprox.

**The included data are not all owned or licensed by the EucApprox author.** Credit and comply with each originating data provider's applicable terms. See the repository's [`DATASETS.md`](https://github.com/DArabadjis/EucApprox/blob/main/DATASETS.md) for provenance and processing details.

## MNIST

Original handwritten-digit data: Yann LeCun and Corinna Cortes (derived from NIST data). License described by Keras as [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/). Original images/labels may be supplied here in MATLAB-tabular form (flattened pixels); that format transformation and derived computations are not a claim of original ownership. See https://keras.io/api/datasets/mnist/.

## Fashion-MNIST

Copyright (c) 2017 Zalando SE. Distributed under the [MIT License](https://github.com/zalandoresearch/fashion-mnist/blob/master/LICENSE); when reproducing the dataset, include its original copyright and full license notice. The source-code repository includes a copy: [`THIRD_PARTY_LICENSES/FASHION-MNIST-LICENSE-MIT.txt`](https://github.com/DArabadjis/EucApprox/blob/main/THIRD_PARTY_LICENSES/FASHION-MNIST-LICENSE-MIT.txt). Converted files and experimental results are not a change of copyright ownership of the original images.

## NCI GDC lung RNA-seq

Source: National Cancer Institute Genomic Data Commons, https://portal.gdc.cancer.gov/. Cite: Heath et al., *The NCI Genomic Data Commons*, Nat. Genet. **53**, 257–262 (2021), https://doi.org/10.1038/s41588-021-00791-5.

GDC prohibits reidentification attempts and distinguishes open- from controlled-access data (https://gdc.cancer.gov/analyze-data/data-analysis-policies). The publisher of this archive must verify access permissions for its GDC-derived clinical and expression records; the public code zip alone does **not** establish this status. Do **not** redistribute controlled-access data or data subject to prohibitive project DUAs. Add the project's identifiers and file manifest to `DATASETS.md` before treating GDC data provenance as complete.

## Algorithmic outputs

Graph distances, embeddings, classification indices, and numeric experiment results produced by EucApprox are computational outputs. Their distribution does not change the rights in underlying third-party image, RNA-seq, or clinical data they may contain.
