# Datasets: provenance, attribution, and usage terms

This file documents datasets used in the experiments described by the EucApprox repository and its [NeurIPS 2026 paper](https://arxiv.org/abs/2605.19243). **Dataset licenses and access conditions remain distinct from the copyright status of EucApprox software.** The repository distributes scripts; a separate [research-data archive](https://e.pcloud.link/publink/show?code=kZLnliZ51ghbjdeRGRfIYEh7dVBAywNz6Ly) is linked from its README.

## MNIST handwritten digits

- **Creators / attribution:** Yann LeCun and Corinna Cortes; derivative of original NIST data.
- **Source:** https://yann.lecun.com/exdb/mnist/ (historical dataset homepage); see also https://keras.io/api/datasets/mnist/.
- **License statement:** [Creative Commons Attribution–ShareAlike 3.0](https://creativecommons.org/licenses/by-sa/3.0/) (as documented by Keras).
- **Data used here:** The 10,000-image test split, converted to MATLAB/tabular representation with flattened image pixels, plus labels; see `(f)mnist.mat` / experiment files in the external data archive.
- **Changes:** File-format conversion / flattening into MATLAB arrays and subsequent derived graph/embedding/evaluation outputs. The MNIST source data and labels were not created by EucApprox.
- **Attribution when sharing:** Credit MNIST's copyright holders, identify the CC BY-SA 3.0 terms, link the dataset/source and license, and indicate the conversions. Further adaptations may trigger the license's ShareAlike requirements. These dataset terms do not govern independent EucApprox source files simply because they are stored in the same research collection.

## Fashion-MNIST

- **Provider / copyright:** Zalando SE (2017).
- **Source:** https://github.com/zalandoresearch/fashion-mnist
- **License notice:** [MIT license](https://github.com/zalandoresearch/fashion-mnist/blob/master/LICENSE); a copy is in [`THIRD_PARTY_LICENSES/FASHION-MNIST-LICENSE-MIT.txt`](THIRD_PARTY_LICENSES/FASHION-MNIST-LICENSE-MIT.txt).
- **Data used here:** The 10,000-image test split, flattened/converted into MATLAB arrays with class labels, plus graph/embedding/evaluation outputs.
- **Changes:** MATLAB file-format conversion and the additional derived experimental constructs/results.
- **Redistribution:** Where Fashion-MNIST data are redistributed, accompany them with the upstream copyright and MIT permission notice; do not represent the data as exclusively licensed by EucApprox.
- **Suggested academic citation:** Han Xiao, Kashif Rasul, and Roland Vollgraf, *Fashion-MNIST: a Novel Image Dataset for Benchmarking Machine Learning Algorithms*, arXiv:1708.07747 (2017).

## NCI Genomic Data Commons (GDC) — lung RNA-seq data

- **Data provider:** National Cancer Institute (NCI) Genomic Data Commons: https://portal.gdc.cancer.gov/
- **Access policy:** https://gdc.cancer.gov/analyze-data/data-analysis-policies
- **Citation:** Heath, A. P., Ferretti, V., Agrawal, S., et al. *The NCI Genomic Data Commons*. **Nature Genetics** 53, 257–262 (2021). https://doi.org/10.1038/s41588-021-00791-5. Follow any additional attribution instructions for the originating project(s).
- **Relevant local processing code:** `EucldnEmbdngPublic/readGDCdata.m` combines downloaded expression TSV files with GDC file-to-case metadata and clinical information; `EucldnEmbdngPublic/demo_RNAseq_embddng.m` is the demonstration script.
- **Processed data described in the README:** `GDC-raw_data/GDC-lung-RNAseq/GDC-lung-RNAseq-tpm.mat` and `PBMC_5k_DistEuc_FullDim.mat`; the latter stores expression matrices, diagnosis/classification information, log-transformed values (`log2(1+pX)`), graph data, embeddings, and evaluation outputs.

**Access/redistribution qualification:** GDC's policy distinguishes *open access* from *controlled access*. Its open-access material has no GDC restrictions on analysis or publication other than its prohibition on reidentification attempts; controlled-access data may have project-specific Data Use Agreements (DUAs). **The GitHub repository zip does not contain the source download manifest or file IDs, so it is not possible to confirm from the examined material that every GDC-derived record in the external archive was originally open access.** Do not treat this description as verification that processed clinical data are redistributable. Before distributing/continuing to distribute the external RNA-seq archive, check the access classification and applicable project terms for every contributing file and associated clinical field.

**Provenance details to complete from the original GDC download (not inferable from the code archive):**

- GDC project identifier(s): **TO BE VERIFIED / FILLED IN BY MAINTAINER**.
- Manifest / GDC file UUIDs and GDC data release date: **TO BE VERIFIED / FILLED IN BY MAINTAINER**.
- Each contributing expression and clinical file's access category (`open`/`controlled`): **TO BE VERIFIED / FILLED IN BY MAINTAINER**.
- Expression workflow and precise measurement column supplied as `varRd` to `readGDCdata.m` (e.g., TPM): **TO BE VERIFIED / FILLED IN BY MAINTAINER**.
- Inclusion/exclusion criteria and diagnostic-label mapping: **TO BE VERIFIED / FILLED IN BY MAINTAINER**.

**Privacy restriction:** Do not attempt reidentification of any GDC participant. Any controlled-access content must not be republished without appropriate authorization.

## Synthetic test manifolds and DRToolbox examples

Synthetic experiments are generated by the authors' code and, where indicated, the externally obtained DRToolbox example generators. Follow DRToolbox's non-commercial software restrictions for any redistributed modified generator code; see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

## Separate data archive

When distributing MATLAB copies of datasets outside GitHub, place [`DATA_ARCHIVE_NOTICE.md`](DATA_ARCHIVE_NOTICE.md) (and the applicable license texts/links) **inside that download** so attribution does not depend on readers first visiting the source-code repository.
