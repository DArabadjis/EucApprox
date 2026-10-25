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

## Author-prepared lung RNA-seq dataset derived from NCI GDC open-access files

- **Dataset preparation:** Dimitris Arabadjis (EucApprox project). The lung RNA-seq experimental dataset was assembled and processed by the author from **open-access files that the author downloaded directly** from the [National Cancer Institute Genomic Data Commons (GDC)](https://portal.gdc.cancer.gov/). It is **not** a dataset published or processed by Heath et al.
- **Source data:** The original downloaded GDC files are provided alongside the experimental results, under `experiments_results/GDC-raw_data/` in the separately distributed [experiments/data archive](https://e.pcloud.link/publink/show?code=kZLnliZ51ghbjdeRGRfIYEh7dVBAywNz6Ly). `GDC-raw_data` is the archive's folder name; it does not imply that the files are raw sequencing reads, rather than GDC-provided gene-expression quantifications or other source files.
- **Author's processing code:** `EucldnEmbdngPublic/readGDCdata.m` reads GDC expression TSV files, associates them with downloaded JSON file/case metadata and clinical TSV information, and constructs the aligned expression matrix, gene metadata and diagnostic/classification labels. `EucldnEmbdngPublic/demo_RNAseq_embddng.m` demonstrates the embedding experiment.
- **Prepared dataset:** `GDC-raw_data/GDC-lung-RNAseq/GDC-lung-RNAseq-tpm.mat` contains the assembled expression data, while `PBMC_5k_DistEuc_FullDim.mat` contains prepared inputs and experiment-derived distances, graphs, embeddings and evaluations as described in the repository README. The original downloaded files and the author's processed `.mat` outputs have distinct provenance.
- **GDC acknowledgment (not authorship of the processed dataset):** A. P. Heath, V. Ferretti, S. Agrawal, et al., *The NCI Genomic Data Commons*, **Nature Genetics** 53, 257–262 (2021), https://doi.org/10.1038/s41588-021-00791-5. The GDC [explicitly requests this reference](https://gdc.cancer.gov/content/how-do-i-cite-nci-gdc) when researchers use data from its portal. Its citation acknowledges the GDC infrastructure and does not attribute the author's downstream curation, transformations, analyses, or embeddings to the paper's authors.
- **Project-specific acknowledgment:** If the downloaded source files belong to TCGA, include the [TCGA Research Network acknowledgment](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/citing); for other contributing GDC projects, follow their attribution policies where applicable. Specific project IDs should be recorded from the original download metadata rather than guessed.
- **Use conditions:** The dataset was assembled from **GDC open-access source files**, as confirmed by the dataset preparer. GDC permits analysis and publication of open-access data with proper accreditation, but [prohibits reidentification attempts](https://gdc.cancer.gov/analyze-data/data-analysis-policies). This statement does not authorize redistribution of any controlled-access files; only the open-access source files and properly derived material should be provided in the public archive.

**Suggested dataset attribution:** "Lung RNA-seq experimental dataset compiled and processed by Dimitris Arabadjis for EucApprox from publicly available NCI Genomic Data Commons files (original downloaded files and processed outputs supplied in the EucApprox experiments archive)." Cite the GDC infrastructure separately as above, and the original project where applicable.

## Synthetic test manifolds and DRToolbox examples

Synthetic experiments are generated by the authors' code and, where indicated, the externally obtained DRToolbox example generators. Follow DRToolbox's non-commercial software restrictions for any redistributed modified generator code; see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

## Separate data archive

When distributing MATLAB copies of datasets outside GitHub, place [`DATA_ARCHIVE_NOTICE.md`](DATA_ARCHIVE_NOTICE.md) (and the applicable license texts/links) **inside that download** so attribution does not depend on readers first visiting the source-code repository.
