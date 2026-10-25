# Third-party software notices

This document describes third-party components in, or used by, the **EucApprox** research repository. It does **not** license the independently authored EucApprox code. The absence of a general EucApprox `LICENSE` file means that no general permission to reproduce, modify, or redistribute that original code is granted by this repository; applicable copyright law, exceptions, and GitHub platform terms continue to apply. Third-party components retain their own permissions and restrictions.

## 1. Modified DRToolbox routines — **included source files**

**Upstream:** Laurens van der Maaten, *Matlab Toolbox for Dimensionality Reduction*  
**Project:** https://lvdmaaten.github.io/drtoolbox/  
**Permission statement:** The upstream project permits use, modification and redistribution **for non-commercial purposes only**. See its website and the copyright/permission headers in the files listed below.

Files in this repository derived from DRToolbox:

- `IntrnscVrs/hlle_intrnsc.m`
- `IntrnscVrs/laplacian_eigen_intrnsc.m`
- `IntrnscVrs/lle_intrnsc.m`
- `IntrnscVrs/ltsa_intrnsc.m`

**Modification disclosure:** These are modified, *intrinsic-distance-graph* variants of DRToolbox's corresponding manifold-learning algorithms. The routines have been adapted to operate on graph distances/neighborhood information rather than requiring the original feature-based inputs; additional changes appear in each source file. Original authorship and permission notices are retained in the source headers.

**Important:** The upstream non-commercial restriction continues to apply to these derived files. Public distribution of EucApprox does not expand those permissions, nor does the absence of an EucApprox license remove the upstream permissions.

## 2. HYPRE sparse solvers — **included binary files**

**Upstream:** HYPRE Project, Lawrence Livermore National Laboratory and other contributors  
**Project:** https://github.com/hypre-space/hypre  
**Upstream licensing:** `Apache-2.0 OR MIT` (recipient may choose either).  
**MIT notice:** [`THIRD_PARTY_LICENSES/HYPRE-LICENSE-MIT.txt`](THIRD_PARTY_LICENSES/HYPRE-LICENSE-MIT.txt)  
**Upstream notices:** [COPYRIGHT](https://github.com/hypre-space/hypre/blob/master/COPYRIGHT), [NOTICE](https://github.com/hypre-space/hypre/blob/master/NOTICE), [LICENSE-APACHE](https://github.com/hypre-space/hypre/blob/master/LICENSE-APACHE).

EucApprox contains optional, precompiled MATLAB MEX interfaces associated with HYPRE's algebraic-multigrid preconditioner:

- `EucldnEmbdngPublic/hypre_amg_mex.mexa64`
- `EucldnEmbdngPublic/hypre_amg_mex.mexw64`

They are **optional** and the repository's default setting leaves AMG disabled. The MIT license text is reproduced for the upstream HYPRE component; this is **not** a grant of an MIT license for EucApprox code.

**Binary provenance verification still required:** The repository archive does not identify the build date/version, author and license of the MATLAB MEX *wrapper*, or all libraries statically/dynamically incorporated into these two binaries. Before redistributing binaries, the maintainer must verify that HYPRE's licensing statement applies to the compiled version and that the wrapper and any bundled dependency licenses permit redistribution. The notice above does not independently establish those rights. If the provenance cannot be verified, remove the binaries from the public release and provide separate installation/build instructions instead.

## 3. MATLAB UMAP — **external download, not bundled**

**Distribution:** Connor Meehan, Jonathan Ebrahimian, Wayne Moore, and Stephen Meehan, *Uniform Manifold Approximation and Projection (UMAP)*, MATLAB Central File Exchange.  
**Project/download:** https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap  
**Original UMAP algorithm:** Leland McInnes, John Healy and James Melville.

UMAP is a separately installed comparison dependency; its source and binaries are **not included** in this repository archive. Users must obtain it from its original distributor and comply with the terms accompanying the particular version they download. Do not assume that the original Python UMAP project's license necessarily governs every file in this MATLAB distribution.

## 4. MATLAB and other external requirements

MATLAB and any required MathWorks toolboxes are proprietary third-party products, separately licensed by MathWorks. EucApprox does not distribute MATLAB. Other baseline tools installed by users retain their respective source licenses.

## 5. Datasets are separately documented

MNIST, Fashion-MNIST, and NCI Genomic Data Commons (GDC) data are **not licensed by the EucApprox copyright holder**. Their attribution, reuse conditions, download sources, and processed-data status are documented in [`DATASETS.md`](DATASETS.md). Some processed files are distributed through a **separate external data archive**, not in the GitHub source repository. That archive should contain its own copy of [`DATA_ARCHIVE_NOTICE.md`](DATA_ARCHIVE_NOTICE.md) and any accompanying required notices.

## Questions or corrections

Licensing statements are based on the upstream project notices linked above and the contents of the archived public repository. If another upstream license or copyright attribution applies to a specific file, that file's actual license takes precedence. Direct questions about permission to reuse the independently authored EucApprox implementation to the repository maintainer.
