# Third-party software notices

The independently authored EucApprox implementation is not distributed under a general open-source license. Copyright and permission information for that code appears in [LICENSE.md](LICENSE.md). The third-party components identified below retain their respective copyrights and licensing conditions.

## Modified DRToolbox routines

**Original software:** Matlab Toolbox for Dimensionality Reduction (DRToolbox), Laurens van der Maaten, Delft University of Technology.  
**Source and permission terms:** https://lvdmaaten.github.io/drtoolbox/

The following files are adapted from original DRToolbox routines:

- `IntrnscVrs/hlle_intrnsc.m`
- `IntrnscVrs/laplacian_eigen_intrnsc.m`
- `IntrnscVrs/lle_intrnsc.m`
- `IntrnscVrs/ltsa_intrnsc.m`

They implement intrinsic-distance-graph variants of the corresponding dimensionality-reduction algorithms. The original DRToolbox attribution and permission statements are retained in the source files.

DRToolbox permits use, modification, and redistribution **for non-commercial purposes only**. These conditions continue to govern the adapted DRToolbox code and are not superseded by the copyright notice covering independently authored EucApprox code.

## HYPRE

**Upstream project:** HYPRE Project, Lawrence Livermore National Laboratory and other HYPRE contributors.  
**Project:** https://github.com/hypre-space/hypre  
**Licensing:** Apache License 2.0 **or** MIT License, at the recipient's option.  
**MIT license notice:** [THIRD_PARTY_LICENSES/HYPRE-LICENSE-MIT.txt](THIRD_PARTY_LICENSES/HYPRE-LICENSE-MIT.txt).

The repository includes the following optional, precompiled MATLAB interfaces associated with the HYPRE algebraic-multigrid preconditioner:

- `EucldnEmbdngPublic/hypre_amg_mex.mexa64`
- `EucldnEmbdngPublic/hypre_amg_mex.mexw64`

The HYPRE copyright and license conditions apply to the HYPRE portions of those binaries. The HYPRE license does not automatically extend to independently authored interface code or other separately licensed components. Upstream HYPRE notices are available in its [COPYRIGHT](https://github.com/hypre-space/hypre/blob/master/COPYRIGHT) and [NOTICE](https://github.com/hypre-space/hypre/blob/master/NOTICE) files. The default EucApprox configuration does not enable this optional preconditioner.

## MATLAB UMAP

**Implementation:** Connor Meehan, Jonathan Ebrahimian, Wayne Moore, and Stephen Meehan, *Uniform Manifold Approximation and Projection (UMAP)*, MATLAB Central File Exchange.  
**Source and license:** https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap

The MATLAB UMAP implementation is an external comparison dependency and is **not included** in this repository. Its distribution and use are governed by the license accompanying the separately obtained implementation.

## MATLAB

MATLAB and MathWorks toolboxes are separate products governed by their applicable MathWorks licenses. They are not redistributed in this repository.

## Experimental datasets

MNIST, Fashion-MNIST, and NCI Genomic Data Commons data are documented in [DATASETS.md](DATASETS.md). The original data providers' terms are separate from those governing EucApprox software.
