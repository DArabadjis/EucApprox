# EucApprox

This repository contains the official implementation of the *EucApprox* method for determining the functional form of the optimal Euclidean embedding of data. The provided code implements the methodology as presented in the ongoing preprint ["Euclidean Embedding of Data Using Local Distances"](https://arxiv.org/abs/2605.19243)

Previous versions : -

Present version : 04/2026 - **v0** - Implementation of the NeurIPS 2026 paper 

## Requirements

MATLAB code is developed and tested in MATLAB R2023b.

By default, when iterative sparse Laplacian inversion is enabled (```wRcrsSlv``` is ```true```), the execution is dynamically parallelized thread-wise, using the Parallel Computing Toolbox, if the number of embedding dimensions is large enough to require more than 1sec / embedding iteration. This can be deactivated by setting the last input ```wPrllQ``` of ```EcldnEmbdngDatSgndNewCalc_Edge``` to ```false```.

Exploitation of a multigrid preconditioner for the iterative sparse Laplacian inversion builds upon [*Livermore’s "hypre" library of linear solvers*](https://github.com/hypre-space/hypre). Only precompiled mex files of the preconditioner constructor are provided: ```hypre_amg_mex.mexw64``` (win), ```hypre_amg_mex.mexa64``` (linux). Currently, the multigrid preconditioner incorporation is **experimental** and it may obstruct the algorithm's convergence. By default, this preconditioner is opted out and can be selected by setting the flag ```wAMG``` to ```true``` (the iterative inversion flag ```wRcrsSlv``` should also be ```true```).

Before running comparison experiments, one should:

1. Download Laurens van der Maaten's [Toolbox for Dimensionality Reduction](https://lvdmaaten.github.io/drtoolbox/) and [UMAP](https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)

2. ```
   addpath(genpath(<path to Dimensionality Reduction Toolbox>)); addpath(genpath(<path to UMAP>));
   ```

3. Place the intrinsic versions of the embedding methods, included in the repository's ```IntrnscVrs``` folder, in the ```techniques``` folder within the ```<path to Dimensionality Reduction Toolbox>```

## Code usage pattern

### Getting the embedding

The embedding algorithm needs as input: 

a) the (sparse) adjacency matrix of the distance graph of the data neighborhoods, say ```dstMtrx```,

b) the embedding dimension, say ```dmEstm```.

Then, the $ℝ^{\mathrm{dmEstm}}$ embedding, ```pXQ```, is obtained by the following archetypical call:

```eval
[pXEdge,dstMtrxEdge]=CmptLclBase(dstMtrx);  
wVlmFrm=false;wRcrsSlv=true;wAMG=false;wPlot=true;tlrnc=[];wErrOut=[];wPrllQ=[];  
[pXQ,errs]=EcldnEmbdngDatSgndNewCalc_Edge(dstMtrxEdge,pXEdge(:,1:dmEstm+2),dmEstm,wVlmFrm,wRcrsSlv,wAMG,wPlot,tlrnc,wErrOut,wPrllQ);
```

The flags set before the iterative optimal embedding determination ```EcldnEmbdngDatSgndNewCalc_Edge``` are explained within this function's m-file. The given syntax requires convergence up to machine-epsilon (```tlrnc=[];```) and (due to ```wErrOut=[];```) returns into ```errs``` the pairs $(\mathrm{error~functional} , \mathrm{convergence})$, per iteration.

### Getting the neighborhoods' distance graph

The repository contains functions to estimate the embedding dimensionality ```dmEstm``` and prepare a corresponding distance matrix ```dstMtrx```.

As it is, the code presumes that 

a) the data points are given in a tabular form, say in an ``` observations x features ``` matrix ```pX``` and 

b) that there is a function handler ```fncD``` fitting the pattern ```@(x,y)<distance function>(x,y)```, where ```x``` stands for the data point in hand, ```y``` for its candidate neighbors (e.g the whole data set) and ```<distance function>``` for a function or a formula that computes the distances of all data points of ```y``` from ```x```.

The ```<distance function>``` is supposed to be accurate only locally, within each data point's neighborhood.

The functions ```CmptDstSprs```, ```MnmlNghbrhdsODim``` and ```mkCnnctd```, used to construct the distance graph, exploit the tabular form of the dataset, ```pX```, only to determine the number of observations and iterate on them. If needed, one can easily modify these functions to fit any structure of ```pX```.

#### v0 - pipeline

1. Estimate the data set's intrinsic dimension using the minimal neighborhood method of [*Facco, E., d’Errico, M., Rodriguez, A. et al. Estimating the intrinsic dimension of datasets by a minimal neighborhood information. Sci Rep 7, 12140 (2017)*](https://doi.org/10.1038/s41598-017-11873-y)
   ```eval
   dstMtrx_2NN=CmptDstSprs(pX,fncD,false,3,false);
   dmEstm=round(estmDim2NNs(dstMtrx_2NN,false))+1;
   ```

2. Using ```dmEstm```, two options are provided:  
   a) Determine the minimal distance graph whose neighborhoods are of rank ```dmEstm```
      ```eval
      dstMtrx=CmptDstSprs(pX,fncD,false,dmEstm,false);
      ```
   b) Determine the minimal distance graph whose neighborhoods have estimated rank  as close as possible to ```dmEstm```  
      ```eval
      dstMtrx=MnmlNghbrhdsODim(false,dmEstm-0.5,pX,fncD);
      ```

3.  Connect the distance graph's connected components
      ```eval
      dstMtrx=mkCnnctd(pX,dstMtrx,fncD);
      ```


## Evaluation

### Experiments and demo code

```batchSynthetic.m``` - For experiments on drtoolbox synthetic datasets

```demo_RNAseq_embddng.m``` - For experiments on RNAseq dataset  
   - CmptDstSprs is slow here due to the 60660 dims of ```pXL```. Better preload ```PBMC_5k_DistEuc_FullDim.mat```
         
```batchRealMNIST.m``` - k-fold experiments on MNIST and FMNIST test sets 

```prcss_MNIST_FMNIST.m``` - embedding experiment on the total MNIST and FMNIST test sets 

### Experiments' data

Data are provided in  
https://e.pcloud.link/publink/show?code=kZLnliZ51ghbjdeRGRfIYEh7dVBAywNz6Ly

#### DNA-seq

Dataset collected from https://portal.gdc.cancer.gov/ and its raw data folder is in ```GDC-raw_data```

The extracted features are in ```GDC-raw_data/GDC-lung-RNAseq/GDC-lung-RNAseq-tpm.mat```

The extraction function is ```readGDCdata.m```

**Results**:
```PBMC_5k_DistEuc_FullDim.mat```: The results on the lung cancer dataset.
- ```pX``` - original RNAseq data
   
- ```pY``` - cases' classification
   
- ```Diagnwseis``` - diagnosis of each class
   
- ```pXL``` - log2(1+pX)
   
- ```fncD ``` - distance function
   
- ```dmEstm_dmFull``` - intrinsic data dimension
   
- ```dstMtrxEdgeSparse``` - the graph adjacency matrix
   
- ```pXEdgeSparse``` - the metric trivializing local transports
   
- ```pX_Euc_short_short_sparse_dmEstm``` - the Euclidean embedding
   
- ```pX_Euc_VlmFrm_short_short_sparse_dmEstm``` - the Euclidean embedding computed after volume normalization of the local metric

- ```datOut_dmEstm_*``` - <*> method's classification capability  
   - ```prbsAcc``` - accuracy  
	- ```prbs_ARI``` - adjusted rand index  
	- ```prbs_NMI``` - normalized mutual information

- ```datOutQuality_dmEstm_*``` - <*> method's embedding quality

**Visualization** via ```demo_RNAseq_visualize.m```

#### MNIST & FMNIST test set
```(f)mnist_test-10fold-short-short-sparse-dmFull-NewCalc.mat```: 10-fold results - The test set is randomly split into 10 independent parts to evaluate how methods perform in this high-dim space if the full distribution of 10k images is sparsely sampled.
