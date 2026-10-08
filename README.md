# EucApprox

The official implementation of the paper ["Euclidean Embedding of Data Using Local Distances"](https://arxiv.org/abs/2605.19243). 

>📋  Optional: include a graphic explaining your approach/main result, bibtex entry, link to demos, blog posts and tutorials

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


## Evaluation

To evaluate my model on ImageNet, run:

```eval
python eval.py --model-file mymodel.pth --benchmark imagenet
```

>📋  Describe how to evaluate the trained models on benchmarks reported in the paper, give commands that produce the results (section below).

## Pre-trained Models

You can download pretrained models here:

- [My awesome model](https://drive.google.com/mymodel.pth) trained on ImageNet using parameters x,y,z. 

>📋  Give a link to where/how the pretrained models can be downloaded and how they were trained (if applicable).  Alternatively you can have an additional column in your results table with a link to the models.

## Results

Our model achieves the following performance on :

### [Image Classification on ImageNet](https://paperswithcode.com/sota/image-classification-on-imagenet)

| Model name         | Top 1 Accuracy  | Top 5 Accuracy |
| ------------------ |---------------- | -------------- |
| My awesome model   |     85%         |      95%       |

>📋  Include a table of results from your paper, and link back to the leaderboard for clarity and context. If your main result is a figure, include that figure and link to the command or notebook to reproduce it. 


## Contributing

>📋  Pick a licence and describe how to contribute to your code repository. 
