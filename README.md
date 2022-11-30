# rMultiNet: An R Package For Multilayer Networks Analysis

## Overview
rMultiNet introduces an extension R package that includes a variety of traditional and state-of-the-art tensor decomposition methods for mixture multilayer networks analysis.The package is developed with the modular pipeline mode: generative modeling, embedding algorithms and visualization.  rMultiNet aims to help study complex networks, especially mixture multilayer networks.


![](https://github.com/ChenyuzZZ73/rMultiNet/blob/main/framework.png)

## Getting Started

### Prerequisites
Before using the rMultiNet, you need to install the following packages.
```
require(rTensor)
require(ggplot)
require(Matrix)
require(plotly)
```

### Installation

Install the stable version of R-package from Git rMultiNet package directly with:
```sh
devtools::install_github("ChenyuzZZ73/rMultiNet")
bulid()
library(rMultiNet)
```

### Usage example
The package is developed with the modular pipeline mode: generative modeling, embedding algorithms and visualization. 
#### Generative modeling
```sh
library(rMultiNet)
GenerateMMSBM(n, m, L, K, d = NULL, r = NULL)
GenerateMMLSM(n, m, L, rank, U mean= 0.5, cmax =1, d, int type = ‘Uniform’, kernel fun = ‘logit’, scale par=1)
```
Also, rMultiNet provides three datasets for study：human malaria parasite gene network, worldwide food trading network and UN Commodity trading network.
```sh
load("~/Desktop/rMultiNet/data/malariagene/malaria.RData")
```

#### Emdedding algorithms
Take one of the algorithms for example.
```sh
InitializationMMSBM(tnsr, ranks=NULL)
PowerIteration(tnsr, ranks=NULL, type=”TWIST”, U 0 list, delta1=1000, delta2=1000, max iter = 25, tol = 1e-05)
```

#### Visulization
```sh
Embedding_network(network membership,L, paxis=2)
Community_cluster_km(embedding,type,cluster number)
Community_cluster_dbscan(embedding,type,eps value =.05, pts_value=5)
```




## Citation
If you use `rMultiNet` or reference our tutorials in a presentation or publication, we would appreciate citations of our library.
>



## Reference
- Bing-Yi Jing, Ting Li, Zhongyuan Lyu, and Dong Xia. Community detection on mixture
multilayer networks via regularized tensor decomposition. The Annals of Statistics, 49
(6):3181–3205, 20
- Zhongyuan Lyu, Dong Xia, and Yuan Zhang. Latent space model for higher-order networks and generalized tensor decomposition. arXiv preprint arXiv:2106.16042, 2021.
- Xiaowen Dong, Pascal Frossard, Pierre Vandergheynst, and Nikolai Nefedov. Clustering
with multi-layer graphs: A spectral perspective. IEEE Transactions on Signal Processing,
60(11):5820–5831, 2012
