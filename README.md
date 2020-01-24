HAPOD - Hierarchical Approximate Proper Orthogonal Decomposition
================================================================

* HAPOD - Hierarchical Approximate POD
* version: 3.0 (2020-01-24)
* by: C. Himpe (0000-0003-2194-6754), S. (0000-0003-0439-7212)
* under: BSD 2-Clause License (opensource.org/licenses/BSD-2-Clause)
* summary: Fast distributed or incremental POD computation.

## Scope

* Proper Orthogonal Decomposition (POD)
* Singular Value Decomposition (SVD)
* Principal Component Analysis (PCA)
* Empirical Orthogonal Functions (EOF)
* Karhunen-Loeve Transformation (KLT)
* Empirical Eigenfunctions (EEF)

## Applications

* Dimension Reduction
* Model Reduction
* Low-Rank Approximation
* Data Compression

## Features

* Error-driven
* Single pass (each data vector is needed only once)
* Column-wise data partition
* Custom SVD backends

## Functionality

* Standard POD
* Incremental HAPOD -> for memory-limited environments: single-board-computers
* Distributed HAPOD -> for distributd memory environments: super-computers
* Distributed-of-Incremental HAPOD

## Algorithm

C. Himpe, T. Leibner and S. Rave.
"[Hierarchical Approximate Proper Orthogonal Decomposition](http://hdl.handle.net/21.11116/0000-0002-5342-6)".
SIAM Journal on Scientific Computing, 40(5): A3267--A3292, 2018.

## Compatibility

* GNU Octave >= 4.2
* Mathworks MATLAB >= 2013b

## Basic Usage

```
[svec,sval,meta] = hapod(data,bound,topo,relax,meta,depth,mysvd)
```

## Arguments

* `data`   {cell}  - snapshot data set, partitioned by column (blocks)
* `bound` {scalar} - mean L_2 projection error bound
* `topo`  {string} - tree topology (see **Topology**)
* `relax` {scalar} - relaxation parameter in (0,1] (see **Relaxation**)
* `depth` {scalar} - total number of levels in tree (only required for `incr_1`)
* `meta`  {struct} - meta information structure (see **Meta-Information**)
* `mysvd` {handle} - custom SVD backend (see **Custom SVD**) 

## Return Values

* `svec` {matrix} POD modes (column vectors)
* `sval` {vector} Singular values (column vector)
* `meta` {struct} Meta-information structure

## Topology

The HAPOD is computed based on a tree topology `topo`, with the data partitions
at the tree's leafs. The following topologies are available:

* `'none'`   Standard POD
* `'incr'`   Incremental HAPOD (Complete)
* `'incr_1'` Incremental HAPOD (Child nodes)
* `'incr_r'` Incremental HAPOD (Root node)
* `'dist'`   Distributed HAPOD (Complete)
* `'dist_1'` Distributed HAPOD (Child nodes)
* `'dist_r'` Distributed HAPOD (Root node)

If all data partitions can be passed as the data argument, the types: `none` 
(standard POD), `incr`(emental) HAPOD or `dist`(ributed) HAPOD are applicable.
In case only a single partition is passed, the types: `incr_1` and `dist_1`
should be used for the child nodes of the associated HAPOD tree, while the
types: `incr_r` and `dist_r` should be used for the root nodes. The returned
meta-information structure (or a cell-array thereof) has to be passed to the
parent node in the associated HAPOD tree. 

## Relaxation

The relaxation parameter `w` (`0 < w < 1`) balances accuracy versus speed.
Larger `w` near one means be more accurate, while `w` near zero means faster
computation. The default value is `w = 0.5`.

## Meta-Information

The `meta` structure contains the following meta-information of the completed
sub-tree:

* `nSnapshots` - Number of data columns passed to this hapod and its children.
* `nModes`     - Number of intermediate modes
* `tNode`      - Computational time at this hapod's branch

The argument `meta` only needs to be passed for topology argument `incr_r`,
`dist_r` and `incr_1` unless it is first leaf. This means especially the user
never has to create such a structure, as if it is required it is given as a
previous return value.

## Custom SVD

Via the `mysvd` argument a custom SVD function can be provided via a function
handle with the following signature:

```
[U,d] = mysvd(X)
```

for a data matrix `X`, and return left singular vectors in matrix `U` and
singular values in column vector `d`. By default (or `mysvd` = `eco`) a standard
rank-revealing SVD is used. Additionally, by `mysvd` = `mos` the method of
snapshots can be selected.

## Getting Started

Run the sample code:

```
RUNME()
```

which demonstrates the different implemented HAPOD variants and can be used
as a template.

## Cite As

C. Himpe, T. Leibner and S. Rave.
"[Hierarchical Approximate Proper Orthogonal Decomposition](https://doi.org/10.1137/16M1085413)".
SIAM Journal on Scientific Computing, 40(5): A3267--A3292, 2018.

## Used In

* P. Benner, C. Himpe.
"[Cross-Gramian-Based Dominant Subspaces](https://doi.org/10.1007/s10444-019-09724-7)".
Advances in Computational Mathematics, 45(5): 2533--2553, 2019.

* C. Himpe, T. Leibner, S. Rave, J. Saak.
"[Fast Low-Rank Empirical Cross Gramians](https://doi.org/10.1002/pamm.201710388)".
Proceedings in Applied Mathematics and Mechanics, 17: 841--842, 2017.

* C. Himpe, T. Leibner, S. Rave.
"[Comprehensive Memory-Bound Simulations on Single Board Computers](https://doi.org/10.5281/zenodo.814497)".
Extended Abstract, 2nd Conference on Power Aware Computing (PACO), 2017.

## See Also

* C. Himpe, T. Leibner, S. Rave.
"[HAPOD - Fast, Simple and Reliable Distributed POD Computation](https://doi.org/10.11128/arep.55.a55283)".
ARGESIM Report 55 (MATHMOD 2018 Volume): 119--120, 2018.

* C. Himpe and S. Rave.
"[HAPOD - Hierarchical Approximate POD](https://himpe.science/poster/rave16_morml.pdf)".
Data-Driven Model Order Reduction and Machine Learning (MORML), 2016.
