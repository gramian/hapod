HAPOD - Hierarchical Approximate Proper Orthogonal Decomposition
================================================================

* HAPOD - Hierarchical Approximate POD
* version: 2.0 ( 2019-01-21 )
* by: Christian Himpe (0000-0003-2194-6754), Stephan Rave (0000-0003-0439-7212)
* under: BSD 2-Clause License (open-source)
* summary: Distributed or incremental POD / SVD computation

## Scope

* Proper Orthogonal Decomposition (POD)
* (Truncated) Singular Value Decomposition (SVD)
* (Sparse) Principal Compoenent Analysis (PCA)
* Empirical Orthogonal Functions (EOF)
* Empirical Eigenfunctions
* Karhunen-Loeve Decomposition
* Dimension Reduction
* Model Reduction | Model Order Redction

## Features

* Standard POD
* Incremental (HA)POD
* Distributed (HA)POD
* Custom SVD backend
* Method of Snapshots

## Compatibility

* GNU Octave >= 4.2
* Mathworks MATLAB >= 2013b

## Basic Usage

```
[svec,sval,snfo] = hapod(data,bound,type,relax,config,mysvd)
```

## Documentation

### Arguments

* `data` {cell array}
* `bound` {scalar} L2 mean projection error bound
* `type` {string} HAPOD tree type
    * `'incr'` Incremental HAPOD (Complete)
    * `'incr_1'` Incremental HAPOD (Non-root node)
    * `'incr_r'` Incremental HAPOD (Root node)
    * `'dist'` Distributed HAPOD (Complete)
    * `'dist_1'` Distributed HAPOD (Non-root node)
    * `'dist_r'` Distributed HAPOD (Root node)
    * `'none'` Standard POD
* `relax` {scalar} Relaxation parameter in (0,1] by default: 0.5.
* `config` {structure} Configuration structure by default empty.
    * `nLevels` - Total number of levels in tree
    * `nSnapshots` - Number of data snapshots (columns) seen at each node
    * `nModes` - Number of modes resulting at each node
    * `tNode` - Time consumed at each node
* `mysvd` {handle} Function handle to a custom POD method.
    * by default an economic SVD is used.
    * alternatively the method-of-snapshots can be used via `'mos'`.
    * otherwise a function handle can be provided.

### Return Values

* `svec` {matrix} POD modes (column vectors)
* `sval` {vector} Singular values
* `snfo` {structure} Information structure

### Usage

If all data partitions can be passed as the data argument, the types: `none` 
(standard POD), `incr`(emental) HAPOD or `dist`(ributed) HAPOD are applicable.
In case only a single partition can be passed, the types: `incr_1` and `dist_1`
should be used for the non-root nodes of the associated HAPOD tree, while the
types: `incr_r` and `dist_r` should be used for the root node. The returned
information structure (or a cell-array thereof) can be passed to the parent
nodes in the associated HAPOD tree. 

### Information and Configuration Structure

* `nLevels` - Total number of levels in tree
* `nSnapshots` - Number of data columns passed to this hapod and its children.
* `nModes` - Number of intermediate modes
* `tNode` - computational time at this hapod's branch

Only for `incr_1`, the number of levels `nLevels` in the tree needs to be
provided in a structure. The other fields are filled by the `hapod` function.

### Custom SVD Backend

Signature, arguments and return values as for MATLAB's svd function.

```
[U,D,V] = mysvd(X)
```

### Getting Started

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
"[Cross-Gramian-Based Dominant Subspaces](https://arxiv.org/abs/1809.08066)".
arXiv, math.OC: 1809.08066, 2018.

* T. Taddei, A.T. Patera.
"[A Localization Strategy for Data Assimilation; Application to State Estimation and Parameter Estimation](https://doi.org/10.1137/17M1116830)".
SIAM Journal on Scientific Computing 40(2): B611--B636, 2018.

* C. Himpe, T. Leibner, S. Rave.
"[HAPOD - Fast, Simple and Reliable Distributed POD Computation](https://doi.org/10.11128/arep.55.a55283)".
ARGESIM Report 55 (MATHMOD 2018 Volume): 119--120, 2018.

* B.J. Beach.
"[An Implementation-Based Exploration of HAPOD: Hierarchical Approximate Proper Orthogonal Decomposition](http://hdl.handle.net/10919/81938)".
Master Thesis, Virginia Tech, 2018.

* C. Himpe, T. Leibner, S. Rave, J. Saak.
"[Fast Low-Rank Empirical Cross Gramians](https://doi.org/10.1002/pamm.201710388)".
Proceedings in Applied Mathematics and Mechanics 17: 841--842, 2017.

* C. Himpe, T. Leibner, S. Rave.
"[Comprehensive Memory-Bound Simulations on Single Board Computers](https://doi.org/10.5281/zenodo.814497)".
Extended Abstract, 2nd Conference on Power Aware Computing (PACO), 2017.
