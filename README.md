HAPOD - Hierarchical Approximate Proper Orthogonal Decomposition
================================================================

* HAPOD - Hierarchical Approximate POD ( http://git.io/hapod )
* version: 1.3 ( 2018-02-09 )
* by: Christian Himpe (0000-0003-2194-6754), Stephan Rave (0000-0003-0439-7212)
* under: BSD 2-Clause License (open-source)
* summary: Distributed or incremental POD / SVD computation

# Scope

* (Truncated) Singular Value Decomposition (SVD)
* Proper Orthogonal Decomposition (POD)
* Principal Compoenent Analysis (PCA)
* Empirical Orthogoanl Functions (EOF)
* Empirical Eigenfunctions
* Karhunen Loeve Decomposition
* Model Reduction | Model Order Redction

# Features

* Standard POD
* Incremental HAPOD
* Distributed HAPOD
* Custom SVD backend

# Basic Usage

```
[svec,sval,snfo] = hapod(data,bound,topo,relax,config,mysvd)
```

## Arguments

* `svec` {matrix} POD modes (column vectors)
* `sval` {vector} Singular values
* `snfo` {structure} Information structure

## Return Values

* `data` {cell array} 
* `bound` {scalar} L2 mean projection error bound
* `topo` {string} HAPOD graph topology
    * `'none'` Standard POD
    * `'incr'` Incremental HAPOD 
    * `'incr_0'` Incremental HAPOD (First node only)
    * `'incr_1'` Incremental HAPOD (One intermediary node only)
    * `'incr_r'` Incremental HAPOD (Root node only)
    * `'dist'` Distributed HAPOD
    * `'dist_0'` Distributed HAPOD (First node only)
    * `'dist_1'` Distributed HAPOD (One intermediary node only)
    * `'dist_r'` Distributed HAPOD (Root node only)
* `relax` {scalar} Relaxation parameter in [0,1] by default 0.5.
* `config` {structure} Configuration structure by default empty
* `mysvd` {handle} Function handle to a custom POD method

# Documentation

## Information and Configuration Structure

* `nSnapshots` - Number of data columns passed to this hapod and its children. 
* `nModes` - Number of intermediate modes
* `tNode` - computational time at this hapod's branch

## Custom SVD Backend

Signature, arguments and return values as for Matlab's svd function.

# Cite As

C. Himpe, T. Leibner and S. Rave.
"[Hierarchical Approximate Proper Orthogonal Decomposition](http://arxiv.org/abs/1607.05210)".
Preprint, arXiv math.NA: 1607.05210, 2017.
