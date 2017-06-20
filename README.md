HAPOD - Hierarchical Approximate Proper Orthogonal Decomposition
================================================================

* HAPOD - Hierarchical Approximate Proper Orthogonal Decomposition ( http://git.io/hapod )
* version: 1.2 (2017-20-06)
* by: Christian Himpe (0000-0003-2194-6754), Stephan Rave (0000-0003-0439-7212)
* under: BSD 2-Clause License (open-source)
* summary: 

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
* Custom POD backend

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
* `bound` {scalar} L2 projection error bound
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

## Custom POD Backend

## Distributed Memory Parallelization

# Cite As

C. Himpe, T. Leibner and S. Rave.
"Hierarchical Approximate Proper Orthogonal Decomposition".
Preprint, arXiv math.NA: 1607.05210, 2016.
