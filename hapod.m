function [svec,sval,snfo] = hapod(data,bound,type,relax,config,mysvd)
%%% project: hapod - Hierarchical Approximate POD ( https://git.io/hapod )
%%% version: 2.0 ( 2019-01-21 )
%%% authors: C. Himpe ( 0000-0003-2194-6754 ), S. Rave ( 0000-0003-0439-7212 )
%%% license: BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%%% summary: Distributed or incremental POD / SVD computation
%
% SYNTAX:
%  [svec,sval,snfo] = hapod(data,bound,type,relax,config,mysvd)
%
% DESCRIPTION:
%  The hierarhical approximate proper orthogonal decomposition is a tree-based
%  algorithm to compute low-rank representations of column partitioned data
%  matrices, of which the special cases of Incremental HAPOD and Distributed
%  HAPOD are implemented.
%
% ABOUT:
%  Compatible with OCTAVE and MATLAB.
%
% ARGUMENTS:
%    (cell) data   - snapshot data set, partitioned by column (blocks)
%  (scalar) bound  - mean L_2 projection error bound
%  (string) type   - tree type, default: 'none'
%                  * 'incr' - incremental HAPOD
%                  * 'incr_1' - incremental HAPOD (non-root node)
%                  * 'incr_r' - incremental HAPOD (root node)
%                  * 'dist' - distributed HAPOD
%                  * 'dist_1' - distributed HAPOD (non-root node)
%                  * 'dist_r' - distributed HAPOD (root node)
%                  * 'none' - standard POD
%  (scalar) relax  - relaxation parameter in (0,1], default: 0.5
%  (struct) config - partial HAPOD configuration, default: {}
%                  * nLevels - total number of levels in tree
%                  * nSnapshots - cell array of data column count per node
%                  * nModes - cell array of mode column count per node
%                  * tNode - cell array of per node timings
%  (handle) mysvd  - custom SVD backend, default: []
%                  * signature: [leftvec,singval] = mysvd(data)
%
% RETURNS:
%  (matrix) svec - POD modes of the data
%  (vector) sval - singular values to the data
%  (struct) snfo - information structure
%
% USAGE:
%  If all data partitions can be passed as the data argument, the types: none 
%  (standard POD), incr(emental) HAPOD or dist(ributed) HAPOD are applicable.
%  In case only a single partition can be passed, the types: incr_1 and dist_1
%  should be used for the non-root nodes of the associated HAPOD tree, while 
%  the types: incr_r and dist_r should be used for the root node. The returned
%  information structure (or a cell-array thereof) can be passed to the parent
%  nodes in the associated HAPOD tree. 
%
% CITATION:
%  C. Himpe, T. Leibner and S. Rave.
%  "Hierarchical Approximate Proper Orthogonal Decomposition".
%  SIAM Journal on Scientific Computing, 40(5): A3267--A3292, 2018.
%
% SEE ALSO:
%  svd, svds, princomp
%
% KEYWORDS:
%  POD, BPOD, SVD, PCA, Model Reduction, Dimension Reduction
%
% Further information: https://git.io/hapod

    if(strcmp(data,'version')), svec = 2.0; return; end%if

    if( (nargin<3) || isempty(type)  ), type  = 'none'; end%if
    if( (nargin<4) || isempty(relax) ), relax = 0.5; end%if
    if( nargin<5 ), config = []; end%if
    if( nargin<6 ), mysvd = []; end%if

    % Precompute common quantities
    scaledBound = bound * sqrt(1.0 - relax * relax);
    relaxBound = bound * relax;

%% Compute HAPOD

    switch(lower(type))

        case 'incr_1' % Incremental HAPOD (Non-Root Node)

            [svec,sval,snfo] = incr_1(data,scaledBound,config,mysvd);

        case 'incr_r' % Incremental HAPOD (Root Node Only)

            [svec,sval,snfo] = incr_r(data,relaxBound,config,mysvd);

        case 'incr'   % Incremental HAPOD (Complete)

            svec = [];
            snfo.nLevels = size(data,2);

            for k = 1:size(data,2)-1

                [svec,sval,snfo] = incr_1({data{k},svec},scaledBound,snfo,mysvd);
            end%for

            [svec,sval,snfo] = incr_r({data{end},svec},relaxBound,snfo,mysvd);

        case 'dist_1' % Distributed HAPOD (Non-Root Node)

            [svec,sval,snfo] = dist_1(data,scaledBound,config,mysvd);

        case 'dist_r' % Distributed HAPOD (Root Node Only)

            [svec,sval,snfo] = dist_r(data,relaxBound,config,mysvd);

        case 'dist'   % Distributed HAPOD (Complete)

            for k = 1:size(data,2)

                [svec{k},sval,snfo{k}] = dist_1(data(k),scaledBound,config,mysvd);
            end%for

            [svec,sval,snfo] = dist_r(svec,relaxBound,snfo,mysvd);

        otherwise % Standard POD

            tId = tic();
            snfo.nSnapshots = sum(cellfun(@(M) size(M,2),data));
            snfo.nLevels = 1;
            normBound = sqrt(snfo.nSnapshots) * bound;
            [svec,sval] = pod(data,normBound,mysvd);
            snfo.nModes = size(svec,2);
            snfo.tNode = toc(tId);
    end%switch
end

% LOCAL FUNCTION: incr_1
function [svec,sval,config] = incr_1(data,scaledBound,config,mysvd)
% summary: Incremental HAPOD: Non-root node only

    if(not(isfield(config,'nSnapshots')))
        config.nSnapshots{1} = size(data{1},2);
    else
        config.nSnapshots{end+1} = size(data{1},2) + size(data{2},2) - config.nModes{end};
    end%if

    tId = tic();
    nodeBound = scaledBound * sqrt(sum(cell2mat(config.nSnapshots)) / (config.nLevels - 1));
    [svec,sval] = pod(data,nodeBound,mysvd);
    svec = bsxfun(@times,svec,sval);

    if(not(isfield(config,'nModes')))
        config.nModes{1} = size(svec,2);
        config.tNode{1} = toc(tId);
    else
        config.nModes{end+1} = size(svec,2);
        config.tNode{end+1} = toc(tId);
    end%if
end

% LOCAL FUNCTION: incr_r
function [svec,sval,config] = incr_r(data,relaxBound,config,mysvd)
% summary: Incremental HAPOD: Root node only

    tId = tic();
    config.nSnapshots{end+1} = size(data{1},2) + size(data{2},2) - config.nModes{end};
    rootBound = sqrt(sum(cell2mat(config.nSnapshots))) * relaxBound;
    [svec,sval] = pod(data,rootBound,mysvd);
    config.nModes{end+1} = size(svec,2);
    config.tNode{end+1} = toc(tId);
end

% LOCAL FUNCTION: dist_1
function [svec,sval,snfo] = dist_1(data,scaledBound,config,mysvd)
% summary: Distributed HAPOD: Non-root node only

    if(iscell(data))
        snfo.nSnapshots = size(data{1},2);
    else
        snfo.nSnapshots = size(data,2);
    end%if

    tId = tic();
    nodeBound = sqrt(snfo.nSnapshots) * scaledBound;
    [svec,sval] = pod(data,nodeBound,mysvd);
    svec = bsxfun(@times,svec,sval);
    snfo.nLevels = 2;
    snfo.nModes = size(svec,2);
    snfo.tNode = toc(tId);
end

% LOCAL FUNCTION: dist_r
function [svec,sval,snfo] = dist_r(data,relaxBound,config,mysvd)
% summary: Distributed HAPOD: Root node only

    tId = tic();
    snfo.nSnapshots = cellfun(@(c) {c.nSnapshots},config);
    snfo.nSnapshots{end+1} = sum(cell2mat(snfo.nSnapshots));
    rootBound = sqrt(snfo.nSnapshots{end}) * relaxBound;
    [svec,sval] = pod(data,rootBound,mysvd);
    snfo.nLevels = 2;
    snfo.nModes = [cellfun(@(c) {c.nModes},config),size(svec,2)];
    snfo.tNode = [cellfun(@(c) {c.tNode},config),toc(tId)];
end

% LOCAL FUNCTION: pod
function [svec,sval] = pod(data,bound,mysvd)
% summary: Generic proper orthogonal decomposition

    X = cell2mat(data);

    if(isempty(mysvd))          % Plain economic (rank-revealing) SVD
        [U,D,V] = svd(X,'econ');
        D = diag(D);
    elseif(strcmp(mysvd,'mos')) % Method of snapshots
        [E,V] = eig(X'*X,'vector');
        [D,I] = sort(sqrt(abs(V)),'descend');
        U = X * bsxfun(@rdivide,E(:,I),D');
    else                        % Custom SVD
        [U,D,V] = mysvd(X);
        if(size(D,2) > 1), D = diag(D); end%if
    end%if

    d = flipud(cumsum(flipud(D.^2)));

    K = find(d <= (bound*bound),1);
    if(isempty(K)), K = size(X,2) + 1; end%if
    sval = D(1:K-1)';
    svec = U(:,1:K-1);
end
