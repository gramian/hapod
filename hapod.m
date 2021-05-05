function [svec,sval,meta] = hapod(data,bound,topo,relax,meta,depth,mysvd)
%%% project: hapod - Hierarchical Approximate POD ( https://git.io/hapod )
%%% version: 3.2 (2021-05-05)
%%% authors: C. Himpe (0000-0003-2194-6754), S. Rave (0000-0003-0439-7212)
%%% license: BSD 2-Clause License (opensource.org/licenses/BSD-2-Clause)
%%% summary: Fast distributed or incremental POD computation.
%
% SYNTAX:
%  [svec,sval,meta] = hapod(data,bound,topo,relax,depth,meta,mysvd)
%
% DESCRIPTION:
%  The hierarchical approximate proper orthogonal decomposition (HAPOD) is a
%  tree-based algorithm to compute low-rank representations of column-wise
%  partitioned data matrices, of which the special cases of incremental HAPOD
%  and distributed HAPOD are implemented. The HAPOD is an error-driven POD
%  method with low communication for distributed memory systems as well as
%  for memory-limited shared memory systems.
%
% ABOUT:
%  Compatible with OCTAVE and MATLAB.
%
% ARGUMENTS:
%   {cell}  data  - snapshot data set, partitioned by column (blocks)
%  {scalar} bound - mean L_2 projection error bound
%  {string} topo  - tree topology, default: 'none'
%                   * 'incr'   - incremental HAPOD (unbalanced binary tree)
%                   * 'incr_1' - incremental HAPOD (child nodes)
%                   * 'incr_r' - incremental HAPOD (root node)
%                   * 'dist'   - distributed HAPOD (star)
%                   * 'dist_1' - distributed HAPOD (child nodes)
%                   * 'dist_r' - distributed HAPOD (root node)
%                   * 'none'   - standard POD
%  {scalar} relax - relaxation parameter in (0,1], default: 0.5
%  {struct} meta  - meta information structure, default: []
%                   * required only for: incr_1, incr_r, dist_r
%  {scalar} depth - total number of levels in tree, default: 2
%                   * required only for: incr_1
%  {handle} mysvd - custom SVD backend, default: 'eco'
%                   * 'eco' - economic (rank revealing) SVD
%                   * 'mos' - method of snapshots
%                   * handle to function with signature: [U,d] = mysvd(data)
%
% RETURNS:
%  {matrix} svec  - POD modes of the data
%  {vector} sval  - singular values to the data
%  {struct} meta  - meta information structure
%                   * nSnapshots - cell array of data column count per node
%                   * nModes - cell array of mode column count per node
%                   * tNode - cell array of timings per node
%
% USAGE:
%  If all data partitions can be passed as the data argument, the topologies:
%  none (standard POD), incr(emental) HAPOD or dist(ributed) HAPOD are
%  applicable. In case only a single partition can be passed, the topologies:
%  incr_1 and dist_1 should be used for the child nodes of the associated HAPOD
%  tree, while the topologies: incr_r and dist_r should be used for the root
%  node. The returned information structure (or a cell-array thereof) can be
%  passed to the parent nodes in the associated HAPOD tree.
%
% CITE AS:
%  C. Himpe, T. Leibner, S. Rave:
%  "Hierarchical Approximate Proper Orthogonal Decomposition".
%  SIAM Journal on Scientific Computing, 40(5): A3267--A3292, 2018.
%
% SEE ALSO:
%  svd, svds, princomp
%
% KEYWORDS:
%  POD, SVD, PCA, KLT, EEF, Model Reduction, Dimension Reduction, Data Science
%
% COPYRIGHT: Christian Himpe, Stephan Rave
%
% Further information: https://git.io/hapod

    if strcmp(data,'version'), svec = 3.2; return; end%if

    % Default arguments
    if nargin<3 || isempty(topo),  topo  = 'none'; end%if
    if nargin<4 || isempty(relax), relax = 0.5; end%if
    if nargin<6 || isempty(depth), depth = 2;  end%if
    if nargin<5, meta  = []; end%if

    % Decode custom SVD
    if nargin<7 || isempty(mysvd) || strcmp(mysvd,'eco')

        mysvd = @eco;
    elseif strcmp(mysvd,'mos')

        mysvd = @mos;
    elseif not(ishandle(mysvd))

        error('hapod: bad mysvd!');
    end%if

    % Argument validation
    assert(isnumeric(bound) && isscalar(bound) && bound>0,'hapod: bad bound!');

    assert(isnumeric(relax) && isscalar(relax) && relax>0 && relax<=1,'hapod: bad relax!');

    assert(not(strcmp(topo,'incr_1')) || isscalar(depth) || depth<2 || (mod(depth,1)==0),'hapod: bad depth!');

    assert(not(strcmp(topo,'incr_r')) || not(strcmp(topo,'dist_r')) || isempty(meta),'hapod: missing meta!');

    % Precompute common quantities
    nodeBound = bound * sqrt(1.0 - relax^2);
    rootBound = bound * relax;

%% Compute HAPOD

    switch lower(topo)

        case 'incr_1' % Incremental HAPOD (Child Nodes)

            [svec,sval,meta] = incr_1(data,nodeBound / sqrt(depth - 1),meta,mysvd);

        case 'incr_r' % Incremental HAPOD (Root Node)

            [svec,sval,meta] = incr_r(data,rootBound,meta,mysvd);

        case 'incr'   % Incremental HAPOD (Complete)

            svec = [];

            for k = 1:size(data,2)-1

                [svec,~,meta] = incr_1({data{k},svec},nodeBound / sqrt(size(data,2) - 1),meta,mysvd);
            end%for

            [svec,sval,meta] = incr_r({data{end},svec},rootBound,meta,mysvd);

        case 'dist_1' % Distributed HAPOD (Child Nodes)

            [svec,sval,meta] = dist_1(data,nodeBound,meta,mysvd);

        case 'dist_r' % Distributed HAPOD (Root Node)

            [svec,sval,meta] = dist_r(data,rootBound,meta,mysvd);

        case 'dist'   % Distributed HAPOD (Complete)

            svec = cell(1,numel(data));
            meta = cell(1,numel(data));

            for k = 1:size(data,2)

                [svec{k},~,meta{k}] = dist_1(data(k),nodeBound,meta,mysvd);
            end%for

            [svec,sval,meta] = dist_r(svec,rootBound,meta,mysvd);

        case 'none'   % Standard POD

            tId = tic();
            meta.nSnapshots = sum(cellfun(@(M) size(M,2),data));
            [svec,sval] = pod(data,sqrt(meta.nSnapshots) * bound,mysvd);
            meta.nModes = size(svec,2);
            meta.tNode = toc(tId);

        otherwise

            error('hapod: unknown topology!');
    end%switch
end

%% LOCAL FUNCTION: incr_1
function [svec,sval,meta] = incr_1(data,incrBound,meta,mysvd)
% summary: Incremental HAPOD: Child nodes

    if isempty(meta)

        meta.nSnapshots{1} = size(data{1},2);
    else

        meta.nSnapshots{end+1} = size(data{1},2) + size(data{2},2) - meta.nModes{end};
    end%if

    tId = tic();
    nodeBound = incrBound * sqrt(sum(cell2mat(meta.nSnapshots)));
    [svec,sval] = pod(data,nodeBound,mysvd);
    svec = bsxfun(@times,svec,sval');

    if not(isfield(meta,'nModes'))

        meta.nModes{1} = size(svec,2);
        meta.tNode{1} = toc(tId);
    else

        meta.nModes{end+1} = size(svec,2);
        meta.tNode{end+1} = toc(tId);
    end%if
end

%% LOCAL FUNCTION: incr_r
function [svec,sval,meta] = incr_r(data,rootBound,meta,mysvd)
% summary: Incremental HAPOD: Root node

    tId = tic();
    meta.nSnapshots{end+1} = size(data{1},2) + size(data{2},2) - meta.nModes{end};
    nodeBound = sqrt(sum(cell2mat(meta.nSnapshots))) * rootBound;
    [svec,sval] = pod(data,nodeBound,mysvd);
    meta.nModes{end+1} = size(svec,2);
    meta.tNode{end+1} = toc(tId);
end

%% LOCAL FUNCTION: dist_1
function [svec,sval,meta] = dist_1(data,distBound,dummy,mysvd)
% summary: Distributed HAPOD: Child nodes

    if iscell(data)

        meta.nSnapshots = size(data{1},2);
    else

        meta.nSnapshots = size(data,2);
    end%if

    tId = tic();
    nodeBound = sqrt(meta.nSnapshots) * distBound;
    [svec,sval] = pod(data,nodeBound,mysvd);
    svec = bsxfun(@times,svec,sval');
    meta.nModes = size(svec,2);
    meta.tNode = toc(tId);
end

%% LOCAL FUNCTION: dist_r
function [svec,sval,meta] = dist_r(data,rootBound,cnfo,mysvd)
% summary: Distributed HAPOD: Root node

    tId = tic();
    meta.nSnapshots = cellfun(@(c) {c.nSnapshots},cnfo);
    meta.nSnapshots{end+1} = sum(cell2mat(meta.nSnapshots));

    nodeBound = sqrt(meta.nSnapshots{end}) * rootBound;
    [svec,sval] = pod(data,nodeBound,mysvd);

    meta.nModes = [cellfun(@(c) {c.nModes},cnfo),size(svec,2)];
    meta.tNode = [cellfun(@(c) {c.tNode},cnfo),toc(tId)];
end

%% LOCAL FUNCTION: pod
function [svec,sval] = pod(data,bound,mysvd)
% summary: Generic proper orthogonal decomposition

    [U,D] = mysvd(cell2mat(data));
    d = flipud(cumsum(flipud(D.^2)));
    K = find(d > bound^2,1,'last');
    sval = D(1:K);
    svec = U(:,1:K);
end

%% LOCAL FUNCTION: eco
function [U,d] = eco(X)
% summary: Standard economic SVD

    [U,D,~] = svd(X,'econ');
    d = diag(D);
end

%% LOCAL FUNCTION: mos
function [U,d] = mos(X)
% summary: Method of snapshots

    [E,V] = eig(X' * X);
    [d,L] = sort(sqrt(abs(diag(V))),'descend');
    U = X * bsxfun(@rdivide,E(:,L),d');
end

