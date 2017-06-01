function [svec,sval,snfo] = hapod(data,bound,topo,relax,config,mysvd)
%%% project: hapod - Hierarchical Approximate POD ( http://git.io/hapod )
%%% version: 1.1 ( 2017-06-01 )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 ),
%%%          Stephan Rave ( 0000-0003-0439-7212 )
%%% license: BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%
%% SYNTAX:
%   [svec,sval,snfo] = hapod(data,bound,topo,relax,config,mysvd)
%
%% SUMMARY:
%   Compatible with OCTAVE and MATLAB.
%
%% ARGUMENTS:
%     (cell) data   - snapshot data set, partitioned by column (blocks)
%   (scalar) bound  - mean L_2 projection error bound
%   (string) topo   - tree topology, default: 'none'
%                   * 'incr' - incremental HAPOD
%                   * 'incr_0' - incremental HAPOD
%                   * 'incr_1' - incremental HAPOD
%                   * 'incr_r' - incremental HAPOD
%                   * 'dist' - distributed HAPOD
%                   * 'dist_1' - distributed HAPOD
%                   * 'dist_r' - distributed HAPOD
%                   * 'none' - standard POD
%   (scalar) relax  - relaxation parameter in (0,1], default: 0.5
%   (struct) config - partial HAPOD configuration, default: []
%                   * nSets - number of snapshot sets
%                   * nSnapshots - vector of data column count per node
%                   * nModes - vector of mode column count per node
%                   * nodeIndex - current nodes global index
%                   * tNode - vector of per node timings
%   (handle) mysvd  - custom SVD backend, default: []
%                   * signature: [leftvec,singval] = mysvd(data)
%
%% RETURNS:
%   (matrix) svec - POD modes of the data
%   (vector) sval - singular values to the data
%   (struct) snfo - information structure (nSets, nSnapshots, nModes)
%
%% CITATION:
%   C. Himpe, T. Leibner and S. Rave.
%   "Hierarchical Approximate Proper Orthogonal Decomposition".
%   Preprint, arXiv math.NA: 1607.05210, 2016.
%
%% SEE ALSO:
%   svd, svds, princomp
%
%% KEYWORDS:
%   POD, BPOD, SVD, tSVD, PCA, MOR, Model Reduction
%
% Further information: <http://git.io/hapod>
%*
    if(strcmp(data,'version')), svec = 1.1; return; end;

    if( (nargin<3) || isempty(topo)  ), topo  = 'none'; end;
    if( (nargin<4) || isempty(relax) ), relax = 0.5; end;
    if( nargin<6 ), mysvd = []; end;

    nSets  = size(data,2);
    nSnapshots  = zeros(nSets,1);
    nModes = zeros(nSets,1);
    tNode = zeros(nSets,1);
    scaledBound = bound * sqrt(1.0 - relax^2);

    switch(lower(topo))

        case 'incr' % Incremental HAPOD

            leafBase = [];

            for nodeIndex = 1:nSets-1

                tId = tic();
                nSnapshots(nodeIndex) = size(data{nodeIndex},2);
                leafBound = scaledBound * sqrt(sum(nSnapshots) / ((nSets - 1) * nSnapshots(nodeIndex)));
                [leafBase,leafValues] = pod({leafBase,data{nodeIndex}},leafBound,mysvd);
                leafBase = leafBase.*leafValues;
                nModes(nodeIndex) = size(leafBase,2);
                tNode(nodeIndex) = toc(tId);
            end

            tId = tic();
            nSnapshots(end) = size(data{end},2);
            rootBound = bound * relax * sqrt(sum(nSnapshots)/nSnapshots(nodeIndex));
            [svec,sval] = pod({leafBase,data{end}},rootBound,mysvd);
            tNode(nSets) = toc(tId);
            snfo = struct('nSets',nSets,'nSnapshots',nSnapshots,'nModes',nModes,'tNode',tNode);

        case 'incr_0' % Incremental HAPOD (First Node Only)

            config.nSnapshots = zeros(1,config.nSets);
            config.nSnapshots(1) = size(data{1},2);
            config.nModes = zeros(1,config.nSets);
            config.tNode = zeros(1,config.nSets);

            tId = tic();
            leafBound = scaledBound * sqrt(1.0 / (config.nSets - 1));
            [svec,sval] = pod(data,leafBound,mysvd);
            svec = svec.*sval;
            config.nModes(1) = size(svec,2);
            config.tNode(1) = toc(tId);
            snfo = config;

        case 'incr_1' % Incremental HAPOD (One Node Only)

            tId = tic();
            config.nSnapshots(config.nodeIndex) = size(data{1},2) + size(data{2},2) - config.nModes(config.nodeIndex-1);
            leafBound = scaledBound * sqrt(sum(config.nSnapshots) / ((config.nSets - 1) * config.nSnapshots(config.nodeIndex)));
            [svec,sval] = pod(data,leafBound,mysvd);
            svec = svec.*sval;
            config.nModes(config.nodeIndex) = size(svec,2);
            config.tNode(config.nodeIndex) = toc(tId);
            snfo = config;

        case 'incr_r' % Incremental HAPOD (Root Node Only)

            tId = tic();
            config.nSnapshots(config.nodeIndex) = size(data{1},2) + size(data{2},2) - config.nModes(config.nodeIndex-1);
            rootBound = bound * relax * sqrt(sum(config.nSnapshots)/config.nSnapshots(config.nodeIndex));
            [svec,sval] = pod(data,rootBound,mysvd);
            config.nModes(config.nodeIndex) = size(svec,2);
            config.tNode(config.nodeIndex) = toc(tId);
            snfo = config;

        case 'dist' % Distributed HAPOD

            leafBases = cell(1,nSets);

            for nodeIndex = 1:nSets

                tId = tic();
                nSnapshots(nodeIndex) = size(data{nodeIndex},2);
                [leafBases{nodeIndex},leafValues] = pod(data(nodeIndex),scaledBound,mysvd);
                leafBases{nodeIndex} = leafBases{nodeIndex}.*leafValues;
                nModes(nodeIndex) = size(leafBases{nodeIndex},2);
                tNode(nodeIndex) = toc(tId);
            end

            [svec,sval] = pod(leafBases,bound * relax,mysvd);
            snfo = struct('nSets',nSets,'nSnapshots',nSnapshots,'nModes',nModes,'tNode',tNode);

        case 'dist_1' % Distributed HAPOD (One Node Only)

            config.nSnapshots = size(data{1},2);

            tId = tic();
            [svec,sval] = pod(data,scaledBound,mysvd);
            svec = svec.*sval;
            config.nModes = size(svec,2);
            config.tNode = toc(tId);
            snfo = config;

        case 'dist_r' % Distributed HAPOD (Root Node Only)

            tId = tic();
            [svec,sval] = pod(data,bound * relax,mysvd);
            snfo = struct('nSnapshots',[config.nSnapshots],'nModes',[config.nModes],'tNode',[config.tNode]);

        case 'none' % Standard POD

            tId = tic();
            [svec,sval] = pod(data,bound,mysvd);
            tNode = toc(tId);
            snfo = [];

        otherwise

            error('hapod: unknown HAPOD topology!');
    end
end

function [svec,sval] = pod(data,bound,mysvd)

    X = cell2mat(data);

    if(isempty(mysvd))
        [U,D,V] = svd(X,'econ');
    else
        [U,D,V] = mysvd(X);
    end

    D = diag(D);
    d = flipud(cumsum(flipud(D.^2)));

    K = find(d <= size(X,2)*bound^2,1);
    if(isempty(K)), K = size(X,2); end;
    sval = D(1:K)';
    svec = U(:,1:K);
end
