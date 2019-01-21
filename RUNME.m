function RUNME()
%%% project: hapod - Hierarchical Approximate POD ( https://git.io/hapod )
%%% version: 2.0 ( 2019-01-21 )
%%% authors: C. Himpe ( 0000-0003-2194-6754 ), S. Rave ( 0000-0003-0439-7212 )
%%% license: BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%%% summary: basic test for incremental HAPOD and distributed HAPOD

%% Generate test data

    randn('seed',1009);			   % seed random number generator
    n = 16;				   % set number of partitions
    N = n*n;				   % set test problem size
    [a,b,c] = svd(randn(N,N));		   % compute SVD of normal random matrix
    s = a*diag(logspace(0,-16,N))*b';	   % reassign singular values
    S = mat2cell(s,size(s,1),n*ones(n,1)); % split data matrix into partitions
    w = 0.1;				   % relaxation parameter omega
    E = 1e-8;				   % target mean L2 projection error
    disp(' ');

    % Compute default POD
    disp('REFERENCE POD:');
    [U,D,C] = hapod(S,E,'none');
    [ug,dg,cg] = hapod(S,E*w,'none');
    prescribed_mean_L2 = E
    mean_L2 = norm(s-U*U'*s,'fro')/sqrt(N)
    global_mode_bound = size(ug,2)
    disp(' ');

%% Incremental HAPOD

    % Test bulk incremental HAPOD
    disp('INCREMENTAL HAPOD (BULK):');
    [U,D,C] = hapod(S,E,'incr',w);
    mean_L2 = norm(s-U*U'*s,'fro')/sqrt(N)
    num_global_modes = size(U,2)
    max_local_modes = max(cell2mat(C.nModes(1:end-1)))
    disp(' ');

    % Test chunk incremental HAPOD
    disp('INCREMENTAL HAPOD (CHUNK):');
    u1 = [];
    c1 = struct('nLevels',n);

    for k=1:n-1

        [u1,d1,c1] = hapod({S{k},u1},E,'incr_1',w,c1);
    end

    [U,D,C] = hapod({S{n},u1},E,'incr_r',w,c1);

    mean_L2 = norm(s-U*U'*s,'fro')/sqrt(N)
    num_global_modes = size(U,2)
    max_local_modes = max(cell2mat(C.nModes(1:end-1)))
    disp(' ');

%% Distributed HAPOD

    % Test bulk distributed HAPOD
    disp('DISTRIBUTED HAPOD (BULK):');
    [U,D,C] = hapod(S,E,'dist',w);
    mean_L2 = norm(s-U*U'*s,'fro')/sqrt(N)
    num_global_modes = size(U,2)
    max_local_modes = max(cell2mat(C.nModes(1:end-1)))
    disp(' ');

    % Test chunk distributed HAPOD
    disp('DISTRIBUTED HAPOD (CHUNK):');
    u2 = cell(1,n);

    for k = 1:n

        [u2{k},d2,c2{k}] = hapod(S(k),E,'dist_1',w);
    end

    [U,D,C] = hapod(u2,E,'dist_r',w,c2);

    mean_L2 = norm(s-U*U'*s,'fro')./sqrt(N)
    num_global_modes = size(U,2)
    max_local_modes = max(cell2mat(C.nModes(1:end-1)))
    disp(' ');

%% Distributed-of-Incremental HAPOD

    % Test dist-of-inc HAPOD
    disp('DIST-OF-INC HAPOD (CHUNK):');
    strand = {[1:ceil(n/2)],ceil(n/2)+1:n};
    M = numel(strand);
    u3 = repmat({[]},1,M);
    c3 = repmat({struct('nLevels',n+1)},1,M);

    for m = 1:M

        for k = strand{m}

            [u3{m},d3,c3{m}] = hapod({S{k},u3{m}},E,'incr_1',w,c3{m});
        end

        c3{m}.nSnapshots =  sum([c3{m}.nSnapshots{:}]);
        c3{m}.nModes = max([c3{m}.nModes{:}]);
        c3{m}.tNode = sum([c3{m}.tNode{:}]);
    end

    [U,D,C] = hapod(u3,E,'dist_r',w,c3);

    mean_L2 = norm(s-U*U'*s,'fro')./sqrt(N)
    num_global_modes = size(U,2)
    max_local_modes = max(cell2mat(C.nModes(1:end-1)))
    disp(' ');
end
