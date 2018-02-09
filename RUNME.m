function RUNME()
%%% project: hapod - Hierarchical Approximate POD ( http://git.io/hapod )
%%% version: 1.3 ( 2018-02-09 )
%%% authors: C. Himpe ( 0000-0003-2194-6754 ), S. Rave ( 0000-0003-0439-7212 )
%%% license: BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%%% summary: basic test for incremental HAPOD and distributed HAPOD
%

%% Generate test data

    randn('seed',1009);
    n = 16;
    [a,b,c] = svd(randn(n*n,n*n));
    s = a*diag(logspace(0,-16,n*n))*b';
    S = mat2cell(s,size(s,1),n*ones(n,1));

    w = 0.1;
    E = 1e-8;

    % Compute default POD
    [U,D,C] = hapod(S,E,'none');
    PRESCRIBED_MEAN_L2 = E
    POD_MEAN_L2 = norm(s-U*U'*s,'fro')/sqrt(n*n)
    disp('');

%% Incremental HAPOD

    % Compute global mode bound POD for incremental HAPOD
    [ug,dg,cg] = hapod(S,E*w,'none');
    iHAPOD_GLOBAL_MODE_BOUND = size(ug,2)

    % Compute local mode bound POD for incremental HAPOD
    [ul,dl,cl] = hapod(S,E*sqrt(1.0-w^2)/sqrt(n*n-1),'none');
    iHAPOD_LOCAL_MODE_BOUND = size(ul,2)
    disp('');

    % Test full incremental HAPOD
    [U,D,C] = hapod(S,E,'incr',w);
    iHAPOD_MEAN_L2 = norm(s-U*U'*s,'fro')/sqrt(n*n)
    iHAPOD_MODES = size(U,2)
    iHAPOD_MAX_MODES = max(C.nModes)
    disp('');

    % Test partial incremental HAPOD
    c1.nSets = n;
    [u1,d1,c1] = hapod(S(1),E,'incr_0',w,c1);

    for k=2:n-1

        c1.nodeIndex = k;
        [u1,d1,c1] = hapod({S{1},u1},E,'incr_1',w,c1);
    end

    c1.nodeIndex = n;
    [u1,d1,c1] = hapod({S{1},u1},E,'incr_r',w,c1);

    iHAPOD_MEAN_L2 = norm(s-U*U'*s,'fro')/sqrt(n*n)
    iHAPOD_MODES = size(U,2)
    iHAPOD_MAX_MODES = max(C.nModes)
    disp('');

%% Distributed HAPOD

    % Compute global mode bound POD for incremental HAPOD
    [ug,dg,cg] = hapod(S,E*w,'none');
    dHAPOD_GLOBAL_MODE_BOUND = size(ug,2)

    % Compute local mode bound POD for incremental HAPOD
    [ul,dl,cl] = hapod(S,E*sqrt(1.0-w^2),'none');
    dHAPOD_LOCAL_MODE_BOUND = size(ul,2)
    disp('');

    % Test full distributed HAPOD
    [U,D,C] = hapod(S,E,'dist',w);
    dHAPOD_MEAN_L2 = norm(s-U*U'*s,'fro')/sqrt(n*n)
    dHAPOD_MODES = size(U,2)
    dHAPOD_MAX_MODES = max(C.nModes)
    disp('');

    % Test partial distributed HAPOD 
    u2t = cell(1,n);
    d2t = cell(1,n);

    nSnapshots = 0;
    nModes = 0;
    tNode = 0;

    for k=1:n
 
        [u2t{k},d2t{k},c2t] = hapod(S(k),E,'dist_1',w);
        nSnapshots = nSnapshots + c2t.nSnapshots;
        nModes = max(nModes,c2t.nModes);
        tNode = tNode + c2t.tNode;
    end

    [u2,d2,c2] = hapod(u2t,E,'dist_r',w,struct('nSnapshots',nSnapshots,'nModes',nModes,'tNode',tNode));

    dHAPOD_MEAN_L2 = norm(s-u2*u2'*s,'fro')./sqrt(n*n)
    dHAPOD_MODES = size(u2,2)
    dHAPOD_MAX_MODES = max(c2.nModes)
end
