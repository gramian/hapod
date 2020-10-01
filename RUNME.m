function RUNME()
%%% project: hapod - Hierarchical Approximate POD ( https://git.io/hapod )
%%% version: 3.1 (2020-10-01)
%%% authors: C. Himpe (0000-0003-2194-6754), S. Rave (0000-0003-0439-7212)
%%% license: BSD 2-Clause License (opensource.org/licenses/BSD-2-Clause)
%%% summary: Basic tests for incremental HAPOD and distributed HAPOD

%% Generate test data

    randn('seed',1009);			% seed random number generator
    n = 32;					% set number of partitions
    N = n*n;					% set test problem size
    [a,~,c] = svd(randn(N,N));			% compute SVD of random matrix
    b = logspace(0,-16,N)';			% artificial singular values
    s = a*diag(b)*c';	 			% reassign singular values
    S = mat2cell(s,size(s,1),n*ones(n,1));	% split data into partitions
    w = 0.5;					% relaxation parameter omega
    E = sqrt(eps);				% target mean L2 error

    meanl2 = @(U) norm(s-U*(U'*s),'fro') / sqrt(N);	% tested mean L2 error

    disp(' ');
    HAPOD_VERSION = hapod('version')
    disp(' ');
    prescribed_mean_L2 = E
    disp(' ');

    % Compute POD
    disp('REFERENCE POD (ECO):');
    [Uref,Dref,~] = hapod(S,E,'none');
    mean_L2 = meanl2(Uref)
    error_bound_OK = mean_L2 <= E
    num_global_modes = size(Uref,2)
    disp(' ');

    % Compute POD
    disp('REFERENCE POD (MOS):');
    [Umos,Dmos,~] = hapod(S,E,'none',[],[],[],'mos');
    mean_L2 = meanl2(Umos)
    error_bound_OK = mean_L2 <= E
    num_global_modes = size(Umos,2)
    disp(' ');

    % Determine mode bound
    [Uloc,~,~] = hapod(S,E*sqrt(1-w^2)/sqrt(n-1),'none');
    local_mode_bound = size(Uloc,2)
    disp(' ');

%% Incremental HAPOD

    % Test bulk incremental HAPOD
    disp('INCREMENTAL HAPOD (BULK):');
    [Uincr,Dincr,Cincr] = hapod(S,E,'incr',w);
    mean_L2 = meanl2(Uincr)
    error_bound_OK = mean_L2 <= E
    num_global_modes = size(Uincr,2)
    max_local_modes = max(cell2mat(Cincr.nModes(1:end-1)))
    mode_bound_OK = max_local_modes <= local_mode_bound
    disp(' ');

    % Test chunk incremental HAPOD
    disp('INCREMENTAL HAPOD (CHUNK):');
    u1 = [];
    c1 = [];

    for k = 1:n-1

        [u1,~,c1] = hapod({S{k},u1},E,'incr_1',w,c1,n);
    end%for

    [Uincr,Dincr,Cincr] = hapod({S{n},u1},E,'incr_r',w,c1,n);

    mean_L2 = meanl2(Uincr)
    error_bound_OK = mean_L2 <= E
    num_global_modes = size(Uincr,2)
    max_local_modes = max(cell2mat(Cincr.nModes(1:end-1)))
    mode_bound_OK = max_local_modes <= local_mode_bound
    disp(' ');

%% Distributed HAPOD

    % Test bulk distributed HAPOD
    disp('DISTRIBUTED HAPOD (BULK):');
    [Udist,Ddist,Cdist] = hapod(S,E,'dist',w);
    mean_L2 = meanl2(Udist)
    error_bound_OK = mean_L2 <= E
    num_global_modes = size(Udist,2)
    max_local_modes = max(cell2mat(Cdist.nModes(1:end-1)))
    mode_bound_OK = max_local_modes <= local_mode_bound
    disp(' ');

    % Test chunk distributed HAPOD
    disp('DISTRIBUTED HAPOD (CHUNK):');
    u2 = cell(1,n);
    c2 = cell(1,n);

    for k = 1:n

        [u2{k},~,c2{k}] = hapod(S(k),E,'dist_1',w);
    end%for

    [Udist,Ddist,Cdist] = hapod(u2,E,'dist_r',w,c2);

    mean_L2 = meanl2(Udist)
    error_bound_OK = mean_L2 <= E
    num_global_modes = size(Udist,2)
    max_local_modes = max(cell2mat(Cdist.nModes(1:end-1)))
    mode_bound_OK = max_local_modes <= local_mode_bound
    disp(' ');

%% Distributed-of-Incremental HAPOD

    % Test dist-of-incr HAPOD
    disp('DIST-OF-INC HAPOD (CHUNK):');
    strand = {1:ceil(n/2),ceil(n/2)+1:n};
    M = numel(strand);
    u3 = repmat({[]},1,M);
    c3 = repmat({[]},1,M);

    for m = 1:M

        for k = strand{m}

            [u3{m},~,c3{m}] = hapod({S{k},u3{m}},E,'incr_1',w,c3{m},n);
        end%for

        c3{m}.nSnapshots =  sum([c3{m}.nSnapshots{:}]);
        c3{m}.nModes = max([c3{m}.nModes{:}]);
        c3{m}.tNode = sum([c3{m}.tNode{:}]);
    end%for

    [Udofi,Ddofi,Cdofi] = hapod(u3,E,'dist_r',w,c3);

    mean_L2 = meanl2(Udofi)
    error_bound_OK = mean_L2 <= E
    num_global_modes = size(Udofi,2)
    max_local_modes = max(cell2mat(Cdofi.nModes(1:end-1)))
    mode_bound_OK = max_local_modes <= local_mode_bound
    disp(' ');

%% Plot Results

    frobsv = @(D) sqrt(flipud(cumsum(flipud(b(1:numel(D)).^2))));

    % Plot Singular Values
    figure('Name','HAPOD Test')
    subplot(1,2,1);
    semilogy(1:numel(Dref),b(1:numel(Dref)),'LineWidth',2);
    hold on;
    semilogy(1:numel(Dref),Dref,'--','LineWidth',2);
    semilogy(1:numel(Dmos),Dmos,'--','LineWidth',2);
    semilogy(1:numel(Dincr),Dincr,'--','LineWidth',2);
    semilogy(1:numel(Ddist),Ddist,'--','LineWidth',2);
    semilogy(1:numel(Ddofi),Ddofi,'--','LineWidth',2);
    hold off;
    xlim([1,max([numel(Dref),numel(Dmos),numel(Dincr),numel(Ddist),numel(Ddofi)])]);
    ylabel('Singular Values');
    legend('Exact','Reference POD','Method of Snapshots','Incremental HAPOD','Distributed HAPOD','Distributed-of-Incremental HAPOD','Location','SouthOutside');

    % Plot Frobenius Norm
    subplot(1,2,2);
    semilogy(1:numel(Dref),frobsv(b(1:numel(Dref))),'LineWidth',2);
    hold on;
    semilogy(1:numel(Dref),frobsv(Dref),'--','LineWidth',2);
    semilogy(1:numel(Dmos),frobsv(Dmos),'--','LineWidth',2);
    semilogy(1:numel(Dincr),frobsv(Dincr),'--','LineWidth',2);
    semilogy(1:numel(Ddist),frobsv(Ddist),'--','LineWidth',2);
    semilogy(1:numel(Ddofi),frobsv(Ddofi),'--','LineWidth',2);
    hold off;
    xlim([1,max([numel(Dref),numel(Dmos),numel(Dincr),numel(Ddist),numel(Ddofi)])]);
    ylabel('Frobenius Norm');
    legend('Exact','Reference POD','Method of Snapshots','Incremental HAPOD','Distributed HAPOD','Distributed-of-Incremental HAPOD','Location','SouthOutside');
end
