function MAPRED()
%%% project: hapod - Hierarchical Approximate POD ( https://git.io/hapod )
%%% version: 3.2 (2021-05-05)
%%% authors: C. Himpe (0000-0003-2194-6754), S. Rave (0000-0003-0439-7212)
%%% license: BSD 2-Clause License (opensource.org/licenses/BSD-2-Clause)
%%% summary: Basic test of (parallel) distributed HAPOD via MapReduce

%% Generate Test Data

    randn('seed',1009);
    n = 32;
    N = n*n;
    [a,~,c] = svd(randn(N,N));
    b = logspace(0,-16,N)';
    S = a*diag(b)*c';

    E = sqrt(eps);
    w = 0.5;

%% MapReduce Setup

    % Define datastore
    ds = arrayDatastore(S,'IterationDimension',2,'ReadSize',n);

    % Define mapper
    function hapod_mapper(data,info,intermKVStore)

        [u,~,c] = hapod(data',E,'dist_1',w);
        add(intermKVStore,'leaf',{u,c});
    end

    % Define reducer
    function hapod_reducer(intermKey,intermValIter,outKVStore)

        u = {};
        c = {};

        while hasnext(intermValIter)

            value = getnext(intermValIter);
            u{end+1} = value{1};
            c{end+1} = value{2};
        end%while

        [U,D,C] = hapod(u,E,'dist_r',w,c);
        addmulti(outKVStore,{'root_singvec','root_singval','root_info'},{U,D,C});
    end

    % Apply map and reduce
    mr = mapreduce(ds,@hapod_mapper,@hapod_reducer).readall();

    singvals = mr.Value{find(strcmp(mr.Key,'root_singval'))};

%% Plot Results

    figure;
    semilogy(1:numel(singvals),b(1:numel(singvals)),'LineWidth',2);
    hold on;
    semilogy(1:numel(singvals),singvals,'LineWidth',2,'LineStyle','--');
    hold off;
    xlim([1,numel(singvals)]);
    ylabel('Singular Values');
    legend('Exact','Distributed HAPOD','Location','SouthOutside');
end

