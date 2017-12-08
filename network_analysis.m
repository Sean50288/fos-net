function [fosnet] = network_analysis()

% Downloaded from here: https://sites.google.com/site/bctnet/
addpath(genpath('BCT'));

%% Load weighted undirected fos networks
%%% Load in data
alldata = tdfread('data/Fos counts - clean.csv',',');
conditions = {'ALL' 'NJ' 'Nad' 'SJ' 'SA'};

%%% Get data in tidier format
cols = fieldnames(alldata);
fosdata = [];
behavdata = [];
for icol = 1:length(cols)
    if icol>=5 && icol<10
        behavdata = [behavdata alldata.(cols{icol})];
    end
    if icol>=10
        fosdata = [fosdata alldata.(cols{icol})];
    end
end
fosnet.nodelabels = cols(10:end);
fosnet.behavlabels = cols(5:9);

%%% Region info
regions = tdfread('data/region_info.csv',','); % contains labels etc
fosnet.nodegroups = cellstr(regions.AnatGroup); % confirmed that labels in same order
fosnet.labels = cellfun(@horzcat,fosnet.nodegroups,repmat({'-'},29,1),fosnet.nodelabels,'UniformOutput',false);

%% Calculate correlation matrices
%--- Set threshold for creating binary matrices
%--- If empty, matrix will be thresholded by uncorrected p<.05
fosnet.thresh = [];

fprintf('Starting analysis...\n');

%--- Separately for each condition
for con = conditions
    fprintf('%s\n',con{1});
    
    %--- Calculate weighted correlation matrix
    mask = strcmp(alldata.Condition,con);
    rawdata = fosdata(mask,:);
    behdata = behavdata(mask,:);
    
    %--- Options for bootstrapping, if desired
    % num_iter = 1 for no bootstrapping
    %%%% Note that for the paper we do not bootstrap
    num_iter = 1;
    
    for iter=1:num_iter
        
        if mod(iter,50)==0
            fprintf('%s ',num2str(iter));
        end
        
        %--- Randomly sample with replacement
        %%%% Note that for the paper we do not bootstrap
        samp = randsample(size(rawdata,1),size(rawdata,1),true);
        if iter==1 % first iteration - use what's already there
            sampdata = rawdata;
        else
            sampdata = rawdata(samp,:);
        end
        
        %--- Calculate correlation
        [A,p] = corr(sampdata,'type','Kendall','rows','pairwise');
        A(eye(size(A,1))==1) = 0; % set node-to-self relationships to zero for BCT
        
        %--- Calculate weighted binary matrix by thresholding
        %--- If not comparing conditions, can base threshold on significance
        %%%% Note that for the paper, we never use Abin
        Abin = A;
        if ~isempty(fosnet.thresh)
            Abin(A<fosnet.thresh) = 0;
        else
            Abin(p>.05) = 0;
        end
        
        %--- Scale matrices by the mean non-zero weights
        %--- Necessary for comparing across conditions
        %%%% Note that for the paper, we do not compare networks directly
        %%%% across conditions
        Aorig = A;
        A = A./(mean(mean(A(A~=0))));
        Abin = Abin./(mean(mean(Abin(Abin~=0))));
        
        %--- Store all of the variables
        fosnet.(con{1}).condition = con{1};
        %------ original data
        fosnet.(con{1}).rawdata = rawdata;
        fosnet.(con{1}).behavdata = behdata;
        
        if iter==1 % only need to save these out for real matrix
            
            %------ weighted, undirected matrix - all values, before scaling
            fosnet.(con{1}).Aorig = Aorig;
            %------ weighted, undirected matrix - all values
            fosnet.(con{1}).A = A;
            %------ weighted, undirected matrix - thresholded
            fosnet.(con{1}).Abin = Abin;
            %------ uncorrected p-values for each pairwise correlation
            fosnet.(con{1}).pvalues = p;
            %------ number of nodes in network
            fosnet.(con{1}).N = size(A,1);
            
        end
        
        
        %% Calculate graph metrics
        
        %%%%%%%%%%%%%% Perform community detection - modularity functions
        %%% Newman's spectral community detection (essentially deterministic)
        %---- gamma setting: 1 for classic, >1 for smaller mods, <1 for larger mods
        gamma = 1;
        %---- Ci = community structure, net_Q = modularity
        [fosnet.(con{1}).Ci(:,iter),fosnet.(con{1}).net_Q(iter)] = modularity_und(Aorig,gamma);
        
        %%%%%%%%%%%%%% Calculate participation coefficient
        %--- uses community structure Ci from above - use structure from
        %ALL
        [fosnet.(con{1}).node_Ppos(:,iter),fosnet.(con{1}).node_Pneg(:,iter)] = participation_coef_sign(Aorig,fosnet.ALL.Ci(:,iter)); % only works if ALL run first
        
    end
    
end

fname = 'data/all_results_allconds.mat';
fprintf('\nSaving results to %s\n',fname);
save(fname,'fosnet');

%% Write out results for plotting in R

%---- First save out correlation matrix - original (unscaled)
for con=conditions
    origmat = {'seed' 'target' 'connection' 'connection.thresh' 'seed.group' 'target.group' 'seed.module' 'target.module'};
    counter = 1;
    for i=1:size(fosnet.(con{1}).Aorig,1)
        for j=1:size(fosnet.(con{1}).Aorig,2)
            
            counter = counter + 1;
            
            origmat{counter,1} = fosnet.nodelabels{i};
            origmat{counter,2} = fosnet.nodelabels{j};
            origmat{counter,3} = fosnet.(con{1}).Aorig(i,j);
            origmat{counter,4} = fosnet.(con{1}).Abin(i,j);
            origmat{counter,5} = fosnet.nodegroups{i};
            origmat{counter,6} = fosnet.nodegroups{j};
            origmat{counter,7} = fosnet.(con{1}).Ci(i);
            origmat{counter,8} = fosnet.(con{1}).Ci(j);
            
        end
    end
    
    cell2csv(['data/fosnet_matrix_' con{1} '.csv'],origmat);
    
    
    %--- Then save out mean participation coefficient
    
    nodevals = {'metric' 'node' 'nodegroup' 'meanval' 'confmin' 'confmax'};
    counter = 1;
    metric = {'node_Ppos'};
    
    allvals = fosnet.(con{1}).(metric{1});
    meanvals = mean(allvals,2);
    
    for i=1:size(allvals,1)
        tmp = allvals';
        tmp = tmp(:,i);
        tmp = sortrows(tmp);
        if length(tmp)>1 % if multiple iterations
            confmin(i,1) = tmp(round(.025*length(tmp)));
            confmax(i,1) = tmp(round(.975*length(tmp)));
        else % if no iterations
            confmin(i,1) = NaN;
            confmax(i,1) = NaN;
        end
    end
    
    for i=1:size(meanvals,1)
        counter = counter + 1;
        nodevals{counter,1} = metric{1};
        nodevals{counter,2} = fosnet.nodelabels{i};
        nodevals{counter,3} = fosnet.nodegroups{i};
        nodevals{counter,4} = meanvals(i);
        nodevals{counter,5} = confmin(i);
        nodevals{counter,6} = confmax(i);
    end
    
    cell2csv(['data/fosnet_nodevals_' con{1} '.csv'],nodevals);
end

end % main function

