function [NND, MM_Eta, MM_T, MM_R, nComp_best] = model_nnd(vCat,NNDpar,Model,OptArgs)
%
%   Perfomrs NND analysis of the catalog of events
%   vCat - the catalogue of events
%   Model - the structure that defines the seismicity model
%   OprArgs.FitModel - which mixture models to fit: 'GMM', 'Weilbul' 
%   OptArgs.FMD - FMD to use: 'GR', 'Exp', 'TapTapPareto', 'TapParetoPareto'
%   OptArgs.ThresholdType - 'intersect' the thresholds are computed as intersection of GMM compenets;
%                           'minsaddle' the thresholds are computed as the local minima for each saddle if present
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 25 May 2022
%   ...
%   version 2.4.0, 27 October 2024
%
    arguments
        vCat double
        NNDpar struct
        Model struct
        OptArgs.FMD char = 'GR'                  % frequency-magnitude distribution to rescale \eta, T, and R: 'GR', 
        OptArgs.NNDmethod char = 'OriginalNND'   % which NND declustering method to use: 'OriginalNND', 'ZBZ2020'
        OptArgs.Distance char = 'Radial'         % distance between events: 'Radial', 'Euclidean', 'Euclidean3D'   
        OptArgs.TimeUnit char = 'day'
        OptArgs.Display char = 'on'              % 'on'/'off' to plot the distributions or not 
    end
    
    vCat_NND = vCat;
    if strcmp(OptArgs.TimeUnit,'year')
        vCat_NND(:,1) = vCat(:,1)/365.0; % time used is in years
    elseif strcmp(OptArgs.TimeUnit,'day')
    end

    % perform the NND analysis
    tic
    NND  = nnd(vCat_NND,NNDpar,'FMD',OptArgs.FMD,'Distance',OptArgs.Distance);
    toc

    vEta = NND.vEta;
    vT   = NND.vT;
    vR   = NND.vR;
    %DiG = nnd_digraph(NND,vCat);
    %plot_nnd_clusters(DiG,vCat,Model);
    %save_nnd_values(NND,'nnd_values_without_threshold.dat');

    % limit the ranges for fitting mixture models only
    indx = (vEta >= NNDpar.EtaLim(1)) & (vEta <= NNDpar.EtaLim(2)); % limit \eta values between bounds
    indx = indx & (vT >= NNDpar.TLim(1)) & (vT <= NNDpar.TLim(2)); % limit T values between bounds
    vEta = vEta(indx);
    vT   = vT(indx);
    vR   = vR(indx);
    
    % find the best fitting mixture model (MM)
    [vLog10EtaThresh, MM_Eta, MM_T, MM_R, nComp_best, MM_models] = fit_mm(vT,vR,vEta,NNDpar);
    vEtaThresh = 10.^(vLog10EtaThresh);

    if ~isfield(Model,'EtaThresh')
        Model.EtaThresh = vEtaThresh;
    end
    if strcmp(OptArgs.Display,'on')
        disp(['eta = ',num2str(vEtaThresh)])
        plot_R_T_Eta_distributions(vEta,vT,vR,MM_Eta,MM_T,MM_R,Model);
        plot_RT_distribution(vT,vR,Model);
        plot_Eta_GMM_distributions(vEta,MM_models,Model);
    end

    tic
    if strcmp(OptArgs.NNDmethod,'OriginalNND')
        % perform declustering based on the \eta thresholds
        NND = nnd_decluster_eta(NND,vEtaThresh,NNDpar.dt_max,NNDpar.dr_max);
    elseif strcmp(OptArgs.NNDmethod,'ZBZ2020')
        NND = nnd_decluster_zbz2020(vCat,NND,NNDpar,vEtaThresh,OptArgs.Distance);
    end
    toc
    NND.vLog10EtaThresh = vLog10EtaThresh;
    %save_nnd_values(NND,'nnd_values_eta_threshold.dat');

    if strcmp(OptArgs.Display,'on')
        DiG = nnd_digraph(NND,vCat);
        NNDStats = nnd_cluster_stats(DiG,vCat);
        plot_cluster_stats(NNDStats,vCat,Model);
        plot_nnd_clusters(DiG,NNDStats,vCat,Model);
    end
end

% create a digraph from NND structure
function DiG = nnd_digraph(NND,vCat)
%
    G_attr = nnd_graph_attributes(NND,vCat); % creates graph attributes from the NND structure
    DiG = digraph(G_attr.vEdges(:,1),G_attr.vEdges(:,2),G_attr.vWeights,G_attr.NodeTable,'omitselfloops');
    leafidx = find(DiG.outdegree == 0);  % which nodes are leaves
    DiG.Nodes.Type(leafidx) = {'Leaf'};
    %DiG = rmedge(DiG, 1:numnodes(DiG), 1:numnodes(DiG)); % remove all self-loops
end

%
function NNDStats = nnd_cluster_stats(DiG,vCat)
%
    [nComp_indx, nComp_size] = conncomp(DiG,'Type','weak'); % identify connected components
    NNDStats.nComp_indx = nComp_indx;
    NNDStats.nComp_size = nComp_size;
    NNDStats.ClstProp = nnd_cluster_properties(nComp_indx,nComp_size,DiG,vCat);
    NNDStats.cClstComps = categorical(nComp_size); % create the categorical array based on the cluster component size
    NNDStats.nClst = length(nComp_size); % number of clusters
    NNDStats.nSingle = numel(find(nComp_size == 1)); % number of singleton events    
end

function save_nnd_values(NND,sFileName)
%
    nEq = length(NND.vT);
    fid = fopen(sFileName,'w');
    fprintf(fid,'parent, ind,   Dt,         Ddist,        T,        R,         eta,        parent mag,  mag,        cluster \n');
    %for j = 2:nEq-1 % to compare with usgs code
    for j = 1:nEq
        fprintf(fid,'%4d %4d %13f %11f %11f %11f %11f %11f %11f %4d\n',...
            NND.vP_indx(j),j,NND.vDt(j),NND.vDist(j),NND.vT(j),NND.vR(j),NND.vEta(j),NND.vP_mag(j),NND.vCh_mag(j),NND.vClust(j));
    end
    fclose(fid);
end
