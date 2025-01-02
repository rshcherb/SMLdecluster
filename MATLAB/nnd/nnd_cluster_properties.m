function ClstProp = nnd_cluster_properties(nComp_indx,nComp_size,DiG,vCat)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 5 June 2022
%   ...
%   version 1.0.0, 5 June 2022
%
    nClst = length(nComp_size);
    ClstProp.Td = zeros(nClst,1);
    ClstProp.dm = zeros(nClst,1);
    ClstProp.avdepth = zeros(nClst,1);
    ClstProp.delta = zeros(nClst,1);
    ClstProp.BI = zeros(nClst,1);
    for nc = 1:nClst
        idx = nComp_indx == nc; % mark nodes belonging to the selected cluster
        if nComp_size(nc) > 1
            vCat_clst = vCat(idx,:);        % select only events belonging to the cluster
            SG_clst = subgraph(DiG,idx);    % subgraph of the full graph containing the cluster
            
            [avdepth, delta, BI] = graph_prop(SG_clst);

            ClstProp.Td(nc) = vCat_clst(end,1) - vCat_clst(1,1); % the duration of the cluster
            ClstProp.dm(nc) = abs(diff(maxk(vCat_clst(:,4),2))); % difference in magnitude between the largest and second largest event
            ClstProp.avdepth(nc) = avdepth;
            ClstProp.delta(nc) = delta;
            ClstProp.BI(nc) = BI;
        end
    end
end

function [avdepth, delta, BI] = graph_prop(G)
%
    N = numnodes(G); % number of nodes
    E = numedges(G); % number of edges
    P = length(find(~strcmp(G.Nodes.Type,'Leaf'))); % number of parent events 

    sources = G.Nodes.Name(strcmp(G.Nodes.Type,'Root'));
    targets = G.Nodes.Name(strcmp(G.Nodes.Type,'Leaf'));
    d = distances(G,sources,targets); % find distances between the root node and all the leaves
    %d = distances(G,sources,targets,'Method','unweighted'); % find distances between the root node and all the leaves

    avdepth = mean(d);
    delta = avdepth/sqrt(N);
    BI = double(P)/double(E);
end
