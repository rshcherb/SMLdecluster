function G_attr = nnd_graph_attributes(NND,vCat)
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 3 June 2022
%   ...
%   version 1.2.0, 5 June 2022
%
    nEq = size(NND.vEta,1);
    %nthres = length(vEta_thresh);
    vEdges     = zeros(nEq,2);
    vNodes     = zeros(nEq,1);
    vWeights   = zeros(nEq,1);
    NodeNumber = cell(nEq,1);
    NodeMag    = cell(nEq,1);
    NodeType   = cell(nEq,1);
    for j = 1:nEq
        pind = NND.vP_indx(j); % index of a parent node
        if j ~= pind
            vEdges(j,:) = [pind, j];
            vWeights(j) = radialdist(vCat(pind,2),vCat(pind,3),vCat(j,2),vCat(j,3));
            NodeType{j} = 'Node';
        else
            NodeType{j} = 'Root';
        end
        NodeNumber{j} = num2str(j); % event number
        NodeMag{j} = num2str(sprintf('%.2f',vCat(j,4))); % magnitude
    end
    NodeTable = table(NodeNumber,NodeMag,NodeType,'VariableNames',{'Name' 'Magnitude' 'Type'});
    G_attr.vEdges = vEdges;
    G_attr.vNodes = vNodes;
    G_attr.vWeights = vWeights;
    G_attr.NodeTable = NodeTable;
end

