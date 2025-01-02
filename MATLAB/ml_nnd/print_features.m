function [sFeaturesOut] = print_features(sFeatures,OptArgs)
%
%   Author: Sid Kothari
%            
%   version 1.0.0, 24 October 2022
%   ...
%   version 1.1.0, 31 October 2024
%
    arguments
        sFeatures string
        OptArgs.delimiter char = ' | '
    end

    N_NeighboursEta = sum(count(sFeatures,"eta"+digitsPattern()));
    N_NeighboursT = sum(count(sFeatures,"T"+digitsPattern()));
    N_NeighboursR = sum(count(sFeatures,"R"+digitsPattern()));
    
    sFeaturesOut = string.empty;
    sFeaturesOut(end+1) = num2str(N_NeighboursT) + " NNDs";
%     if N_NeighboursEta > 1
%         sFeaturesOut(end+1) = "NND_Eta1 - " + num2str(N_NeighboursEta);
%     elseif N_NeighboursEta > 0
%         sFeaturesOut(end+1) = "NND_Eta" + num2str(N_NeighboursEta);
%     end
%     if N_NeighboursT > 1
%         sFeaturesOut(end+1) = "NND_T1 - " + num2str(N_NeighboursT);
%     elseif N_NeighboursT > 0
%         sFeaturesOut(end+1) = "NND_T" + num2str(N_NeighboursT);
%     end
%     if N_NeighboursR > 1
%         sFeaturesOut(end+1) = "NND_R1 - R" + num2str(N_NeighboursR);
%     elseif N_NeighboursR > 0
%         sFeaturesOut(end+1) = "NND_R" + num2str(N_NeighboursR);
%     end

    
    l = N_NeighboursEta + N_NeighboursT + N_NeighboursR + 1;
    
    for i = l:length(sFeatures)
        sFeaturesOut(end+1) = sFeatures(i);
    end
    
    sFeaturesOut = replace(sFeaturesOut,"_","\_");
    sFeaturesOut = join(sFeaturesOut,OptArgs.delimiter);

end