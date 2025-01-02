function [vLog10EtaThresh, MM_Eta, MM_T, MM_R, nComp_best, MM_models] = fit_mm(vT,vR,vEta,NNDpar,OptArgs)
%
%   Fit a mixture model to the NND distribution for the catalog of events
%   NND - the nearest-neighbour distance structure
%   NNDpar - structure NND parameters
%   OprArgs.MixtureModel - which mixture models to fit: 'GMM', 'Weilbul' 
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 27 October 2024
%   ...
%   version 1.0.0, 27 October 2024
%
    arguments
        vT double
        vR double
        vEta double
        NNDpar struct
        OptArgs.MixtureModel char = 'GMM'   % which mixture model to use: 'GMM', 'WMM' 
        OptArgs.Display char = 'off'        % 'on'/'off' to plot the distributions or not 
    end

    % find the best fitting mixture model (MM)
    if strcmp(OptArgs.MixtureModel,'GMM')
        % fitting is done in log10 space
        MM_T = gmm_best(log10(vT),NNDpar.nGMMmax);
        MM_R = gmm_best(log10(vR),NNDpar.nGMMmax);
        [MM_Eta, nComp_best, vLog10EtaThresh, MM_models] = gmm_best(log10(vEta),NNDpar.nGMMmax,'ThresholdType',NNDpar.ThreshMethod);
        if strcmp(OptArgs.Display,'on')
            disp(['log10(eta) = ',num2str(vLog10EtaThresh)])
        end
    elseif strcmp(OptArgs.MixtureModel,'WMM') % Weibull mixture model
    end
end

