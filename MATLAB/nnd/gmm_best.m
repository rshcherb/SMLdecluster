function [GMM_best, nComp_best, vThresh, GMM_models] = gmm_best(vData,nmax,OptArgs)
%
%   Fits 1:nmax GMMs to the data and finds the best model with nComp_best components
%   vData - one dimensuonal data
%   nmax - the maximum number of components to fit
%   GMM_best - the best GMM model that fits data returned as gmm object
%   nComp_best - the number of components of the best GMM model
%   vThresh - returns threshold values to separate GMM componets that is used in NND analysis: 
%       OptArgs.ThresholdType = 'intersect' the thresholds are computed as intersection of GMM compenets;
%       OptArgs.ThresholdType = 'minsaddle' the thresholds are computed as the local minima for each saddle if present
%   OptArgs.MinCriterion - specifies wether to use 'AIC' or 'BIC' criterion to find the best model
%
%   Author: Dr. Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version 1.0.0, 27 May 2022
%   ...
%   version 2.0.0, 11 October 2024
%
    arguments
        vData double
        nmax {mustBeInteger}
        OptArgs.MinCriterion char = 'AIC'
        OptArgs.ThresholdType char = 'intersect'
    end

    AIC = zeros(1,nmax);
    BIC = zeros(1,nmax);
    GMM_models = cell(1,nmax);
    options = statset('MaxIter',10000,'TolFun',1e-9);
    for k = 1:nmax
        GMM_models{k} = fitgmdist(vData,k,'Options',options);
        %GMM_models{k} = fitgmdist(vData,k,'RegularizationValue',0.1,'Options',options);
        %GMM_models{k} = fitgmdist(vData,k,'CovarianceType','diagonal','Options',options);
        AIC(k)= GMM_models{k}.AIC;
        BIC(k)= GMM_models{k}.BIC;
    end
    if strcmp(OptArgs.MinCriterion,'AIC')
        [minAIC,nComp_best] = min(AIC);
    elseif strcmp(OptArgs.MinCriterion,'BIC')
        [minBIC,nComp_best] = min(BIC);
    end
    GMM_best = GMM_models{nComp_best};

    vThresh = [];
    if nComp_best > 1 % find thresholds to separate components
        if strcmp(OptArgs.ThresholdType,'intersect')
            vThresh = zeros(1,nComp_best-1);
            [~, indx] = sort(GMM_best.mu); % find indexes of the components sorted according to mean values 
            for k = 1:nComp_best-1
            %for k = indx(1:end-1) % go over components sorted in 
                wk = GMM_best.ComponentProportion(indx(k));
                muk = GMM_best.mu(indx(k));
                sigmak = sqrt(GMM_best.Sigma(indx(k)));
                wk1 = GMM_best.ComponentProportion(indx(k+1));
                muk1 = GMM_best.mu(indx(k+1));
                sigmak1 = sqrt(GMM_best.Sigma(indx(k+1)));
                x0 = GMM_best.mu(indx(k+1));
                %fun = @(x) (wk*normpdf(x,gmmb.mu(k),sqrt(gmmb.Sigma(k))) - wk1*normpdf(x,gmmb.mu(k+1),sqrt(gmmb.Sigma(k+1))));
                wnormpdf = @(x,w,mu,sigma) w*normpdf(x,mu,sigma);
                fun = @(x) ( wnormpdf(x,wk,muk,sigmak) - wnormpdf(x,wk1,muk1,sigmak1) );
                vThresh(k) = fzero(fun,x0);
            end
        elseif strcmp(OptArgs.ThresholdType,'minsaddle')
            gmPDF = @(x) pdf(GMM_best,x);
            x = linspace(min(vData),max(vData),1000)';
            lmindx = islocalmin(gmPDF(x));
            vThresh = x(lmindx);
            % disp(vThresh)
            % figure()
            % plot(x,gmPDF(x),x(lmindx),gmPDF(x(lmindx)),'rs');
        end
    end
end

