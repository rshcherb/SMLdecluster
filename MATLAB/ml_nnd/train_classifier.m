function [trainedClassifier, classificationEnsemble] = train_classifier(Classifier,featCat,ETASpar,N_NND,OptArgs)
%
%   Author: Sid Kothari
%            
%   version 1.0.0, 24 October 2022
%   ...
%   version 1.1.0, 31 October 2024
%
    arguments
        Classifier char ...
            {matlab.system.mustBeMember(Classifier,{'Boosted','Bagged'})}
        featCat table
        ETASpar table
        N_NND (1,1) double
        OptArgs.dMag logical = 1
        OptArgs.FAApar logical = 0
        OptArgs.Bkgrd_Int logical = 1
        OptArgs.Bkgrd_Prob logical = 0
    end

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
inputTable = featCat;
reqVars = repmat(["eta","T","R"],N_NND,1);
reqVars = reqVars(:)'+repmat(1:N_NND,1,3);
% reqVars = ["Lat","Long",reqVars];
% reqVars = repmat(["T","R"],N_NND,1);
% reqVars = reqVars(:)'+repmat(1:N_NND,1,2);
if OptArgs.dMag
    reqVars = [reqVars,"dMag"];
end
if OptArgs.FAApar
    reqVars = [reqVars,"Np","Nc"];
end
if OptArgs.Bkgrd_Int
    if matches("Bkgrd_Int",featCat.Properties.VariableNames)
        reqVars = [reqVars,"Bkgrd_Int"];
    elseif matches("Bkgrd_Intensity",featCat.Properties.VariableNames)
        reqVars = [reqVars,"Bkgrd_Intensity"];
    end
end
if OptArgs.Bkgrd_Prob
    reqVars = [reqVars,"Bkgrd_Prob"];
end
predictorNames = reqVars;
predictors = inputTable(:, predictorNames);
response = inputTable.Type2;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
if strcmpi(Classifier,'Bagged')
    template = templateTree(...
        'MaxNumSplits', 7120);
    classificationEnsemble = fitcensemble(...
        predictors, ...
        response, ...
        'Method', 'Bag', ...
        'NumLearningCycles', 30, ...
        'Learners', template, ...
        'ClassNames', [0; 1]);
    %siteName = extractBetween(featCat.Properties.CustomProperties.folder,...
    %    '.\','\etas_train_cat');
    %siteName = [upper(siteName{1}(1)),siteName{1}(2:end)];
    trainedClassifier.type = 'Bagged Tree ML';
elseif strcmpi(Classifier,'Boosted')
    template = templateTree(...
        'MaxNumSplits', 20);
    classificationEnsemble = fitcensemble(...
        predictors, ...
        response, ...
        'Method', 'AdaBoostM1', ...
        'NumLearningCycles', 30, ...
        'Learners', template, ...
        'LearnRate', 0.1, ...
        'ClassNames', [0; 1]);
    trainedClassifier.type = 'Boosted Tree ML';
end
imp = predictorImportance(classificationEnsemble);
figure()
bar(imp)
ax=gca;
ax.XTick = 1:length(classificationEnsemble.PredictorNames);
ax.XTickLabel = classificationEnsemble.PredictorNames;
ax.XTickLabelRotation = 45;
% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.Features = reqVars;
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
% trainingCatalogs = unique(featCat.CatalogName);
% [~,idx] = sort(str2double(extract(trainingCatalogs,digitsPattern())));
% trainedClassifier.trainingCatalogs = trainingCatalogs(idx);
trainedClassifier.trainingCatalogs = unique(featCat.Catalog);
trainedClassifier.ETASpar = ETASpar;
trainedClassifier.N = size(featCat,1);
trainedClassifier.N_NND = N_NND;

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 5);

% Compute validation predictions
% [validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

fprintf('%s trained on %d synth catalogs - Accuracy %.3f\n',...
    trainedClassifier.type,length(trainedClassifier.trainingCatalogs),...
    100*validationAccuracy)
end
