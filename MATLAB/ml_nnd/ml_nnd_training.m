function Classifier = ml_nnd_training(tTrainCats,tTrainPars,options)
%
%   Authors: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%            Sid Kothari
%
%   version: 1.0.0, 17 October 2024
%   ...
%   version: 1.1.0, 30 October 2024
%
    arguments
        tTrainCats table
        tTrainPars table
        options.Classifier char {matlab.system.mustBeMember(options.Classifier,{'Boosted','Bagged'})} = 'Bagged'
        options.N_NND (1,1) double = 10
        options.dMag logical = 1
        options.FAApar logical = 0
        options.Bkgrd_Int logical = 1
        options.Bkgrd_Prob logical = 1
    end
    % The number of synthetic catalogs with which to train each classifier
    nTrainCats = height(tTrainPars);
    % The length of vTrainers is the number of unique classifiers trained on nTrainers randomly chosen training catalogs
    vTrainers = nTrainCats*ones(1,1);
    % Array of unique training catalog identifier numbers used to acquire a randomly selected subset to train and test
    vTrainCatNums = unique(tTrainCats.Catalog);
    disp('Training classifiers...')
    % Loop to train each classifier with a random selection of 
    for i = length(vTrainers):-1:1
        % Random sample of nTrainCats from vTrainCatNums
        vRandTrainCats = randperm(length(vTrainCatNums),vTrainers(i));
        % Train classifier on the subsampled vRandTrainCats training catalogs
        [Classifier(i),~] = train_classifier(options.Classifier,...
            tTrainCats(ismember(tTrainCats.Catalog,vRandTrainCats),:),...
            tTrainPars(ismember(tTrainPars.Catalog,vRandTrainCats),:),options.N_NND,...
            'Bkgrd_Int',options.Bkgrd_Int,'Bkgrd_Prob',options.Bkgrd_Prob);
    end
end
