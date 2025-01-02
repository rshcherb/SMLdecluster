function [targetFiles, parFiles] = get_feature_file_names(sDir,sFilePattern)
%
%   Author: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%
%   version: 1.0.0, 6 September 2024
%   ...
%   version: 1.0.0, 6 September 2024
%
    targetDir = dir(fullfile(sDir));
    targetFiles = targetDir(contains({targetDir.name},sFilePattern{1}));
    parFiles = targetDir(contains({targetDir.name},sFilePattern{2}));

    if numel(targetFiles) > 1
        [~,idx] = sort(str2double(extract(extractfield(targetFiles,'name'),digitsPattern())));
        targetFiles = targetFiles(idx);
        parFiles = parFiles(idx);
        fprintf('\nAvailable Catalogs:\n')
        for i = 1:length(targetFiles)
            fprintf('%d. %s\n',i,targetFiles(i).name)
        end
    end
    % parFile = targetDir(contains({targetDir.name},erase(targetFile.name,'.dat')) & ...
    %     endsWith({targetDir.name},'par.dat') & ~contains({targetDir.name},'init'));
end
