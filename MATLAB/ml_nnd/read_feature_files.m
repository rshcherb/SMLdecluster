function [featCat, ETASpar] = read_feature_files(targetFile,parFile,TimeStart)
%
%   Author: Sid Kothari
%            
%   version 1.0.0, 24 October 2022
%   ...
%   version 1.1.0, 31 October 2024
%
    featCat = table();
    featCat = addprop(featCat,{'name','folder'},{'table','table'});
    featCat.Properties.CustomProperties.folder = targetFile(1).folder;
    for i = 1:length(targetFile)
        NewCat = readtable(fullfile(targetFile(i).folder,targetFile(i).name));
        % RS
        indx = NewCat.Time >= TimeStart;
        NewCat = NewCat(indx,:);
        %
        catNum = str2double(extract(targetFile(i).name,digitsPattern()));
        if isempty(catNum)
            catNum = 1;
        end
        if ismember('Np',NewCat.Properties.VariableNames)
            NewCat.Np = [];
            NewCat.Nc = [];
        end
        if ismember('Bkgrd_Intensity',NewCat.Properties.VariableNames)
            NewCat = renamevars(NewCat,'Bkgrd_Intensity','Bkgrd_Int');
        end
        NewCat.Catalog = repmat(catNum,size(NewCat,1),1);
        NewCat.CatalogName = repmat(string(targetFile(i).name),size(NewCat,1),1);
        featCat = [featCat;NewCat];
        featCat.Properties.CustomProperties.name{i,1} = targetFile(i).name;
        if length(targetFile) < 20
            fprintf('\n%s: %d events\n',targetFile(i).name,size(NewCat,1))
        end
        clear NewCat
    end
    fprintf('\n')
    
    ETASpar = table();
    ETASpar = addprop(ETASpar,{'name','folder'},{'table','table'});
    ETASpar.Properties.CustomProperties.folder = targetFile(1).folder;
    %targetDir = dir(targetFile(1).folder);
    for i = 1:length(targetFile)
        NewPar = readtable(fullfile(parFile(i).folder,parFile(i).name));
        NewPar = NewPar(1,:);
        parNum = str2double(extract(parFile(i).name,digitsPattern()));
        if isempty(parNum)
            parNum = 1;
        end
        NewPar.Catalog = parNum;
        ETASpar = [ETASpar;NewPar];
        ETASpar.Properties.CustomProperties.name{i,1} = parFile(i).name;
        clear NewPar
    end
end

