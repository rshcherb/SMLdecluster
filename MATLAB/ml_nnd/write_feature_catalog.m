function write_feature_catalog(sFileName,Model,vCat,vPar,vParErr,options)
%
%   Authors: Robert Shcherbakov, e-mail: rshcherb@uwo.ca
%            Sid Kothari
%
%   version: 1.0.0, 17 October 2024
%   ...
%   version: 1.1.0, 30 October 2024
%
    arguments
        sFileName char
        Model struct
        vCat double
        vPar double
        vParErr double = []
        options.N_NND = 10
        options.dMag logical = 1
        options.FAApar logical = 1
        options.Bkgrd_Int = []
        options.Bkgrd_Prob = []
        options.Distance char = 'Radial'
    end
    
    % Save ETAS parameters
    [vBvalEst, vBvalErr, ~, ~] = gr_fit(vCat(:,4),Model.fMc,Model.fMmax,Model.fDm);
    fid_par = fopen([sFileName,'_feat_par.dat'],'w');
    formSpecNamePar = "%6s      %6s      %6s      %6s      %6s      %6s      %6s      %6s      %6s      %5s\n";
    varNamesPar     = "mu       A        alpha    c        p        d        q        gamma    bval     N";
    formSpecDatPar =  "%6.4f    %6.4f    %6.4f    %6.4f    %6.4f    %6.4f    %6.4f    %6.4f    %6.4f    %5d";
    ParData = [vPar(1),vPar(2),vPar(3),vPar(4),vPar(5),vPar(6),vPar(7),vPar(8),vBvalEst(1),size(vCat,1)];
    if ~isempty(vParErr)
        ParErrData = [vParErr,vBvalErr(1),0];
    else
        ParErrData = [zeros(1,length(vPar)),vBvalErr(1),0];
    end
    fprintf(fid_par,formSpecNamePar,split(varNamesPar).');
    fprintf(fid_par,formSpecDatPar,ParData);
    fprintf(fid_par,"\n"+formSpecDatPar,ParErrData);
    fclose(fid_par);
        
    % Time Lat Lon Mag
    formSpecName(1) = "%12s     %12s     %12s     %3s   ";
    varNames(1)     = "Time     Lat      Lon      Mag";
    formSpecDat(1)  = "%12.3f   %12.3f   %12.3f   %4.2f   "; % Time Lat Lon Mag
    featCat = [vCat(:,1), vCat(:,2), vCat(:,3), vCat(:,4)];
    % NND
    NND = nnd_NN(vCat,Model.NNDpar,'NN',options.N_NND,'Distance',options.Distance);

    reqSpecNames = repmat(["%10s   ","%10s   ","%10s   "],options.N_NND,1);
    formSpecName(end+1) = strjoin(reqSpecNames(:)','');
    reqVars = repmat(["eta","T","R"],options.N_NND,1);
    varNames(end+1) = strjoin(reqVars(:)'+repmat(1:options.N_NND,1,3));
    reqSpecDats = repmat(["%10.6f   ","%10.6f   ","%10.6f   "],options.N_NND,1);
    formSpecDat(end+1) = strjoin(reqSpecDats(:)','');
    featCat = [featCat,log10(NND.vEta),log10(NND.vT),log10(NND.vR)];

    % dMag
    if options.dMag
        formSpecName(end+1) = "%4s   ";
        varNames(end+1)     = "dMag";
        formSpecDat(end+1)  = "%4.1f   ";
        featCat = [featCat,NND.vDmag];
    end
    
    if options.FAApar
        formSpecName(end+1) = "%3s   %3s   ";
        varNames(end+1)     = "Np    Nc";
        formSpecDat(end+1)  = "%3d   %3d";
        [Nc,~] = histcounts(NND.vP_indx(:,1),size(vCat,1));
%         Nc = [Nc]';%;zeros(length(Nc)-floor(bins(end)),1)];
        Np = Nc(NND.vP_indx(:,1));
        featCat = [featCat,Np',Nc'];
    end
        
    % Background intensity
    if ~isempty(options.Bkgrd_Int)
        formSpecName(end+1) = "%12s   " ;
        varNames(end+1)     = "Bkgrd_Int";
        formSpecDat(end+1)  = "%15.6g   ";
        featCat = [featCat,options.Bkgrd_Int];
    end
    
    % Background probabilities
    if ~isempty(options.Bkgrd_Prob)
        formSpecName(end+1) = "%12s   " ;
        varNames(end+1)     = "Bkgrd_Prob";
        formSpecDat(end+1)  = "%12.6f   ";
        featCat = [featCat,options.Bkgrd_Prob];
    end
    
    % Type Type2
    formSpecName(end+1) = "%4s   %5s\n";
    varNames(end+1)     = "Type  Type2";       % Type2: 0 - background, 1 - aftershocks
    formSpecDat(end+1)  = "%4d   %5d\n";
    if size(vCat,2) > 4
        featCat = [featCat,vCat(:,5) - 1,vCat(:,5) ~= 1];
    else
        featCat = [featCat,nan(size(vCat,1),1),nan(size(vCat,1),1)];
    end
    
    % Concatenate feature specs
    formSpecName = strjoin(formSpecName,'');
    varNames = split(strjoin(varNames)).';
    formSpecDat = strjoin(formSpecDat,'');
    
    % Write to file
    fid = fopen([sFileName,'_feat.dat'],'w');
    fprintf(fid,formSpecName,varNames);
    
    for n = options.N_NND+1:size(featCat,1)
        fprintf(fid,formSpecDat,featCat(n,:));
    end
    fclose(fid);
end
